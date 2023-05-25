
########################################################################
## The goal of this script is to import historical clinical data      ##
## for ECF index calculations, as per the format given by Nick Ndiwa. ##
## The result is a neatly formatted CSV file.                         ##
########################################################################

library(tidyverse)
library(lubridate) # not automatically loaded by library(tidyverse), and fun to play with dates
library(magrittr)
library(stringr)

## HARDCODED SECTION
setwd("~/Work/DRMU/Shiny_stuff/ECF_Index/")
input_file <- "Original_Datasets/All experiments raw data.xlsx"
output_csv <- "historical_data_compiled.csv"
## END HARDCODED SECTION


# The DAILYDATA sheet contains most of the stuff (i.e. the clinical observations)
readxl::read_xlsx(input_file, sheet = "DAILYDATA") %>% select(-c(Treatment, Remarks)) %>% rename(experimentName = Expt_nane, date = Date, schizontsLocal = Schizont_local, schizontsContra = Schizont_contra, piroplasms = Piroplasms, animalID = Anim_2, temp = Rectal_temp, Hgb = HGB, `PCV%` = PCV) -> sampledat
# variable renaming just so that we can work with standard names

unique(sampledat$experimentName) # we need to fix some experiment names
sampledat[sampledat$experimentName == "p67c_nanopartide", "experimentName"] <- "p67C nanoparticle"
sampledat[sampledat$experimentName == "Viralvector", "experimentName"] <- "Viral vector"

# We also import information regarding the group and the death or euthanasia of the animals:
# Warning! Column for groups and animal IDs were swapped for P67c15350:
# I fixed that before importing.
readxl::read_xlsx(input_file, sheet = "BASIC DATA") %>% select(experimentName = Expt_name, animalID = Anim_2, group = Group, exitType = PM_examin, exitDate = Death_DT) -> metadata

metadata[metadata$experimentName == "p67c_nanopartide", "experimentName"] <- "p67C nanoparticle"
metadata[metadata$experimentName == "Viralvector", "experimentName"] <- "Viral vector"

metadata %>% left_join(sampledat, ., by = c("experimentName", "animalID")) -> sampledat
rm(metadata)

### separating the three distinct "Viral vector" experiments
### Viral vector #1: 11/11/2014 to 25/11/2014, animal IDs VV00[1-5]
### Viral vector #2: 02/12/2014 to 16/12/2014, animal IDs [1-5]
### Viral vector #3: 28/04/2015 to 11/05/2015, animal IDs BKxxx
sampledat[sampledat$experimentName == "Viral vector" &
            str_detect(sampledat$animalID, "^VV00") , "experimentName"] <- "Viral vector #1"

sampledat[sampledat$experimentName == "Viral vector" &
            str_detect(sampledat$animalID, "^BK") , "experimentName"] <- "Viral vector #3"

# then the remaining ones:
sampledat[sampledat$experimentName == "Viral vector", "experimentName"] <- "Viral vector #2"

#### CHECKING FOR DATA INTEGRITY ISSUES ####


## data integrity: temperatures

boxplot(sampledat$temp) # some hypothermia-like points...
hist(sampledat$temp)
sampledat %>% filter(temp <= 36) # Lowest temp is 27.8°C for VLP04 in the P67c15350 experiment. All other above 30°C.

# let's see that VLP04 datapoint
# sampledat %>% filter(experimentName == "P67c15350" & animalID == "VLP04") %>% arrange(date) %>% View() # obviously, this "27.8" was actually a "37.8". We fix.
sampledat[!is.na(sampledat$temp) & sampledat$temp == 27.8, "temp"] <- 37.8
# the above after having double-checked there was only one such observation.

# second lowest temp is in MinimumDose, 30°C for animal LA9:
# sampledat %>% filter(experimentName == "MinimumDose" & animalID == "LA9") %>% arrange(date) %>% View() # there is a mistake here. Writing 38.45 as the mean of the previous and
# the following observations.
sampledat[!is.na(sampledat$temp) & sampledat$temp == 30.0, "temp"] <- 38.45
# the above after having double-checked there was only one such observation.

# next is a "30.9" temp, in MinimumDose again:
# sampledat %>% filter(experimentName == "MinimumDose" & animalID == "LA17") %>% arrange(date) %>% View() # looks this 30.9 is a 40.9, with high piroplasms (12):
sampledat[!is.na(sampledat$temp) & sampledat$temp == 30.9, "temp"] <- 40.9
# the above after having double-checked there was only one such observation.

# next again: values of 35.5...
# sampledat %>% filter(experimentName == "combinedExper" & animalID == "LS35") %>% arrange(date) %>% View()
# rewriting into 38.5:
sampledat[!is.na(sampledat$temp) & sampledat$temp == 35.5 &
            sampledat$animalID == "LS35", "temp"] <- 38.5
# the above after having double-checked there was only one such observation.

# sampledat %>% filter(experimentName == "P67c15350" & animalID == "VLP07") %>% arrange(date) %>% View() # this last one is actually also a 38.5:
sampledat[!is.na(sampledat$temp) & sampledat$temp == 35.5 , "temp"] <- 38.5
# we double-check there are no remaining outliers remaining in the boxplot of temps


## data integrity: schizonts
# sampledat %>% filter(schizontsLocal < 0 | schizontsLocal > 3) %>% View()
# we see a "-3" in "MinimumDose", otherwise a bunch of +4 in a "Viral vector".
# we will fix the "-3" into "3" and leave the +4s
# sampledat %>% filter(experimentName == "MinimumDose" & animalID == "LA22") %>% arrange(date) %>% View()
sampledat[!is.na(sampledat$schizontsLocal) & sampledat$schizontsLocal == -3,
          "schizontsLocal"] <- 3

# sampledat %>% filter(schizontsContra < 0 | schizontsContra > 3) %>% View()
# here also, a number of "-3" we transform into "3".
sampledat[!is.na(sampledat$schizontsContra) & sampledat$schizontsContra == -3,
          "schizontsContra"] <- 3

## data integrity: piroplasms
range(sampledat$piroplasms, na.rm = T) # ok, from 1 to 340


#########################################################
## CREATING DAY OF EXPERIMENT AS A NUMERICAL VARIABLE ###
#########################################################

# The below function takes as argument such a tibble as above, an optional column name (default is "Date")
# and an optional start date as a "YYYY-MM-DD" string.
# In case a date is given, that date is taken as the date of inoculation (day 0). If no date is given, the earliest date
# in the tibble is taken as Day 1

# note that the date column should be a "POSIXct" one
# and that a date (e.g. the result of a as_date()) is a day stored as the number of days since 1970-01-01
# At the same time, we are mutating the exitDate into an exitDay
define_start_date <- function (x, date_colname = "date", inoc_date) {
  if (base::missing(inoc_date)) inoc_date = min(`[[`(x, date_colname)) - 1 # defining the inoculation date if missing. It is a POSIXct object.
  inoc_date %<>%  as_date()
  x %>% mutate(dayPostChallenge = as.numeric(as_date(get(date_colname)) - inoc_date), .after = {{date_colname}}) %>% mutate(exitDay = as.numeric(as_date(exitDate) - inoc_date, .after = exitDate))
}

# final alterations:
# add dayPostChallenge, create exitDay, remove exitDate, keep only the
# date component for the column "date"
sampledat %<>% group_by(experimentName) %>%
  group_modify(~ define_start_date(.x)) %>%
  mutate(date = date(date)) %>% select(-exitDate) %>% 
  arrange(experimentName, animalID, dayPostChallenge) 
# very important to arrange as above, because we then have some cumulative sums depending on this ordering


write_csv(sampledat, file = output_csv)
