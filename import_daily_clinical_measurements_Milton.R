# read ECF trial data from Milton's input file
setwd("~/Work/DRMU/Shiny_stuff/ECF_Index/Trials_202304/")

library(tidyverse)
library(stringr)
library(magrittr)

### HARDCODED VALUES BELOW **ONLY**

raw <- readxl::read_xlsx("ECF-TLR Experiment 20230418.xlsx")
# because of the intervening "spacer" lines, all variables are characters
#experiment_name <- "Tick Pick-up April 2023"
#output_filename <- "Tick_pick_up_to_20230418.csv"
experiment_name <- "ECF-TLR April 2023"
output_csv <- "TLR_up_to_20230418.csv"

### HARDCODED VALUES ABOVE **ONLY**


# delete rows that are full of NAs
raw %<>% filter(!if_all(.cols = everything(), .fns = is.na))


# where are the spacer lines? They contain "Y [0-9]" in their first cell
# (case-insensitive).
raw[[1]] %>% str_detect(regex("y[[:space:]]?[0-9]", ignore_case = T)) %>% 
  which() -> spacers

# there are as many animals as the second spacer position
# minus the first spacer position - 1
num_animals <- spacers[2] - spacers[1] - 1

# GENIUS idea :) we stop reading the file (i.e. discard its bottom)
# as soon as the distance between too spacers is below half (arbitrary value!!!)
# of num_animals
chunk_sizes <- lead(spacers) - spacers
which.max(is.na(chunk_sizes) | chunk_sizes < num_animals / 2) -> first_wrong_block
spacers[first_wrong_block] -> boundary
total_num_days <- first_wrong_block - 1

# means that there are fewer than (num_animals/2) rows between
# spacers[first_wrong_block] and spacers[first_wrong_block + 1],
# and therefore, we destroy everything beyond (and including) the boundary row
raw %<>% slice(1:(boundary-1))

# HARD STOP: we stop at this stage in case one of the blocks
# doesn't contain the prescribed number of animals
stopifnot("Error: not all the blocks (days) contain the same number\
          of animals!" = chunk_sizes[1:(first_wrong_block-1)] == num_animals + 1)
# TODO: more informative error message telling where is the "missing animal"

# save the string contained in the first spacer
raw[[spacers[1], 1]] -> initial_day_string

# and get which is the first experimental day
initial_day_string %>%
  str_extract(regex("y[[:space:]]*([0-9]*)[[:space:]]*\\(", ignore_case = T), group = 1) %>%
  as.numeric() -> first_day # here we expect 0 or 1

# and what is the first calendar date? It is within parentheses
initial_day_string %>%
  str_extract("\\((.*)\\)", group = 1) %>% lubridate::dmy() -> first_date


# then write the experiment name and experimental day as first column on the left
# (in the following, we still have the spacer lines,
# therefore each block has (num_animals + 1) rows)
raw %<>% mutate(experimentName = experiment_name, dayPostChallenge = rep(first_day:(first_day + total_num_days - 1), each = num_animals + 1), .before = 1)

# And we add dates:
raw %<>% mutate(date = first_date + dayPostChallenge - first_day, .after = experimentName)

# and then remove all spacer lines
raw %<>% slice(-spacers)

# We finally rename the columns, dropping the ones we don't need.
raw %<>% select(1:12) %>% rename(animalID = "ANIMAL ID", temp = TEMP,
                                 schizontsLocal = REG, schizontsContra = LPG,
                                 piroplasms = BLOOD, WBC = Wbc,
                                 RBC = any_of("Rbc"), Hgb = Hg, toto = any_of("zorg"),
                                 `PCV%` = any_of(c("% PCV", "PCV")))
# We have "PCV" in tick pick-up, "% PCV" in TLR

## recoding the variables corresponding to schizonts and to piroplasms:

# recoding the schizont columns:
# we just base our translation on the detection of strings "+++", "++" or "+"
recode_schizont <- function(x) {
  ifelse(grepl("+++", x, fixed = T), 3L,
         ifelse(grepl("++", x, fixed = T), 2L,
                ifelse(grepl("+", x, fixed = T), 1L, 0L)))
}

raw %<>% mutate(across(starts_with("schizont"), recode_schizont))

# recoding the piroplasm column
recode_piro <- function(x) {
  ifelse(is.na(x) | toupper(x) == "NPS", 0L,
         gsub(pattern = "/1000|<", replacement = "", x = x)) -> y
  return(as.integer(y))
}

raw %<>% mutate(across(piroplasms, recode_piro))


# and finally we change the columns to numeric content
# (using any_of because some like RBC may not be here)
raw %<>% mutate(across(any_of(c("temp", "WBC", "RBC", "Hgb", "PCV%")), as.numeric))

# and we remove rows corresponding to dayPostChallenge <= 0
# and we arrange by (increasing) 1. experimentName (though there shall be
# only one here), 2. animalID and 3. dayPostChallenge
raw %<>% filter(dayPostChallenge > 0) %>% 
  arrange(experimentName, animalID, dayPostChallenge)

# and finally, for compatibility with what we will do later, we add the variables
# corresponding to animal exits
raw %<>% mutate(exitDay = NA_integer_, exitType = NA_character_)

write_csv(x = raw, file = output_csv)
