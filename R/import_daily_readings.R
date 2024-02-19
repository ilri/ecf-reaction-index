# read ECF trial data from a spreadsheet file whose input corresponds to
# the data formatting found in xxxxxx.xlsx (sample file to include in /data)
library(tidyverse)
library(stringr)
library(magrittr)
library(rlang)

### RECODING AND PARSING FUNCTIONS ###

# recoding the schizont columns:
# we just base our translation on the detection of strings "+++", "++" or "+"
recode_schizont <- function(x) {
  ifelse(grepl("+++", x, fixed = T), 3L,
         ifelse(grepl("++", x, fixed = T), 2L,
                ifelse(grepl("+", x, fixed = T), 1L, 0L)))
}

# recoding the piroplasm column
recode_piro <- function(x) {
  ifelse(is.na(x) | toupper(x) == "NPS", 0L,
         gsub(pattern = "/1000|<", replacement = "", x = x)) -> y
  return(as.integer(y))
}

# parsing the cell containing the first day (as a string):
# this function yields a list of two elements, the first being the first day
# (a numeric: 0 or 1) and the second the calendar date corresponding to this day.
get_first_day_first_date <- function(s) {
  # Get which is the first experimental day
  s %>%
    str_extract(regex("y[[:space:]]*([0-9]*)[[:space:]]*\\(", ignore_case = T), group = 1) %>%
    as.numeric() -> first_day # here we expect 0 or 1
  
  # And what is the first calendar date? It is within parentheses.
  s %>%
    str_extract("\\((.*)\\)", group = 1) %>% lubridate::dmy() -> first_date
  return(list(first_day, first_date))
}

### MAIN IMPORT FUNCTION ###

# the function below takes as input: 
# (1) an experiment name (character)
# (2) an input tibble corresponding to a raw export of daily measurements over
#     multiple days (experimental readings)
# In case (1) is not set, the function uses the *name* (symbol) of the input tibble (2).
# The function returns the tibble that is ready to be processed by the routine for ECF RI calculations
# Output tibble has columns:
# "experimentName"   "date"             "dayPostChallenge" "animalID"         "temp"             "schizontsLocal"   "schizontsContra"  "piroplasms"
# "WBC"              "RBC"              "Hgb"              "PCV%"             "exitDay"          "exitType"
import_readings <- function(experiment_data, experiment_name = NULL) {
  
  if (is_null(experiment_name)) experiment_name <- as.character(enexpr(experiment_data))

  raw <- experiment_data
  
  # delete rows that are full of NAs
  raw %<>% filter(!if_all(.cols = everything(), .fns = is.na))
  
  
  # where are the spacer lines? They contain "Y [0-9]" in their first cell
  # (case-insensitive).
  raw[[1]] %>% str_detect(regex("y[[:space:]]?[0-9]", ignore_case = T)) %>% 
    which() -> spacers
  
  # prepare to detect chunks full of NA observations, i.e. NA apart from the firs col:
  raw %<>% mutate(all_na = if_all(.cols = -1, .fns = is.na))
  
  # there are as many animals as the second spacer position
  # minus the first spacer position minus 1
  num_animals <- spacers[2] - spacers[1] - 1
  
  # Follow some calculations of chunk sizes and number of null rows in each chunk.
  led_spacers <- lead(spacers)
  chunk_sizes <- led_spacers - spacers - 1
  
  NA_observations_per_chunk <- numeric()
  for (i in seq_along(spacers)) {
    if (is.na(led_spacers[i])) break # actually fixing the iterator's upper limit! ^^
    if (chunk_sizes[i] == 0) { NA_observations_per_chunk[i] <- 0; next } # empty chunk
    start <- spacers[i] + 1
    end <- led_spacers[i] - 1
    NA_observations_per_chunk[i] <- sum(raw$all_na[start:end])
  }
  NA_observations_per_chunk[length(spacers)] <- NA # just to avoid getting a warning below
  
  # Assuming that missing animals are "intentional" (they correspond to animals who have exited the experiment),
  # the first wrong block is the one with no non-null observations
  
  which.max(is.na(chunk_sizes) | chunk_sizes == 0 | NA_observations_per_chunk == chunk_sizes) -> first_wrong_block
  spacers[first_wrong_block] -> boundary
  total_num_days <- first_wrong_block - 1
  
  # We destroy everything beyond (and including) the boundary row
  raw %<>% slice(1:(boundary-1))
  # so we can also trim neat the chunk sizes vector:
  chunk_sizes <- chunk_sizes[1:total_num_days]
  
  # simple warning in case one of the blocks
  # doesn't contain the prescribed number of animals
  if (any(chunk_sizes != num_animals)) {
    warning(paste("Day", which(chunk_sizes != num_animals), "contains fewer animal observations than expected.\n"))
  }
  
  # the string contained in the first spacer
  # contains e.g. "Day 1 (02/02/2024)"
  get_first_day_first_date(raw[[spacers[1], 1]]) -> tmp
  tmp[[1]] -> first_day ; tmp[[2]] -> first_date
  
  # then write the experiment name and experimental day as first column on the left
  # (in the following, we still have the spacer lines,
  # therefore each block has (chunk_sizes[i] + 1) rows)
  raw %<>% mutate(experimentName = experiment_name, dayPostChallenge = rep(first_day:(first_day + total_num_days - 1), times = chunk_sizes + 1), .before = 1)
  
  # And we add dates:
  raw %<>% mutate(date = first_date + dayPostChallenge - first_day, .after = experimentName)
  
  # and then remove all spacer lines
  raw %<>% slice(-spacers) # no warnings even if many spacers are out of bounds now.
  
  # We finally rename the columns, dropping the ones we don't need.
  raw %<>% select(1:12) %>% rename(animalID = "ANIMAL ID", temp = TEMP,
                                   schizontsLocal = REG, schizontsContra = LPG,
                                   piroplasms = BLOOD, WBC = Wbc,
                                   RBC = any_of("Rbc"), Hgb = Hg, toto = any_of("zorg"),
                                   `PCV%` = any_of(c("% PCV", "PCV", "PCV%", "PCV %")))
  # We have "PCV" in tick pick-up, "% PCV" in TLR
  
  ## recoding the variables corresponding to schizonts and to piroplasms:
  
  
  raw %<>% mutate(across(starts_with("schizont"), recode_schizont))
  
  
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
  
  return(raw)
}


### TEST ###
# Warning: the tibble must be a symbol if we want to call enexpr() on it.
# If not, we have to make use a second argument in the below:
# import_readings(readxl::read_xlsx("IACUC2023_26_20240218.xlsx"), "IACUC2023-26") -> tmp

# alternatively:
# totoro <- readxl::read_xlsx("IACUC2023_26_20240218.xlsx")
# import_readings(totoro) -> tmp
