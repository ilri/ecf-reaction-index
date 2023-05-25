# This script calculates ECF reaction indices from a standard input CSV file,
# which MUST CONTAIN the following columns, with these exact column names
# (the order in which the columns appear does not matter)
# experimentName
# animalID
# date, dayPostChallenge
# temp
# schizontsLocal, schizontsContra
# piroplasms
# WBC
# exitDay, exitType
#
# optional column: group

library(tidyverse)
library(lubridate) # not automatically loaded by library(tidyverse), and fun to play with dates
library(magrittr)

setwd("~/Work/DRMU/Shiny_stuff/ECF_Index/")

## HARDCODED VARIABLES SECTION ##
input_file <- "historical_data_compiled.csv"
output_csv <- "all_experiments_fully_computed_with_old_index.csv"
## END HARDCODED VARIABLES SECTION ##


read_csv(input_file) -> sampledat

# for later on, we sort the observations by:
# (1) Experiment name
# (2) Animal id
# (3) Day post challenge
# and we remove all experimental days equal to 0
sampledat %<>% filter(dayPostChallenge > 0) %>%
  arrange(experimentName, animalID, dayPostChallenge)
# note the above could well be done upstream,
# but to make sure we don't run into errors, we re-arrange here


## Fundamental recoding:
##
# For piroplasms and schizonts, NA is recoded into 0 (absence of observation)
sampledat %<>% mutate(across(
  .cols = c(paste0("schizonts", c("Local", "Contra")), "piroplasms"),
  .fns = ~replace_na(.x, 0)))


###################################################################
## EXTENDING ALL ANIMALS UNTIL LAST DAY FOR EACH TRIAL  ###########
###################################################################

# Here, we extend observations for animals who quit prematurely the trial.
# Assumption is that the column "exitType" contains either:
# "D" (natural death)
# "E" (euthanasia)
# "R" (removal for treatment)
# BEWARE! Sometimes the last clinical observations are on the exitDay,
# sometimes the day before (e.g. in cases of a death occurring overnight or
# early morning).

# Not to confuse with real observations, we add new columns
# prefixed with "virt" for virtual observations. Starting point:
# we still have a tibble sorted by:
# 1. experimentName,
# 2. animalID,
# 3. dayPostChallenge

# we first add a column experimentDuration that is the max(dayPostChallenge):
sampledat %>% ungroup() %>% mutate(experimentDuration = max(dayPostChallenge),
                                   .by = experimentName, .after = experimentName) -> tmp

# Actually, the above wasn't necessary if
# we use tidyr's complete() function to complete with all test cases, for all
# experiments : in each experiment, we want to get all combinations of
# animalID and dayPostChallenge (and rearrange immediately):
tmp %<>% group_by(experimentName) %>% complete(animalID, dayPostChallenge) %>%
  arrange(experimentName, animalID, dayPostChallenge)

# for rows containing virtual (post-exit) "observations",
# the date column remains NA. We use that to populate these rows.

# ancillary: replicate_last_value() accepts a vector, and gives an output that
# is the same vector where the last NA values have been replaced with the last
# non-NA value
replicate_last_value <- function(x, from = 0) {
  # the "from" argument serves the purpose of indicating the first virtual
  # position, i.e. enables to skip some legit NA observations (e.g. for WBC)

  if(is.na(from) | is.null(from)) return(x) # this is "from" is e.g.
  # the result of some match(TRUE, rep(F, 10))
  
  n <- length(x)

  non_NA_pos <- which(!is.na(x))
  if (length(non_NA_pos) == 0) return(x) # nothing to replace with
  
  last_non_NA_pos <- max(non_NA_pos)
  if(last_non_NA_pos == n) return(x) # the last position of the vector is not NA 
  
  last_value <- x[last_non_NA_pos]
  x[max(last_non_NA_pos+1, from):n] <- last_value # not replacing before "from"
  return(x)
}

# ancillary to the ancillary: get the position of the first NA within the final
# repetition of NAs that extends all the way to the end of the vector.
# in case the last element of the vector is not NA, then we return 0
position_last_stretch_of_NAs <- function(x) {
  length(x) -> n
  if(n == 0 | !is.na(x[n])) return(0)
  
  # the vector has at least one element and his last element is NA.
  # beware that is can still be full of NAs!
  rev(is.na(x)) -> to_walk
  match(FALSE, to_walk) -> first_non_NA_from_end
  if(is.na(first_non_NA_from_end)) return(1) # the input was full of NAs
  
  # else, we return the length of the initial vector - first_non_NA_from_end + 2
  return(n - first_non_NA_from_end + 2)
}

tmp %>% group_by(experimentName, animalID) %>% mutate(across(
  .cols = any_of(c("temp", "schizontsLocal", "schizontsContra", "piroplasms", "WBC",
                   "group", "exitType", "exitDay")),
  .fns = ~replicate_last_value(.x, from = position_last_stretch_of_NAs(date)))) -> sampledat
# BEWARE! Sometimes there could be some interspersed date == NA,
# when clerks did not sample the animals on a given day.
# That's e.g. in experimentName = "Viral vector #3"`, `animalID = "BK269"`.
# But with the way we have written things using position_last_stretch_of_NAs,
# we are safe.
rm(tmp)

###################################################################
## CALCULATING THE VARIABLES NECESSARY FOR THE ECF INDEX ##########
###################################################################

# we will calculate the mean of values greater than 0 with an ancillary function.
# Actually, we won't use this one, since intervening values equal to zero
# for piroplasms and schizonts WILL be integrated in the calculation of the mean.
mean_values_gt_0 <- function(x) { # x is a vector
  indices <- which(!is.na(x) & x > 0)
  mean(x[indices])
}

#########
# IMPORTANT NB: all variables should be calculated in a "so far" mode, i.e. we calculate
# the ECF reaction index for day _n_ based only on data collected up to and including day _n_.
# The only "look-ahead" variables are exitDay and exitType
########

# sanity rules when manipulating several experiments:
# (1) we don't assume any grouping at the beginning of a processing chain.
# (2) all calculations take place after a group_by() on experimentName and animalID.
# (3) whenever we calculate summaries, when we join back to the main table, we join on Experiment AND Animal id.
# (4) to calculate summaries through functions called *within groups*, we use the "normal" summarize(). group_modify() is used to map a function to subtables defined by groups (tibble in, tibble out)

# we use lowercase variables as temporary ones, and uppercase e.g. SC1_FST for the final ones.
# for all experimental days BEFORE the first appearance of schizonts on the local side,
# SC1_FST is set to NA.

# From Nagda & Rowlands (and from their unpublished paper), we get that for pyrexia as well as for schizonts,
# days with absence (i.e. temp below 39.4, or absence of detection of schizonts) contributed a zero value
# into the calculation of the means.
# BUT since their calculation was made as "total over length", those zero values we still included
# as datapoints in the calculation of the mean intensities.
# ACTUALLY, in 2004 they say "values can be negative" for the temperature diffs,
# so we will keep all values from the onset of pyrexia.



# getting the first day the schizont was seen on the "local" side:
# VARIABLE #1: "Days post challenge to first recording of schizont parasitaemia in local drainage
# lymph node."

sampledat %>% filter(schizontsLocal > 0) %>% group_by(experimentName, animalID) %>% summarize(first_local = min(dayPostChallenge), last_local = max(dayPostChallenge), .groups = "drop") %>% left_join(sampledat, ., by = c("experimentName", "animalID")) -> tmp
# NB that an animal for which the schizont is never seen on the local side will simply not appear
# in the result of summarize(), leading (appropriately) to NAs in the join.

#ok, now first_local and last_local will be "corrected" in place to remove any look-ahead information:
tmp %<>% mutate(SC1_FST = ifelse(dayPostChallenge < first_local, NA, first_local), .after = schizontsLocal)
# For last_local, the logic is that if last_local is in the future and first_local is in the past or today, then last_local is modified to be the current experimental day.
# If both first_local and last_local are in the future, then last_local is NA.
tmp %<>% mutate(SC1_LST = ifelse(dayPostChallenge < last_local, ifelse(dayPostChallenge < first_local, NA, dayPostChallenge), last_local), .after = SC1_FST)

# testing here: the below must yield empty tibbles
# tmp %>% filter(SC1_FST > dayPostChallenge) 
# tmp %>% filter(SC1_LST > dayPostChallenge)
# tmp %>% filter(SC1_LST < SC1_FST)

# the below two functions are no longer used
# now that all numerical variables (except temp)
# have seen their NA values converted to 0
cumsum_tweaked <- function(x) {
  # x is a vector. We calculate its cumulative sums. A (partial) vector full of NAs yields NA,
  # but after the first non-NA element, the "na.rm = T" fashion applies.
  # NB: match(FALSE, is.na(x)) yields the first position of a non-NA element (and NA if all elements in x are NA)
  n = length(x) ; result = numeric(n)
  if(n == 0) return(result)
  
  first_pos_non_NA <- match(FALSE, is.na(x))
  if(is.na(first_pos_non_NA)) return(rep(NA_real_, n))
  
  result[1] <- x[1] # can be an NA
  for(i in 2:n) if(i < first_pos_non_NA) result[i] = NA else result[i] = sum(result[i-1], x[i], na.rm = T)
  return(result)
}

cummean_tweaked <- function(x) {
  # x is a vector. We calculate its cumulative (progressive) means. A (partial) vector full of NAs yields NA,
  # but after the first non-NA element, the "na.rm = T" fashion applies.
  # NB: match(FALSE, is.na(x)) yields the first position of a non-NA element (and NA if all elements in x are NA)
  n = length(x) ; result = numeric(n)
  if(n == 0) return(result)
  
  first_pos_non_NA <- match(FALSE, is.na(x))
  if(is.na(first_pos_non_NA)) return(rep(NA_real_, n))
  
  result[1] <- x[1] # can be an NA
  for(i in 2:n) if(i < first_pos_non_NA) result[i] = NA else result[i] = mean(x[1:i], na.rm = T)
  return(result)
}

# now we no longer need cumsum_tweaked since we have zeros,
# and the below is the new version of cummean_tweaked, starting to
# accumulate from the first non_zero value

cummean_tweaked_first_non_zero <- function(x) {
  # x is a vector. We calculate its cumulative (progressive) means
  # starting at the first non-zero element.
  # after the last non-zero element, the cumulative mean remains constant
  # NB: match(FALSE, is.na(x)) yields the first position of a non-NA element (and NA if all elements in x are NA)
  n = length(x) ; result = numeric(n)
  if(n == 0) return(result)
  
  # first occurrence of a non-zero element
  first_pos_non_zero <- match(FALSE, x == 0)
  if(is.na(first_pos_non_zero)) return(rep(0, n))
  
  # last occurrence of a non-zero element
  which(x != 0) -> indices
  last_pos_non_zero <- ifelse(length(indices) == 0, length(x), max(indices))
  
  # building the result vector
  result[1] <- x[1] # can be a zero
  for(i in 2:n) {
    if(i < first_pos_non_zero) result[i] = 0 else
      if (i > last_pos_non_zero) result[i] = result[i-1] else
        result[i] = mean(x[first_pos_non_zero:i], na.rm = T)
  } # endfor
  return(result)
}

# now that we have SC1_LST and SC1_FST, we can calculate SC1_LEN and the mean severity over the period:
tmp %>% mutate(SC1_LEN = SC1_LST - SC1_FST + 1, .after = SC1_LST) %>% group_by(experimentName, animalID) %>% mutate(SC1_AVG = cummean_tweaked_first_non_zero(schizontsLocal), .after = SC1_LEN) -> sampledat
rm(tmp)
# works because we have initially sorted the dataset (1) by animal ID and (2) by experimental day.
# and cummean_tweaked() takes care of interspersed days with missing observations, e.g. row 1705 is missing the observation of schizontsLocal, while there is one for the following day!
# EDIT -> this is no longer relevant, a 0 value having been used:
# sampledat[1705,] %>% View()
# sampledat %>% filter(experimentName == "P67c15350" & animalID == "VLP05") %>% View()


# VARIABLE #2: "Duration of schizont parasitaemia (days)"
# WARNING! Nagda&Rowlands2004 is not explicit whether this is the total number of days of parasitaemia
# irrespective of the side, or whether this is calculated merely on the local (inoculation) side.
# Here, we decide to look only at the local side, as reading Rowlands et al 2000 suggests.
# More precisely, we calculate the distance between the first day and the last day of
# a non-zero schizont observation on the local side, allowing for possible intermediate 0's or NAs.
# (SC1_LEN calculated above already)

# VARIABLE #3: "Average intensity of parasitaemia (from 1 to 3) over period that parasitaemia
# is recorded in local drainage lymph node."
# (SC1_AVG calculated above already)


# getting the number of days between first appearances of schizont on the "local" and on the "contra" side:
# VARIABLE #4 (SC_TRAVEL_TIME): "Days from first recording of schizont parasitaemia in local drainage lymph node
# to first recording of schizont parasitaemia in contralateral, pre-scapular lymph node" 
sampledat %>% filter(schizontsContra > 0) %>% group_by(experimentName, animalID) %>% summarize(first_contra = min(dayPostChallenge), last_contra = max(dayPostChallenge), .groups = "drop") %>% left_join(sampledat, ., by = c("experimentName", "animalID")) -> tmp

# for each animal, this number of days (variable #4) is constant, but is NA until the moment when we have passed
# first_local AND first_contra:
tmp %>% mutate(SC_TRAVEL_TIME = ifelse(dayPostChallenge < first_local | dayPostChallenge < first_contra, NA, first_contra - first_local), .after = schizontsContra) -> sampledat


# VARIABLE #5 (SC2_LEN): duration of schizont parasitaemia in the contralateral side
# will be calculated at the same time as SC2_FST, SC2_LST and
# VARIABLE #6 (SC2_AVG): average intensity of Schizont parasitaemia on the contralateral side.
# We copy-paste and modify below the same stuff we have computed for the local side.

#now, first_contra and last_contra will be "corrected" in place to remove any look-ahead information:
sampledat %>% mutate(SC2_FST = ifelse(dayPostChallenge < first_contra, NA, first_contra), .after = schizontsContra) -> tmp
# For last_local, the logic is that if last_local is in the future and first_local is in the past or today, then last_local is modified to be the current experimental day.
# If both first_local and last_local are in the future, then last_local is NA.
tmp %<>% mutate(SC2_LST = ifelse(dayPostChallenge < last_contra, ifelse(dayPostChallenge < first_contra, NA, dayPostChallenge), last_contra), .after = SC_TRAVEL_TIME)

# now that we have SC2_LST and SC2_FST, we can calculate SC2_LEN and the mean severity over the period on the contra side:
tmp %>% mutate(SC2_LEN = SC2_LST - SC2_FST + 1, .after = SC2_LST) %>% group_by(experimentName, animalID) %>% mutate(SC2_AVG = cummean_tweaked_first_non_zero(schizontsContra), .after = SC2_LEN) -> sampledat
rm(tmp)



# VARIABLE #7 (PIRO_FST_REF_SC1_FST): "Days from first recording of schizont parasitaemia in local drainage lymph node
# to first recording of piroplasms"
sampledat %>% filter(piroplasms > 0) %>% group_by(experimentName, animalID) %>% summarize(first_piro = min(dayPostChallenge), last_piro = max(dayPostChallenge), .groups = "drop") %>% left_join(sampledat, ., by = c("experimentName", "animalID")) -> tmp

# now we "correct" into PIRO_FST_REF_SC1_FST
tmp %<>% mutate(PIRO_FST_REF_SC1_FST = ifelse(dayPostChallenge < first_piro | dayPostChallenge < first_local, NA, first_piro - first_local), .after = piroplasms)
# and now PIRO_LST_REF_SC1_LST
tmp %<>% mutate(PIRO_LST_REF_SC1_FST = ifelse(dayPostChallenge < last_piro, ifelse(dayPostChallenge < first_piro, NA, dayPostChallenge - first_local), last_piro - first_local), .after = PIRO_FST_REF_SC1_FST)

# now we can calculate PIRO_LEN and the log max severity of piroplasms:

logmax_so_far <- function(x) {
  # ultra-simple now that we have zeros instead of NAs: log(0+1) = 0
  # x is a vector. We calculate its log(max(x)) in a step-by-step (cumulative-like) fashion.
  # A (partial) vector full of NAs yields NA,
  # but after the first non-NA element, the "na.rm = T" fashion applies.
  # NB: match(FALSE, is.na(x)) yields the first position of a non-NA element (and NA if all elements in x are NA)
  n = length(x) ; result = numeric(n)
  if(n == 0) return(result)
  
  first_pos_non_zero <- match(FALSE, x ==0)
  if(is.na(first_pos_non_zero)) return(rep(0, n))
  
  result[1] <- log(x[1] + 1) # can be 0, no problem
  for(i in 2:n) if(i < first_pos_non_zero) result[i] = 0 else result[i] = max(result[i-1], log(x[i]+1), na.rm = T)
  return(result)
}

tmp %>% mutate(PIRO_LEN = PIRO_LST_REF_SC1_FST - PIRO_FST_REF_SC1_FST + 1, .after = PIRO_LST_REF_SC1_FST) %>% group_by(experimentName, animalID) %>% mutate(MAX_LG_PIRO = logmax_so_far(piroplasms), .after = PIRO_LEN) -> sampledat
rm(tmp)

# VARIABLE #8 (PIRO_LEN) already calculated above
# VARIABLE #9 (MAX_LG_PIRO) already calculated above

# VARIABLE #10 (FEVER_FST_REF_SC1_FST): "Days from first recording of schizont parasitaemia in local drainage lymph node
# to first recording of pyrexia (rectal temp > 39.4). Value can be negative
sampledat %>% filter(temp > 39.4) %>% group_by(experimentName, animalID) %>% summarize(first_fever = min(dayPostChallenge), last_fever = max(dayPostChallenge), .groups = "drop") %>% left_join(sampledat, ., by = c("experimentName", "animalID")) -> tmp
# correcting
tmp %<>% mutate(FEVER_FST_REF_SC1_FST = ifelse(dayPostChallenge < first_fever | dayPostChallenge < first_local, NA, first_fever - first_local), .after = temp)
# and now the last day with fever:
tmp %<>% mutate(FEVER_LST_REF_SC1_FST = ifelse(dayPostChallenge < last_fever, ifelse(dayPostChallenge < first_fever | dayPostChallenge < first_local, NA, dayPostChallenge - first_local), last_fever - first_local), .after = FEVER_FST_REF_SC1_FST)

# VARIABLE #11 (FEVER_LEN): "Duration of observed pyrexia (days)"
# now we calculate the duration of pyrexia (can include days when the temp
# has dropped below the 39.4°C threshold, then gone up again):

fever_average_so_far <-  function(x) { # NO DIVISION BY 10 HERE!
  # x is a vector. We calculate the average of (x-39.4) in a step-by-step (cumulative-like) fashion.
  # A (partial) vector full of NAs yields NA,
  # but after the first non-NA element, the "na.rm = T" fashion applies.
  # NB: match(FALSE, is.na(x)) yields the first position of a non-NA element (and NA if all elements in x are NA)
  n = length(x) ; fever_yes_no = (x > 39.4) ; result = numeric(n)
  if(n == 0) return(result)
  
  first_fever_position <- match(TRUE, fever_yes_no)
  if(is.na(first_fever_position)) return(rep(NA_real_, n))
  last_fever_position <- max(which(fever_yes_no))
  
  # the below relies heavily in the proper looping order from 1 to n. :D
  for(i in 1:n) result[i] <- ifelse(i < first_fever_position,
                                   NA,
                                   ifelse(i > last_fever_position, result[i-1],
                                          (mean(x[first_fever_position:i], na.rm = T) - 39.4)))
  return(result)
}

tmp %>% mutate(FEVER_LEN = FEVER_LST_REF_SC1_FST - FEVER_FST_REF_SC1_FST + 1, .after = FEVER_LST_REF_SC1_FST) %>% group_by(experimentName, animalID) %>% mutate(FEVER_AVG = fever_average_so_far(temp), .after = FEVER_LEN) -> sampledat

# VARIABLE #12 (FEVER_AVG) already calculated above: "Average level of pyrexia (x°C - 39.4°C) over period
# that pyrexia is recorded.


# VARIABLE #13 (AVG_WBC): "Average white blood cell count between days 13-19 post challenge
# NB: NAs galore in this, because animals can be removed before these days, and WBC measurements are not performed
# every day.
sampledat %<>% mutate(wbc_masked = ifelse(dayPostChallenge < 13 | dayPostChallenge > 19, NA, WBC)) %>% group_by(experimentName, animalID) %>% mutate(AVG_WBC = cummean_tweaked(wbc_masked), .after = WBC)



##############################################
## CALCULATING THE ECF INDEX ITSELF ##########
##############################################

# described in Appendix B of Nagda&Rowlands2004:
ECF_RI_p67_lab_13_vars_23_days <- function(
    SC1_FST, SC1_LEN, SC1_AVG, SC_TRAVEL_TIME,
    SC2_LEN, SC2_AVG,
    PIRO_FST_REF_SC1_FST, PIRO_LEN, MAX_LG_PIRO,
    FEVER_FST_REF_SC1_FST, FEVER_LEN, FEVER_AVG,
    AVG_WBC) {
  
  ECF_RI <- -0.177 * ((SC1_FST - 9) / 1.601)
  ECF_RI <- ECF_RI + 0.298 * (SC1_AVG - 1.613) / 0.613
  ECF_RI <- ECF_RI + 0.322 * (SC1_LEN - 10.53) / 5.087
  ECF_RI <- ECF_RI - 0.251 * (SC_TRAVEL_TIME - 4.5) / 2.132
  ECF_RI <- ECF_RI + 0.309 * (SC2_AVG - 1.03) / 0.713
  ECF_RI <- ECF_RI + 0.337 * (SC2_LEN - 5.379) / 4.901
  ECF_RI <- ECF_RI - 0.199 * (PIRO_FST_REF_SC1_FST - 6.427) / 2.092
  ECF_RI <- ECF_RI + 0.339 * (MAX_LG_PIRO - 2.55) / 1.994
  ECF_RI <- ECF_RI + 0.321 * (PIRO_LEN - 7.142) / 4.428
  ECF_RI <- ECF_RI - 0.042 * (FEVER_FST_REF_SC1_FST - 1.341) / 2.089
  ECF_RI <- ECF_RI + 0.312 * (FEVER_LEN - 9.276) / 4.722
  ECF_RI <- ECF_RI + 0.284 * (FEVER_AVG - 0.685) / 0.341
  ECF_RI <- ECF_RI - 0.260 * (AVG_WBC - 3.731) / 2.219
  
  ECF_RI <- (ECF_RI + 6.3) / 1.13
  
  return(ECF_RI)
}

# same for calculations based on 13 vars and 21 days ("RIGHT" formula):
ECF_RI_p67_lab_13_vars_21_days <- function(
    SC1_FST, SC1_LEN, SC1_AVG, SC_TRAVEL_TIME,
    SC2_LEN, SC2_AVG,
    PIRO_FST_REF_SC1_FST, PIRO_LEN, MAX_LG_PIRO,
    FEVER_FST_REF_SC1_FST, FEVER_LEN, FEVER_AVG,
    AVG_WBC) {
  
  # ad-hoc adjustments first (from Appendix E):
  which(is.na(SC1_FST)) -> i
  SC1_FST[i] <- 14; SC1_LEN[i] <- 0; SC1_AVG[i] <- 0
  
  which(is.na(SC_TRAVEL_TIME)) -> i
  SC_TRAVEL_TIME[i] <- 8; SC2_LEN[i] <- 0; SC2_AVG[i] <- 0
  
  which(is.na(PIRO_FST_REF_SC1_FST)) -> i
  PIRO_FST_REF_SC1_FST[i] <- 10; PIRO_LEN[i] <- 0; MAX_LG_PIRO[i] <- 0
  
  which(is.na(FEVER_FST_REF_SC1_FST)) -> i
  FEVER_FST_REF_SC1_FST[i] <- 6; FEVER_LEN[i] <- 0; FEVER_AVG[i] <- 0

  AVG_WBC[is.na(AVG_WBC)] <- 0 # cancelling the WBC component of the PCA before day 13
  AVG_WBC[AVG_WBC > 8] <- 8 # hard-coded cap
  
  ECF_RI <- -0.177 * ((SC1_FST - 9) / 1.601)
  ECF_RI <- ECF_RI + 0.298 * (SC1_AVG - 1.613) / 0.613
  ECF_RI <- ECF_RI + 0.322 * (SC1_LEN - 9.724) / 4.136 # norm. params. adjusted for 21 days
  ECF_RI <- ECF_RI - 0.251 * (SC_TRAVEL_TIME - 4.5) / 2.132
  ECF_RI <- ECF_RI + 0.309 * (SC2_AVG - 1.03) / 0.713
  ECF_RI <- ECF_RI + 0.337 * (SC2_LEN - 4.871) / 4.101 # norm. params. adjusted for 21 days
  ECF_RI <- ECF_RI - 0.199 * (PIRO_FST_REF_SC1_FST - 6.427) / 2.092
  ECF_RI <- ECF_RI + 0.339 * (MAX_LG_PIRO - 2.55) / 1.994
  ECF_RI <- ECF_RI + 0.321 * (PIRO_LEN - 5.685) / 3.211 # norm. params. adjusted for 21 days
  ECF_RI <- ECF_RI - 0.042 * (FEVER_FST_REF_SC1_FST - 1.341) / 2.089
  ECF_RI <- ECF_RI + 0.312 * (FEVER_LEN - 8.724) / 3.518 # norm. params. adjusted for 21 days
  ECF_RI <- ECF_RI + 0.284 * (FEVER_AVG - 0.685) / 0.341
  ECF_RI <- ECF_RI - 0.260 * (AVG_WBC - 3.731) / 2.219
  
  ECF_RI <- (ECF_RI + 6.3) / 1.13 # we retain the same as Appendix B (p67, lab, 13 vars, 23 days)
  
  # hard boundaries between 0 and 10 (not necessarily scientifically sound!)
  ECF_RI[ECF_RI < 0] <- 0
  ECF_RI[ECF_RI > 10] <- 10
  
  return(ECF_RI)
}

# from Nicholas' docx file for calculations based on 13 vars and 21 days
# BUT with the reference day always being the inoculation date:
OLD_ECF_RI_Nick_lab_13_vars_21_days <- function(
    SC1_FST, SC1_LEN, SC1_AVG,
    SC2_FST, SC2_LEN, SC2_AVG,
    PIRO_FST_REF_SC1_FST, PIRO_LEN, MAX_LG_PIRO,
    FEVER_FST_REF_SC1_FST, FEVER_LEN, FEVER_AVG,
    AVG_WBC) {
  
  # Initial variable transformations, which leave some NAs (for non-lookahead) to be treated later
  # Please note that SC2_FST is already calculated in the input table.
  PIRO_FST <- PIRO_FST_REF_SC1_FST + SC1_FST # reversing the "correct" calculation into the "incorrect" one
  FEVER_FST <- FEVER_FST_REF_SC1_FST + SC1_FST # reversing the "correct" calculation into the "incorrect" one
  
  # ad-hoc adjustments first (from Appendix E):
  which(is.na(SC1_FST)) -> i
  SC1_FST[i] <- 14; SC1_LEN[i] <- 0; SC1_AVG[i] <- 0
  
  which(is.na(SC2_FST)) -> i
  SC2_FST[i] <- SC1_FST[i] + 8; SC2_LEN[i] <- 0; SC2_AVG[i] <- 0 ## MODIFIED!
  
  which(is.na(FEVER_FST)) -> i
  FEVER_FST[i] <- SC1_FST[i] + 6; FEVER_LEN[i] <- 0; FEVER_AVG[i] <- 0
  
  which(is.na(PIRO_FST)) -> i
  PIRO_FST[i] <- SC1_FST[i] + 10; PIRO_LEN[i] <- 0; MAX_LG_PIRO[i] <- 0 ## MODIFIED! and weird, because we can have 24 days here (okay, actually...)
  
  
  AVG_WBC[is.na(AVG_WBC)] <- 0 # not too sure about this: when NA, we zero it?? Okay, it's before the first measurement taken into account (day 13)...
  # let's note, by the way, that in many experiments, there are no WBC values on day 13, but one in day 14...
  AVG_WBC[AVG_WBC > 8] <- 8
  
  ECF_RI <- -0.177 * ((SC1_FST - 9) / 1.601)
  ECF_RI <- ECF_RI + 0.298 * (SC1_AVG - 1.613) / 0.613
  ECF_RI <- ECF_RI + 0.322 * (SC1_LEN - 9.724) / 4.136 # last two numbers altered
  ECF_RI <- ECF_RI - 0.251 * (SC2_FST - 4.5) / 2.132
  ECF_RI <- ECF_RI + 0.309 * (SC2_AVG - 1.03) / 0.713
  ECF_RI <- ECF_RI + 0.337 * (SC2_LEN - 4.871) / 4.101 # last two numbers altered
  ECF_RI <- ECF_RI - 0.199 * (PIRO_FST - 6.427) / 2.092
  ECF_RI <- ECF_RI + 0.339 * (MAX_LG_PIRO - 2.55) / 1.994
  ECF_RI <- ECF_RI + 0.321 * (PIRO_LEN - 5.685) / 3.211 # last two numbers altered
  ECF_RI <- ECF_RI - 0.042 * (FEVER_FST - 1.341) / 2.089
  ECF_RI <- ECF_RI + 0.312 * (FEVER_LEN - 8.724) / 3.518 # last two numbers altered
  ECF_RI <- ECF_RI + 0.284 * (FEVER_AVG - 0.685) / 0.341
  ECF_RI <- ECF_RI - 0.260 * (AVG_WBC - 3.731) / 2.219
  
  ECF_RI <- 0.49 + (ECF_RI + 5.63) / 1.12 # same as Appendix E (schizont 9 vars, 21 days)
  
  # hard boundaries between 0 and 10 (not necessarily scientifically sound!)
  ECF_RI[ECF_RI < 0] <- 0
  ECF_RI[ECF_RI > 10] <- 10
  
  return(ECF_RI)
}

# from the MS Access script Nicholas shared, with parameters for 23 days
# and wrong final standardization:
WRONG_ECF_RI_ACCESSDB_lab_13_vars_23_days <- function(
    SC1_FST, SC1_LEN, SC1_AVG, SC_TRAVEL_TIME,
    SC2_LEN, SC2_AVG,
    PIRO_FST_REF_SC1_FST, PIRO_LEN, MAX_LG_PIRO,
    FEVER_FST_REF_SC1_FST, FEVER_LEN, FEVER_AVG,
    AVG_WBC) {
  
  # ad-hoc adjustments first (from Appendix E):
  which(is.na(SC1_FST)) -> i
  SC1_FST[i] <- 14; SC1_LEN[i] <- 0; SC1_AVG[i] <- 0
  
  which(is.na(SC_TRAVEL_TIME)) -> i
  SC_TRAVEL_TIME[i] <- 8; SC2_LEN[i] <- 0; SC2_AVG[i] <- 0
  
  which(is.na(FEVER_FST_REF_SC1_FST)) -> i
  FEVER_FST_REF_SC1_FST[i] <- 6; FEVER_LEN[i] <- 0; FEVER_AVG[i] <- 0
  
  which(is.na(PIRO_FST_REF_SC1_FST)) -> i
  PIRO_FST_REF_SC1_FST[i] <- 10; PIRO_LEN[i] <- 0; MAX_LG_PIRO[i] <- 0
  
  AVG_WBC[is.na(AVG_WBC)] <- 0
  AVG_WBC[AVG_WBC > 8] <- 8
  
  ECF_RI <- -0.177 * ((SC1_FST - 9) / 1.601)
  ECF_RI <- ECF_RI + 0.298 * (SC1_AVG - 1.613) / 0.613
  ECF_RI <- ECF_RI + 0.322 * (SC1_LEN - 10.53) / 5.087 # last two numbers wrong (23-day)
  ECF_RI <- ECF_RI - 0.251 * (SC_TRAVEL_TIME - 4.5) / 2.132
  ECF_RI <- ECF_RI + 0.309 * (SC2_AVG - 1.03) / 0.713
  ECF_RI <- ECF_RI + 0.337 * (SC2_LEN - 5.379) / 4.901 # last two numbers wrong (23-day)
  ECF_RI <- ECF_RI - 0.199 * (PIRO_FST_REF_SC1_FST - 6.427) / 2.092
  ECF_RI <- ECF_RI + 0.339 * (MAX_LG_PIRO - 2.55) / 1.994
  ECF_RI <- ECF_RI + 0.321 * (PIRO_LEN - 7.142) / 4.428 # last two numbers wrong (23-day)
  ECF_RI <- ECF_RI - 0.042 * (FEVER_FST_REF_SC1_FST - 1.341) / 2.089
  ECF_RI <- ECF_RI + 0.312 * (FEVER_LEN - 9.276) / 4.722 # last two numbers wrong (23-day)
  ECF_RI <- ECF_RI + 0.284 * (FEVER_AVG - 0.685) / 0.341
  ECF_RI <- ECF_RI - 0.260 * (AVG_WBC - 3.731) / 2.219
  
  ECF_RI <- (ECF_RI + 6.5) / 1.15 # wrong, though syntactically close to Appendix B
  
  # hard boundaries between 0 and 10 (not necessarily scientifically sound!)
  ECF_RI[ECF_RI < 0] <- 0
  ECF_RI[ECF_RI > 10] <- 10
  
  return(ECF_RI)
}


sampledat %<>% mutate(
  
  oldCALPR1 = WRONG_ECF_RI_ACCESSDB_lab_13_vars_23_days(SC1_FST, SC1_LEN, SC1_AVG, SC_TRAVEL_TIME, SC2_LEN, SC2_AVG, PIRO_FST_REF_SC1_FST, PIRO_LEN, MAX_LG_PIRO, FEVER_FST_REF_SC1_FST, FEVER_LEN, FEVER_AVG, AVG_WBC),
  
  newRI = ECF_RI_p67_lab_13_vars_21_days(SC1_FST, SC1_LEN, SC1_AVG, SC_TRAVEL_TIME, SC2_LEN, SC2_AVG, PIRO_FST_REF_SC1_FST, PIRO_LEN, MAX_LG_PIRO, FEVER_FST_REF_SC1_FST, FEVER_LEN, FEVER_AVG, AVG_WBC)
  
)


#sampledat %>% filter(!is.na(RI)) %>% View()

write_csv(sampledat, file = output_csv)


