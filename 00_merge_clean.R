# =====================================================================
# 00_merge_clean.R
#
# Builds master_snails.csv: morphology + behaviour, merged and cleaned.
#
# Run this ONCE. All downstream scripts read master_snails.csv and do
# no cleaning of their own.
#
# Paths are relative to the RStudio project root:
#   data/raw/        raw source files (never modified)
#   data/processed/  cleaned output written by this script
# =====================================================================

library(dplyr)
library(stringr)

TRIAL_MAX <- 20   # length of each trial, in minutes


# =====================================================================
# THE SETUP
#
# 30 cm cylinder, filled with 3 cm of cue water. Heights are measured
# from the BOTTOM of the cylinder, so 3 cm = the waterline, 30 cm =
# climbed out of the cylinder, 0 cm = never got above the waterline.
# (This is why no height between 0 and 3 exists in the data.)
#
# 60 snails x 3 treatments x 3 repeats = 540 trials.
# =====================================================================


# =====================================================================
# PARSING THE TIME COLUMN
#
# Times were written as min,sec: "6,45" = 6 min 45 s = 6.75 min.
# Reading the comma as a decimal point gives 6.45 min, which is wrong.
#
# Excel stored the column as a number and dropped trailing zeros, so
# "5,10" shows up as 5.1. The seconds are still recoverable: the
# fractional part x 100 IS the seconds.
#
#   5.02 -> 2 s      5.1 -> 10 s      8.4 -> 40 s
#   2.0  -> 0 s      6.45 -> 45 s
#
# One exception: 3.8 would give 80 s, which is impossible. It is a typo
# for "3,08" = 3 min 8 s, and it is the only such value in the data.
# =====================================================================

parse_minsec <- function(x) {

  x <- str_trim(as.character(x))
  x[str_to_upper(x) %in% c("NA", "")] <- NA_character_   # both "Na" and "NA" occur

  v    <- suppressWarnings(as.numeric(gsub(",", ".", x)))
  mins <- floor(v)
  sec  <- round((v - mins) * 100)

  sec <- ifelse(sec >= 60, sec / 10, sec)   # the single "3,8" typo

  stopifnot(all(sec >= 0 & sec < 60, na.rm = TRUE))
  mins + sec / 60
}


# ---------------------------------------------------------------------
# 1. Read the two source files
# ---------------------------------------------------------------------

morph <- read.csv("data/raw/snail_morphology.csv",           # 60 rows
                  header = TRUE, sep = ";",
                  stringsAsFactors = FALSE, fileEncoding = "UTF-8")

behav <- read.csv("data/raw/05012026_row_data.csv",          # 540 rows
                  header = TRUE, sep = ";",
                  stringsAsFactors = FALSE, fileEncoding = "UTF-8")


# ---------------------------------------------------------------------
# 2. Make the snail IDs match (trim spaces, force uppercase)
# ---------------------------------------------------------------------

morph <- morph %>% mutate(Snail_ID = str_to_upper(str_trim(as.character(Snail_ID))))
behav <- behav %>% mutate(Snail_ID = str_to_upper(str_trim(as.character(Snail_ID))))


# ---------------------------------------------------------------------
# 3. Checks before merging  (all must pass, or the script stops)
# ---------------------------------------------------------------------

stopifnot(
  nrow(morph) == 60,                     # one row per snail
  !any(duplicated(morph$Snail_ID)),
  n_distinct(behav$Snail_ID) == 60,
  length(setdiff(behav$Snail_ID, morph$Snail_ID)) == 0,   # IDs match both ways
  length(setdiff(morph$Snail_ID, behav$Snail_ID)) == 0,
  all(table(behav$Snail_ID) == 9)        # 9 trials per snail
)

# the two files must not disagree about who is infected
status_disagree <- morph %>%
  select(Snail_ID, morph_status = Infection.status) %>%
  inner_join(behav %>% select(Snail_ID, behav_status = Infection.status) %>% distinct(),
             by = "Snail_ID") %>%
  filter(morph_status != behav_status)

stopifnot(nrow(status_disagree) == 0)


# ---------------------------------------------------------------------
# 4. Merge
#    Morphology (1 row per snail) is broadcast across that snail's 9
#    trials. Result: 540 rows.
# ---------------------------------------------------------------------

master <- behav %>%
  inner_join(select(morph, -Infection.status), by = "Snail_ID")

stopifnot(nrow(master) == 540)


# ---------------------------------------------------------------------
# 5. Clean + build variables
# ---------------------------------------------------------------------

master <- master %>%
  rename(
    time_to_emerge = Time.to.3cm.line..min.sec.,   # minutes to first cross the waterline
    max_height     = Max.height..cm.               # height at the end of the trial
  ) %>%
  mutate(

    # min,sec -> decimal minutes (see the note above)
    time_to_emerge = parse_minsec(time_to_emerge),

    # max_height uses the comma as a real decimal mark
    max_height = as.numeric(gsub(",", ".", max_height)),

    # snail G23 had Time.1 entered as "00:00"; it should be 3
    Time.1 = ifelse(Snail_ID == "G23" & Day.number == "3" &
                      Replication == 23 & Time.1 == "00:00",
                    "3", Time.1),

    # reference levels set explicitly, so CON and Uninfected are the
    # baselines by design rather than by alphabetical accident
    Treatment        = factor(Treatment, levels = c("CON", "INV", "NAT")),
    Infection.status = factor(Infection.status, levels = c("Uninfected", "Infected")),
    Snail_ID         = as.factor(Snail_ID),
    Day.number       = as.factor(Day.number),

    # Cue.tank in the raw file was filled down in Excel as T2, T3, T2,
    # T3... ignoring the actual treatment, so it carried no information.
    # Rebuilt from Treatment.
    # (Only one tank per cue, so tank and treatment cannot be separated
    #  statistically -- worth a line in the limitations.)
    Cue.tank = c(CON = "T1", NAT = "T2", INV = "T3")[as.character(Treatment)],

    # ---- derived variables ----

    # NA in time_to_emerge = the snail never crossed the waterline in
    # 20 min. That is an observed outcome, not a missing value, so it
    # becomes emerged = 0.
    emerged = ifelse(is.na(time_to_emerge), 0L, 1L),

    # For survival / time-to-event models. These two columns work as a
    # pair: non-emergers get 20, but `censored = "right"` tells the
    # model their true time is GREATER than 20, not equal to it.
    # Nothing is invented, and no snail is dropped.
    time_to_emerge_cens = ifelse(emerged == 1, time_to_emerge, TRIAL_MAX),
    censored            = ifelse(emerged == 1, "none", "right"),
    # which of the 3 repeats, centred (-1, 0, 1) so the model intercept
    # refers to the middle repeat
    rep_c = as.numeric(as.character(Time.1)) - 2
  ) %>%
  select(-any_of(c("Notes", "Comment")))   # free-text columns, not needed


# ---------------------------------------------------------------------
# 6. Checks after merging  (all must pass)
# ---------------------------------------------------------------------

stopifnot(
  nrow(master) == 540,
  all(table(master$Snail_ID) == 9),
  all(table(master$Snail_ID, master$Treatment) == 3),   # 3 repeats per cue
  all(master$time_to_emerge >  0,         na.rm = TRUE),
  all(master$time_to_emerge <= TRIAL_MAX, na.rm = TRUE),
  all(master$max_height >= 0, na.rm = TRUE),
  all(master$max_height == 0 | master$max_height >= 3, na.rm = TRUE),  # 0, or above the waterline
  sum(master$emerged == 0) == sum(is.na(master$time_to_emerge)),
  all(master$rep_c %in% c(-1, 0, 1))
)


# ---------------------------------------------------------------------
# 7. Confirm the time parsing worked
# ---------------------------------------------------------------------

# Under the old parsing, seconds of 10/20/30/40/50 were impossible.
# They should all be present now.
secs <- round((master$time_to_emerge %% 1) * 60)
cat("seconds present:\n"); print(sort(unique(secs[!is.na(secs)])))
cat("\n10/20/30/40/50 present:", all(c(10, 20, 30, 40, 50) %in% secs), "\n")
cat("all values on whole seconds:",
    all(abs(master$time_to_emerge * 60 - round(master$time_to_emerge * 60)) < 1e-6,
        na.rm = TRUE), "\n\n")

summary(master[, c("time_to_emerge", "max_height", "emerged")])
table(master$Treatment, master$emerged, master$Infection.status)


# ---------------------------------------------------------------------
# 8. Save
# ---------------------------------------------------------------------

write.csv(master, "data/processed/master_snails.csv",
          row.names = FALSE, fileEncoding = "UTF-8")

View(master)
# =====================================================================
# COLUMNS ADDED BY THIS SCRIPT
#
#   time_to_emerge        minutes to first cross the waterline.
#                         NA = never crossed it (134 trials).
#   emerged               1 = crossed the waterline at least once
#                         0 = never did
#   time_to_emerge_cens   as above, but non-emergers get 20
#   censored              "none" = time observed
#                         "right" = true time is > 20 min
#   rep_c                 repeat number, centred: -1, 0, 1
#   Cue.tank              rebuilt from Treatment
#
#
# ONE THING TO WATCH: max_height == 0 means two different things.
#   emerged == 0  ->  never reached the waterline        (134 trials)
#   emerged == 1  ->  reached it, ended the trial below   (42 trials)
#
# DECISION (agreed with supervisor): climbing height is analysed on
# all 540 trials. Both kinds of zero are kept as observed zeros; we do
# not split them. Heights reaching the 30 cm cylinder top are
# right-censored at 30 (see 01_setup.R / 02_climbing_height.R).
# =====================================================================
