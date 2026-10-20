# =====================================================================
# 01_setup.R
#
# Loads the cleaned data, sets factor levels, builds the climbing-height
# censoring indicator, defines the colour palette and shared helpers.
#
# Every analysis script (02-06) begins with:  source("01_setup.R")
# so the data and helpers are always prepared the same way.
#
# Paths are relative to the RStudio project root.
# =====================================================================

# ---- 1. Libraries ---------------------------------------------------
library(dplyr)
library(tidyr)
library(ggplot2)
library(lme4)
library(lmerTest)
library(emmeans)
library(car)
library(rptR)
library(patchwork)
library(psych)
library(brms)


# ---- 2. Load cleaned data (built by 00_merge_clean.R) ---------------

total_snails <- read.csv("data/processed/master_snails.csv",
                         header = TRUE, sep = ",",
                         stringsAsFactors = FALSE, fileEncoding = "UTF-8")

names(total_snails)


# ---- 3. Factor levels -----------------------------------------------
# Reference levels set by design: CON and Uninfected are the baselines.
# Treatment order (CON, INV, NAT) is kept identical to 00_merge_clean.R.



total_snails <- total_snails %>%
  mutate(
    across(c(Weight_g, Height_cm, Aperture_cm, Width_cm),
           ~ as.numeric(gsub(",", ".", as.character(.x)))),
    Snail_ID         = as.factor(Snail_ID),
    Infection.status = factor(Infection.status, levels = c("Uninfected", "Infected")),
    Treatment        = factor(Treatment,        levels = c("CON", "INV", "NAT")),
    Day.number       = as.factor(Day.number)
  )

str(total_snails[, c("Weight_g", "Height_cm", "Aperture_cm", "Width_cm")])
# ---- 4. Climbing height: height above the waterline -----------------
# The cylinder was filled to the 3 cm waterline. max_height was recorded
# from the base, so it includes the 3 cm water column.
#
# Conversion to height ABOVE the waterline:
#   max_height == 0  -> snail in the water (never emerged, or emerged
#                       and returned to the water) -> 0
#   max_height == 3  -> just emerged past the waterline, no further
#                       climbing. Recorded as 3.1 (0.1 cm above water)
#                       so these are DISTINCT from the in-water zeros.
#   max_height  > 3  -> emerged; height above water = max_height - 3
#   max_height == 30 -> climbed out of the cylinder (ceiling)
#                       -> 27 cm above water, RIGHT-CENSORED
#
# Censoring indicator uses the ORIGINAL ceiling hit (max_height == 30).

total_snails <- total_snails %>%
  mutate(
    max_height_aw = case_when(
      max_height == 0  ~ 0,          # stayed in water
      max_height == 3  ~ 0.1,        # just above the waterline, no further climb
      max_height == 30 ~ 27,         # ceiling -> censoring point (30 - 3)
      TRUE             ~ max_height - 3
    ),
    h_cens = ifelse(max_height == 30, "right", "none")
  )

# sanity checks
stopifnot(
  all(total_snails$max_height_aw >= 0),
  all(total_snails$max_height_aw <= 27),
  sum(total_snails$max_height == 3)  == sum(total_snails$max_height_aw == 0.1),
  sum(total_snails$max_height == 0)  == sum(total_snails$max_height_aw == 0),
  sum(total_snails$h_cens == "right") == sum(total_snails$max_height == 30)
)

# Height in dm, used only by the individual stability model in 02 (on the
# cm scale that model cannot initialise; see 02, section 6).
total_snails$h_dm <- total_snails$max_height_aw / 10

# ---- 5. Checks ------------------------------------------------------

stopifnot(
  nrow(total_snails) == 540,
  n_distinct(total_snails$Snail_ID) == 60,
  all(table(total_snails$Snail_ID) == 9),
  sum(total_snails$emerged == 0) == sum(is.na(total_snails$time_to_emerge)),
  all(total_snails$time_to_emerge <= 20, na.rm = TRUE),
  all(sapply(total_snails[, c("Weight_g", "Height_cm", "Aperture_cm", "Width_cm")], is.numeric)),
  !any(is.na(total_snails[, c("Weight_g", "Height_cm", "Aperture_cm", "Width_cm")]))
)
# ---- Censoring by treatment -----------------------------------------

# time to emerge: right-censored = did not emerge within 20 min
tab_lat <- with(total_snails, table(Treatment, censored))

# climbing height: right-censored = reached the 30 cm ceiling
tab_hgt <- with(total_snails, table(Treatment, h_cens))

cens_by_treatment <- data.frame(
  Treatment      = levels(total_snails$Treatment),
  n              = as.vector(rowSums(tab_lat)),
  latency_cens   = as.vector(tab_lat[, "right"]),
  latency_pct    = round(100 * tab_lat[, "right"] / rowSums(tab_lat), 1),
  height_cens    = as.vector(tab_hgt[, "right"]),
  height_pct     = round(100 * tab_hgt[, "right"] / rowSums(tab_hgt), 1),
  row.names = NULL
)

print(cens_by_treatment)

all_cens <- total_snails %>%
  group_by(Treatment, Snail_ID) %>%
  summarise(
    lat_all_cens = all(censored == "right"),
    hgt_all_cens = all(h_cens   == "right"),
    .groups = "drop"
  ) %>%
  group_by(Treatment) %>%
  summarise(
    snails_lat_all_cens = sum(lat_all_cens),
    snails_hgt_all_cens = sum(hgt_all_cens),
    .groups = "drop"
  )

print(all_cens)
# Check: latency_cens should sum to 134, height_cens to 92
# ---- 6. Colour palette (consistent across all plots) ----------------

col_infected   <- "#AC2527"  # red
col_uninfected <- "#2B4868"  # dark blue

col_con <- "#61221A"  # deep maroon
col_inv <- "#F8CC5A"  # yellow
col_nat <- "#5C8447"  # green

col_infected_light   <- "#D99A9B"
col_uninfected_light <- "#8FA3B3"


# ---- 7. Output directory for saved figures --------------------------

fig_dir <- "outputs/figures"


# ---- 8. brms sampler settings (shared by all censored models) -------

CHAINS <- 4; ITER <- 4000; WARMUP <- 1000


# ---- 9. Helpers: repeatability and posterior summaries --------------
# R = among-individual variance / total variance, both on the log scale
# for a lognormal model:  R = sd(Snail_ID)^2 / (sd(Snail_ID)^2 + sigma^2)
# Keeps censored trials in, unlike rptR which would silently drop them.
# get_R_draws() returns the posterior draws of R (needed for differences
# between groups); get_R() returns median and 95% CrI, as before.

get_R_draws <- function(model, group = "Snail_ID") {
  draws  <- as_draws_df(model)
  sd_col <- paste0("sd_", group, "__Intercept")
  stopifnot(sd_col %in% names(draws), "sigma" %in% names(draws))
  var_among <- draws[[sd_col]]^2
  var_among / (var_among + draws[["sigma"]]^2)
}

get_R <- function(model, group = "Snail_ID") {
  R <- get_R_draws(model, group)
  c(R   = median(R),
    lwr = quantile(R, 0.025, names = FALSE),
    upr = quantile(R, 0.975, names = FALSE))
}

# Median and 95% CrI of any posterior; for a difference between groups
# also P(difference > 0). Groups are compared by the posterior of the
# difference, not by whether their CrIs overlap.
summ_draws <- function(x) round(c(median = median(x),
                                  lwr = quantile(x, 0.025, names = FALSE),
                                  upr = quantile(x, 0.975, names = FALSE)), 3)
summ_diff  <- function(x) c(summ_draws(x), P_gt0 = round(mean(x > 0), 3))

