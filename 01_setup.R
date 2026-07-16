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


# ---- 3. Factor levels -----------------------------------------------
# Reference levels set by design: CON and Uninfected are the baselines.
# Treatment order (CON, INV, NAT) is kept identical to 00_merge_clean.R.

total_snails <- total_snails %>%
  mutate(
    Snail_ID         = as.factor(Snail_ID),
    Infection.status = factor(Infection.status, levels = c("Uninfected", "Infected")),
    Treatment        = factor(Treatment,        levels = c("CON", "INV", "NAT")),
    Day.number       = as.factor(Day.number)
  )


# ---- 4. Climbing-height censoring indicator -------------------------
# A snail recorded at exactly 30 cm climbed out of the 30 cm cylinder:
# its true height is >= 30 and could not be observed -> right-censored.
# Zeros are NOT censored and are kept as observed values (all 540 trials
# analysed; see 02_climbing_height.R).
#
# The time-to-emerge censoring columns (censored, time_to_emerge_cens)
# are already built in 00_merge_clean.R and are not rebuilt here.

total_snails <- total_snails %>%
  mutate(h_cens = ifelse(max_height == 30, "right", "none"))


# ---- 5. Checks ------------------------------------------------------

stopifnot(
  nrow(total_snails) == 540,
  n_distinct(total_snails$Snail_ID) == 60,
  all(table(total_snails$Snail_ID) == 9),
  sum(total_snails$emerged == 0) == sum(is.na(total_snails$time_to_emerge)),
  all(total_snails$time_to_emerge <= 20, na.rm = TRUE)
)

cat("h_cens (climbing height):\n");  print(table(total_snails$h_cens))          # right = 92, none = 448
cat("\ncensored (time to emerge):\n"); print(table(total_snails$censored))       # right = 134, none = 406


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


# ---- 9. Helper: repeatability (with CI) from a fitted brms model ----
# R = among-individual variance / total variance, both on the log scale
# for a lognormal model:  R = sd(Snail_ID)^2 / (sd(Snail_ID)^2 + sigma^2)
# Keeps censored trials in, unlike rptR which would silently drop them.

get_R <- function(model, group = "Snail_ID") {

  draws <- as_draws_df(model)

  sd_col <- paste0("sd_", group, "__Intercept")
  stopifnot(sd_col %in% names(draws), "sigma" %in% names(draws))

  var_among <- draws[[sd_col]]^2
  var_resid <- draws[["sigma"]]^2

  R <- var_among / (var_among + var_resid)

  c(R   = median(R),
    lwr = quantile(R, 0.025, names = FALSE),
    upr = quantile(R, 0.975, names = FALSE))
}

