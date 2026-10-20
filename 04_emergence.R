# =====================================================================
# 04  EMERGENCE (yes / no)
# ---------------------------------------------------------------------
# Whether a snail crossed the waterline at least once in the 20-min
# trial (emerged = 1) or never did (emerged = 0).
#
# Binomial GLMM. This is the ANCHOR analysis: it does not depend on the
# censoring assumptions used for climbing height (02) and time to emerge
# (03), so it provides a check on the same behavioural signal using a
# model with no censoring involved.
#
# Together vs separate: same rule as in 02. Here only ONE-model analyses
# (means, overall repeatability); no group comparison of repeatability
# (see section 2).
#
# Figures live in 06_figures.R.
# =====================================================================

source("01_setup.R")


# =====================================================================
# 1. MAIN MODEL — Treatment x Infection   [ONE model: means]
# =====================================================================

model_emerged <- glmer(
  emerged ~ Treatment * Infection.status + (1 | Snail_ID),
  data   = total_snails,
  family = binomial
)

summary(model_emerged)
Anova(model_emerged, type = 2)

# Pairwise contrasts = MARGINAL means on the logit scale (emmeans averages
# over the other factor under the interaction model). Add
# type = "response" to pairs() to get odds ratios instead.
emmeans_treatment <- emmeans(model_emerged, ~ Treatment)
pairs(emmeans_treatment)

emmeans_infection <- emmeans(model_emerged, ~ Infection.status)
pairs(emmeans_infection)


# =====================================================================
# 2. REPEATABILITY (overall only)   [ONE model: do snails differ
#    consistently?]
# ---------------------------------------------------------------------
# Emergence is binary and needs no censoring, so rptR is used directly.
# NOT subdivided by treatment or infection: those subsets (n ~ 30-135)
# frequently gave singular fits / estimates indistinguishable from zero,
# reflecting the low information content of a binary measure at these
# sample sizes. Only the overall repeatability is reported.
# rptR gives R on the latent (link) and the original scale; state in the
# paper which one is reported (usually the link scale).
# =====================================================================

rpt_emerged <- rpt(
  emerged ~ Treatment + Infection.status + (1 | Snail_ID),
  grname   = "Snail_ID",
  data     = total_snails,
  datatype = "Binary",
  nboot    = 1000,
  npermut  = 0
)
summary(rpt_emerged)
