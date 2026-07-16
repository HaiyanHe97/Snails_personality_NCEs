# =====================================================================
# 02  CLIMBING HEIGHT — right-censored at 30 cm
# ---------------------------------------------------------------------
# The cylinder is 30 cm tall. A snail recorded at 30 cm climbed out of
# it: its true climbing height is >= 30 and could not be observed. That
# is a right-censored observation, not a measurement of 30. It happened
# in 92 of 540 trials (17%), disproportionately under the native cue
# where snails climb highest — so treating 30 as a real value would
# systematically compress the treatment effect.
#
# Zeros are NOT censored and are kept as observed values. max_height
# records where the snail was at the END of the trial; a snail below the
# waterline was at 0. Zeros are not subdivided: the data cannot tell us
# WHY a snail didn't climb. All 540 trials are kept.
#
# MAIN:        right-censored Gaussian model (brms).
# SENSITIVITY: old lmer / rptR, 30 treated as a real value (kept for
#              comparison / reviewer-requested sensitivity analysis).
#
# Figures live in 06_figures.R. This script fits the models, prints the
# numerical results, and saves the fitted models to outputs/models/.
#
# WARNING: the brms models take a few minutes each; there are several.
# =====================================================================

source("01_setup.R")

# Ceiling hits by cue (context for the censoring)
table(total_snails$h_cens, total_snails$Treatment)


# =====================================================================
# 1. MAIN MODEL — treatment x infection
# =====================================================================

# --- SENSITIVITY (OLD): no censoring, 30 treated as a real value -----
model_height_old <- lmer(max_height ~ Treatment * Infection.status +
                           (1 | Snail_ID), data = total_snails)
summary(model_height_old)
anova(model_height_old, type = 2)

# --- MAIN: right-censored at 30 -------------------------------------
m_height_int <- brm(
  max_height | cens(h_cens) ~ Treatment * Infection.status + (1 | Snail_ID),
  family  = gaussian(),
  data    = total_snails,
  chains  = CHAINS, cores = 4, iter = ITER, warmup = WARMUP,
  control = list(adapt_delta = 0.95), seed = 123
)

m_height_add <- brm(
  max_height | cens(h_cens) ~ Treatment + Infection.status + (1 | Snail_ID),
  family  = gaussian(),
  data    = total_snails,
  chains  = CHAINS, cores = 4, iter = ITER, warmup = WARMUP,
  control = list(adapt_delta = 0.95), seed = 123
)

summary(m_height_int)
summary(m_height_add)

# Is the interaction worth keeping? (the lmer version said no)
loo_compare(loo(m_height_add), loo(m_height_int))

# Posterior predictive check.
# EXPECT THIS TO LOOK IMPERFECT: 176 of 540 heights are exactly 0, and a
# Gaussian cannot produce a spike at zero, so it predicts some negative
# heights. Censoring fixes the ceiling, not the floor. Effects are large
# and consistent, so conclusions do not hinge on this — note in limitations.
pp_check(m_height_add)

# Pairwise contrasts (main effects, from the additive model)
emm_h_trt <- emmeans(m_height_add, ~ Treatment)
pairs(emm_h_trt)

emm_h_inf <- emmeans(m_height_add, ~ Infection.status)
pairs(emm_h_inf)


# =====================================================================
# 2. REPEATABILITY
# =====================================================================

# ---- 2.1 Overall ----------------------------------------------------

# SENSITIVITY (OLD): rptR, no censoring
rpt_height_old <- rpt(max_height ~ Treatment + Infection.status + (1 | Snail_ID),
                      grname = "Snail_ID", data = total_snails,
                      datatype = "Gaussian", nboot = 1000, npermut = 0)
summary(rpt_height_old)

# MAIN: from the censored model
R_height_overall <- get_R(m_height_add)
print(R_height_overall)

# Emergence is binary — 
rpt_emerged <- rpt(emerged ~ Treatment + Infection.status + (1 | Snail_ID),
                   grname = "Snail_ID", data = total_snails,
                   datatype = "Binary", nboot = 1000, npermut = 0)
summary(rpt_emerged)


# ---- 2.2 By infection status ----------------------------------------

fit_h_inf <- function(status) {
  brm(max_height | cens(h_cens) ~ Treatment + (1 | Snail_ID),
      family = gaussian(),
      data   = subset(total_snails, Infection.status == status),
      chains = CHAINS, cores = 4, iter = ITER, warmup = WARMUP,
      control = list(adapt_delta = 0.95), seed = 123)
}

m_h_inf   <- fit_h_inf("Infected")
m_h_uninf <- fit_h_inf("Uninfected")

R_h_inf   <- get_R(m_h_inf)
R_h_uninf <- get_R(m_h_uninf)

print(rbind(Infected = R_h_inf, Uninfected = R_h_uninf))


# ---- 2.3 By treatment ----------------------------------------------

fit_h_trt <- function(trt) {
  brm(max_height | cens(h_cens) ~ Infection.status + (1 | Snail_ID),
      family = gaussian(),
      data   = subset(total_snails, Treatment == trt),
      chains = CHAINS, cores = 4, iter = ITER, warmup = WARMUP,
      control = list(adapt_delta = 0.95), seed = 123)
}

m_h_CON <- fit_h_trt("CON")
m_h_INV <- fit_h_trt("INV")
m_h_NAT <- fit_h_trt("NAT")

R_h_CON <- get_R(m_h_CON)
R_h_INV <- get_R(m_h_INV)
R_h_NAT <- get_R(m_h_NAT)

print(rbind(CON = R_h_CON, INV = R_h_INV, NAT = R_h_NAT))
# Compare with the uncensored rptR estimates (0.28 / 0.32 / 0.46).
# NAT has the most snails hitting the 30 cm ceiling, so its estimate is
# the one most likely to shift.


# =====================================================================
# 3. INDIVIDUAL EFFECTS (BLUPs, descriptive only)
# =====================================================================
# BLUPs are extracted here for the caterpillar plot in 06_figures.R.

blups <- ranef(m_height_add)$Snail_ID[, , "Intercept"] %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Snail_ID") %>%
  rename(fit = Estimate, lwr = Q2.5, upr = Q97.5)

blups$infection <- total_snails$Infection.status[
  match(blups$Snail_ID, total_snails$Snail_ID)]


# =====================================================================
# 4. CROSS-CONTEXT CONSISTENCY
# ---------------------------------------------------------------------
# Every snail has 3 complete height measurements in every treatment, so
# raw per-snail means are valid. Caveat (one sentence in the paper): the
# 92 ceiling values compress the top of each snail's mean, slightly
# attenuating these correlations. The effect is small and the ranking is
# preserved, so we do not refit.
# =====================================================================

snail_mean <- total_snails %>%
  group_by(Snail_ID, Treatment, Infection.status) %>%
  summarise(mean_height = mean(max_height), .groups = "drop")

snail_wide <- snail_mean %>%
  pivot_wider(names_from = Treatment, values_from = mean_height)

cor.test(snail_wide$CON, snail_wide$NAT)
cor.test(snail_wide$CON, snail_wide$INV)


# =====================================================================
# 5. PLASTICITY — random slopes
# =====================================================================

m_height_slope <- brm(
  max_height | cens(h_cens) ~ Treatment + Infection.status +
    (1 + Treatment | Snail_ID),
  family    = gaussian(),
  data      = total_snails,
  chains    = CHAINS, cores = 4, iter = ITER, warmup = WARMUP,
  control   = list(adapt_delta = 0.99, max_treedepth = 12),
  save_pars = save_pars(all = TRUE),   # needed for moment_match LOO
  seed      = 123
)

summary(m_height_slope)
rstan::check_hmc_diagnostics(m_height_slope$fit)

# Do individuals differ in plasticity?  -> sd(TreatmentINV), sd(TreatmentNAT);
#     a CrI clear of zero means yes.
# Plasticity linked to behavioural type? -> cor(Intercept, Treatment*)
# Respond to one cue => respond to the other? -> cor(TreatmentINV, TreatmentNAT)

loo_slope <- loo(m_height_slope, moment_match = TRUE)
loo_compare(loo(m_height_add), loo_slope)

# ---- Fit by infection status (climbing height, censored) -----------
fit_rn_cov <- function(status) {
  brm(bf(max_height | cens(h_cens) ~ Treatment + (1 + Treatment |q| Snail_ID)),
      data = subset(total_snails, Infection.status == status),
      family = gaussian(),
      chains = CHAINS, cores = 4, iter = ITER, warmup = WARMUP,
      control = list(adapt_delta = 0.99), seed = 123)
}
model_reaction_norm_inf   <- fit_rn_cov("Infected")
model_reaction_norm_uninf <- fit_rn_cov("Uninfected")

summary(model_reaction_norm_inf);   rhat(model_reaction_norm_inf)
summary(model_reaction_norm_uninf); rhat(model_reaction_norm_uninf)


# ---- Slope-slope correlation (INV vs NAT) --------------------------
posterior_INVNAT_inf   <- as_draws_df(model_reaction_norm_inf)$cor_Snail_ID__TreatmentINV__TreatmentNAT
posterior_INVNAT_uninf <- as_draws_df(model_reaction_norm_uninf)$cor_Snail_ID__TreatmentINV__TreatmentNAT

cat("\nslope-slope (INV vs NAT):\n")
cat("  Infected:   r =", round(mean(posterior_INVNAT_inf), 3),
    " 95% [", round(quantile(posterior_INVNAT_inf, 0.025), 3), ", ",
    round(quantile(posterior_INVNAT_inf, 0.975), 3), "]\n", sep = "")
cat("  Uninfected: r =", round(mean(posterior_INVNAT_uninf), 3),
    " 95% [", round(quantile(posterior_INVNAT_uninf, 0.025), 3), ", ",
    round(quantile(posterior_INVNAT_uninf, 0.975), 3), "]\n", sep = "")

diff_INVNAT <- posterior_INVNAT_inf - posterior_INVNAT_uninf
cat("  difference (Infected - Uninfected): ", round(mean(diff_INVNAT), 3),
    " 95% [", round(quantile(diff_INVNAT, 0.025), 3), ", ",
    round(quantile(diff_INVNAT, 0.975), 3), "]\n",
    "  P(infected > uninfected) = ", round(mean(diff_INVNAT > 0), 3), "\n", sep = "")


# =====================================================================
# 6. PREDICTABILITY — double-hierarchical model
# ---------------------------------------------------------------------
# EXPLORATORY: this model was unstable / did not converge cleanly in
# earlier runs (was flagged "not available"). Kept here for completeness
# and possible follow-up; do NOT treat as a confirmed result without
# checking convergence diagnostics.
# =====================================================================

m_height_pred <- brm(
  bf(max_height | cens(h_cens) ~ Treatment + Infection.status + (1 |q| Snail_ID),
     sigma ~ Treatment + Infection.status + (1 |q| Snail_ID)),
  family = gaussian(),
  data   = total_snails,
  prior = c(
    prior(normal(2.3, 0.5), class = "Intercept", dpar = "sigma"),  # log(10) ~ 2.3
    prior(normal(0, 0.5),   class = "b",         dpar = "sigma"),
    prior(exponential(2),   class = "sd",        dpar = "sigma")
  ),
  chains = CHAINS, cores = 4, iter = ITER, warmup = WARMUP,
  control = list(adapt_delta = 0.99), seed = 123
)

summary(m_height_pred)

pred_cor_h <- as_draws_df(m_height_pred)$cor_Snail_ID__Intercept__sigma_Intercept

cat("\nlevel x predictability (climbing height)\n",
    "  r    = ", round(median(pred_cor_h), 3), "\n",
    "  95%  = [", round(quantile(pred_cor_h, 0.025), 3), ", ",
    round(quantile(pred_cor_h, 0.975), 3), "]\n",
    "  P(r>0) = ", round(mean(pred_cor_h > 0), 3), "\n", sep = "")


# =====================================================================
# 7. SAVE FITTED MODELS
# ---------------------------------------------------------------------
# 06_figures.R reads these back to build the plots, so the (slow) models
# are fitted once here rather than refitted for every figure edit.
# outputs/models/*.rds is git-ignored (large, regenerated by this script).
# =====================================================================

saveRDS(
  list(
    m_int   = m_height_int,  m_add   = m_height_add,
    m_slope = m_height_slope, m_pred = m_height_pred,
    m_inf   = m_h_inf,       m_uninf = m_h_uninf,
    m_CON   = m_h_CON,       m_INV   = m_h_INV,   m_NAT = m_h_NAT,
    blups   = blups,
    R_overall = R_height_overall,
    R_inf = R_h_inf,  R_uninf = R_h_uninf,
    R_CON = R_h_CON,  R_INV = R_h_INV, R_NAT = R_h_NAT
  ),
  "outputs/models/height_models.rds"
)
saveRDS(
  list(
    rn_inf = model_reaction_norm_inf, rn_uninf = model_reaction_norm_uninf,
    INVNAT_inf = posterior_INVNAT_inf, INVNAT_uninf = posterior_INVNAT_uninf
  ),
  "outputs/models/covariance_models.rds"
)

