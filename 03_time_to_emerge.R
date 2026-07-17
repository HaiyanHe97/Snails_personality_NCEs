# =====================================================================
# 03  TIME TO EMERGE — right-censored at 20 min (all 540 trials)
# ---------------------------------------------------------------------
# rptR cannot handle censored data: feeding it time_to_emerge would
# silently drop the 134 non-emergers, and those are not random (98
# infected, 36 uninfected). So everything is fitted with censored brms
# models (lognormal), which keep all 540 trials — a censored trial
# contributes the information "this snail took longer than 20 min"
# rather than being discarded.
#
# Repeatability is read from the posterior with get_R() (defined in 01),
# using the same definition rptR uses, on the log scale:
#     R = sd(Snail_ID)^2 / (sd(Snail_ID)^2 + sigma^2)
#
# Figures live in 06_figures.R. This script fits the models, prints the
# numerical results, and saves the fitted models to outputs/models/.
#
# WARNING: several brms models, a few minutes each.
# =====================================================================

source("01_setup.R")


# =====================================================================
# 1. MAIN MODEL — treatment + infection (and interaction check)
# =====================================================================

m_tte <- brm(
  time_to_emerge_cens | cens(censored) ~ Treatment + Infection.status +
    (1 | Snail_ID),
  family  = lognormal(),
  data    = total_snails,
  chains  = CHAINS, cores = 4, iter = ITER, warmup = WARMUP,
  control = list(adapt_delta = 0.95), seed = 1
)

summary(m_tte)
pp_check(m_tte)
#pochoc

emm_tte_trt <- emmeans(m_tte, ~ Treatment)
pairs(emm_tte_trt)
emm_tte_inf <- emmeans(m_tte, ~ Infection.status)
pairs(emm_tte_inf)


# Overall repeatability
rpt_tte_overall <- get_R(m_tte)
print(rpt_tte_overall)

# Does the treatment x infection interaction improve the model?
m_tte_int <- brm(
  time_to_emerge_cens | cens(censored) ~ Treatment * Infection.status +
    (1 | Snail_ID),
  family = lognormal(), data = total_snails,
  chains = CHAINS, cores = 4, iter = ITER, warmup = WARMUP,
  control = list(adapt_delta = 0.95), seed = 1
)

loo_compare(loo(m_tte), loo(m_tte_int))   # interaction supported?
fixef(m_tte_int)
hypothesis(m_tte_int, "TreatmentINV:Infection.statusInfected < 0")


# =====================================================================
# 2. REPEATABILITY BY INFECTION STATUS  (treatments pooled)
# =====================================================================

m_tte_inf <- brm(
  time_to_emerge_cens | cens(censored) ~ Treatment + (1 | Snail_ID),
  family = lognormal(),
  data   = subset(total_snails, Infection.status == "Infected"),
  chains = CHAINS, cores = 4, iter = ITER, warmup = WARMUP,
  control = list(adapt_delta = 0.95), seed = 1
)

m_tte_uninf <- brm(
  time_to_emerge_cens | cens(censored) ~ Treatment + (1 | Snail_ID),
  family = lognormal(),
  data   = subset(total_snails, Infection.status == "Uninfected"),
  chains = CHAINS, cores = 4, iter = ITER, warmup = WARMUP,
  control = list(adapt_delta = 0.95), seed = 1
)

rpt_tte_infected   <- get_R(m_tte_inf)
rpt_tte_uninfected <- get_R(m_tte_uninf)

print(rpt_tte_infected)
print(rpt_tte_uninfected)


# =====================================================================
# 3. REPEATABILITY BY TREATMENT  (infection statuses pooled)
# ---------------------------------------------------------------------
# CAUTION: in CON, only 120/180 trials produced an emergence (infected
# snails emerged in just 47/90). Heavy censoring means the CON estimate
# rests on less information — expect a wider CrI, and say so when
# reporting it.
# =====================================================================

fit_by_trt <- function(trt) {
  brm(time_to_emerge_cens | cens(censored) ~ Infection.status + (1 | Snail_ID),
      family = lognormal(),
      data   = subset(total_snails, Treatment == trt),
      chains = CHAINS, cores = 4, iter = ITER, warmup = WARMUP,
      control = list(adapt_delta = 0.95), seed = 1)
}

m_tte_CON <- fit_by_trt("CON")
m_tte_INV <- fit_by_trt("INV")
m_tte_NAT <- fit_by_trt("NAT")

rpt_tte_CON <- get_R(m_tte_CON)
rpt_tte_INV <- get_R(m_tte_INV)
rpt_tte_NAT <- get_R(m_tte_NAT)

print(rbind(CON = rpt_tte_CON, INV = rpt_tte_INV, NAT = rpt_tte_NAT))


# =====================================================================
# 4. INDIVIDUAL EFFECTS (BLUPs, descriptive only)
# ---------------------------------------------------------------------
# The brms equivalent of BLUPs, taken from m_tte (no new model needed).
# Extracted here for the caterpillar plot in 06_figures.R.
#
# Sign convention: on the LOG-TIME scale, a POSITIVE value means the
# snail took LONGER than average to emerge (SHYER) — the opposite
# direction from the climbing-height BLUPs (positive = climbed higher).
# Do not put the two on the same axis without flipping one.
# =====================================================================

re_tte <- ranef(m_tte)$Snail_ID[, , "Intercept"] %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Snail_ID") %>%
  rename(fit = Estimate, lwr = Q2.5, upr = Q97.5)

re_tte$infection <- total_snails$Infection.status[
  match(re_tte$Snail_ID, total_snails$Snail_ID)]


# =====================================================================
# 5. CROSS-CONTEXT CONSISTENCY  (CON vs NAT / INV)
# ---------------------------------------------------------------------
# This CANNOT be done with raw per-snail means (as for climbing height):
# a plain mean would drop the censored trials or treat them as exactly
# 20 min, and the damage is uneven (infected snails in CON emerged in
# only 47/90 trials). Instead, fit one model per treatment with a random
# intercept per snail and correlate the individual estimates, which
# already account for the censoring.
#
# CAVEAT for the paper: correlating point estimates of random effects
# ignores their uncertainty and shrinks them toward the mean, so these
# correlations are approximate — report as descriptive. A proper version
# would fit one multivariate model with treatment-specific random
# intercepts. (The models m_tte_CON/INV/NAT from section 3 are reused.)
# =====================================================================

get_ind <- function(model, label) {
  ranef(model)$Snail_ID[, , "Intercept"] %>%
    as.data.frame() %>%
    tibble::rownames_to_column("Snail_ID") %>%
    select(Snail_ID, !!label := Estimate)
}

ind_wide <- get_ind(m_tte_CON, "CON") %>%
  full_join(get_ind(m_tte_NAT, "NAT"), by = "Snail_ID") %>%
  full_join(get_ind(m_tte_INV, "INV"), by = "Snail_ID") %>%
  mutate(Infection.status = total_snails$Infection.status[
    match(Snail_ID, total_snails$Snail_ID)])

cor.test(ind_wide$CON, ind_wide$NAT)
cor.test(ind_wide$CON, ind_wide$INV)


# =====================================================================
# 6. PLASTICITY — random slopes
# ---------------------------------------------------------------------
# Plasticity = the slope of an individual's response to the cue,
# estimated with a RANDOM SLOPE (not a difference of separately-fitted
# intercepts). This one model answers:
#   Q1 Do individuals differ in plasticity? -> sd(TreatmentINV/NAT) > 0?
#   Q2 Is plasticity linked to behavioural type? -> cor(Intercept, slope)
#   Q3 Does infection change the cue response? -> already answered: the
#      interaction in m_tte_int was not supported by LOO.
# Censoring is handled inside the model, so all 540 trials count.
# =====================================================================

m_tte_slope <- brm(
  time_to_emerge_cens | cens(censored) ~ Treatment + Infection.status +
    (Treatment | Snail_ID),
  family  = lognormal(),
  data    = total_snails,
  chains  = CHAINS, cores = 4, iter = ITER, warmup = WARMUP,
  control = list(adapt_delta = 0.99, max_treedepth = 12),
  seed    = 1
)

summary(m_tte_slope)
rstan::check_hmc_diagnostics(m_tte_slope$fit)   # check divergences first

# Random-slope vs random-intercept only
loo_compare(loo(m_tte), loo(m_tte_slope))

# Individual plasticity estimates (one model -> common scale).
# Negative slope = emerges FASTER under the cue than the average snail.
re <- ranef(m_tte_slope)$Snail_ID

plast <- data.frame(
  Snail_ID       = dimnames(re)[[1]],
  intercept      = re[, "Estimate", "Intercept"],
  plasticity_INV = re[, "Estimate", "TreatmentINV"],
  plasticity_NAT = re[, "Estimate", "TreatmentNAT"]
) %>%
  mutate(Infection.status = total_snails$Infection.status[
    match(Snail_ID, total_snails$Snail_ID)])

# Descriptive only (formal infection test is the interaction in m_tte_int)
plast %>%
  group_by(Infection.status) %>%
  summarise(mean_INV = mean(plasticity_INV),
            mean_NAT = mean(plasticity_NAT),
            sd_INV   = sd(plasticity_INV),
            sd_NAT   = sd(plasticity_NAT),
            .groups  = "drop")

d <- as_draws_df(m_rn_inf)$cor_Snail_ID__TreatmentINV__TreatmentNAT -
  as_draws_df(m_rn_uninf)$cor_Snail_ID__TreatmentINV__TreatmentNAT
cat("diff:", round(mean(d),3), " 95% [", round(quantile(d,.025),3), ",", round(quantile(d,.975),3), "]\n")
# =====================================================================
# 7. REACTION-NORM MODELS BY INFECTION STATUS
# ---------------------------------------------------------------------
# Random-slope models fitted separately within each infection status,
# used for the reaction-norm figure in 06_figures.R.
# =====================================================================

fit_rn <- function(status) {
  brm(time_to_emerge_cens | cens(censored) ~ Treatment + (1 + Treatment | Snail_ID),
      family  = lognormal(),
      data    = subset(total_snails, Infection.status == status),
      chains  = CHAINS, cores = 4, iter = ITER, warmup = WARMUP,
      control = list(adapt_delta = 0.99, max_treedepth = 12),
      seed    = 123)
}

m_rn_inf   <- fit_rn("Infected")
m_rn_uninf <- fit_rn("Uninfected")

summary(m_rn_inf)
summary(m_rn_uninf)


# =====================================================================
# 8. PREDICTABILITY — double-hierarchical model
# ---------------------------------------------------------------------
# EXPLORATORY: distributional (sigma) model; treat with caution and
# check convergence before reporting (kept consistent with the
# climbing-height predictability model in 02).
# =====================================================================

bf_pred_tte <- bf(
  time_to_emerge_cens | cens(censored) ~ Treatment + Infection.status +
    (1 |q| Snail_ID),
  sigma ~ Treatment + Infection.status + (1 |q| Snail_ID)
)

m_pred_tte <- brm(
  bf_pred_tte,
  family  = lognormal(),
  data    = total_snails,
  chains  = CHAINS, cores = 4, iter = ITER, warmup = WARMUP,
  control = list(adapt_delta = 0.99),
  seed    = 123
)

summary(m_pred_tte)

pred_cor_tte <- as_draws_df(m_pred_tte)$cor_Snail_ID__Intercept__sigma_Intercept
report_cor(pred_cor_tte, "Level x predictability (all snails)")

# =====================================================================
# PREDICTABILITY (time to emerge), by infection status
# ---------------------------------------------------------------------
# Double-hierarchical (distributional) model: each snail gets its own
# residual SD (sigma) as well as its own mean. The |q| link estimates
# the among-individual correlation between behavioural level (mean
# latency) and predictability (residual SD).
#
# sigma uses a log link, so a larger sigma-intercept = more variable =
# LESS predictable. A positive level-predictability correlation
# therefore means slower-emerging snails are less predictable.
#
# Fitted separately within each infection status to test whether
# predictability differs between infected and uninfected snails.
#
# WARNING: with only 30 snails / 270 trials per group plus an extra
# sigma random-effect layer, these models may diverge or fail to
# initialise. Check convergence before trusting the estimates.
# =====================================================================

fit_pred_tte_group <- function(status) {
  brm(
    bf(time_to_emerge_cens | cens(censored) ~ Treatment + (1 |q| Snail_ID),
       sigma                                ~ Treatment + (1 |q| Snail_ID)),
    family  = lognormal(),
    data    = subset(total_snails, Infection.status == status),
    chains  = 4, cores = 4, iter = 4000, warmup = 1000,
    control = list(adapt_delta = 0.99), seed = 123
  )
}

m_pred_tte_inf   <- fit_pred_tte_group("Infected")
m_pred_tte_uninf <- fit_pred_tte_group("Uninfected")

# Check convergence FIRST (divergences make the estimates unreliable)
rstan::check_hmc_diagnostics(m_pred_tte_inf$fit)
rstan::check_hmc_diagnostics(m_pred_tte_uninf$fit)

# Extract the level x predictability correlation from each group
pred_cor_tte_inf   <- as_draws_df(m_pred_tte_inf)$cor_Snail_ID__Intercept__sigma_Intercept
pred_cor_tte_uninf <- as_draws_df(m_pred_tte_uninf)$cor_Snail_ID__Intercept__sigma_Intercept

cat("predictability (time to emerge):\n")
cat("  Infected:   r =", round(mean(pred_cor_tte_inf), 3),
    " 95% [", round(quantile(pred_cor_tte_inf, 0.025), 3), ", ",
    round(quantile(pred_cor_tte_inf, 0.975), 3), "]\n", sep = "")
cat("  Uninfected: r =", round(mean(pred_cor_tte_uninf), 3),
    " 95% [", round(quantile(pred_cor_tte_uninf, 0.025), 3), ", ",
    round(quantile(pred_cor_tte_uninf, 0.975), 3), "]\n", sep = "")

# Between-group difference (paired across posterior draws)
d_pred_tte <- pred_cor_tte_inf - pred_cor_tte_uninf
cat("  difference (Infected - Uninfected): ", round(mean(d_pred_tte), 3),
    " 95% [", round(quantile(d_pred_tte, 0.025), 3), ", ",
    round(quantile(d_pred_tte, 0.975), 3), "]\n", sep = "")
# =====================================================================
# 9. SAVE FITTED MODELS
# ---------------------------------------------------------------------
# 06_figures.R reads these back to build the plots (and the joint
# reaction-norm figure that also uses the climbing-height slope model).
# outputs/models/*.rds is git-ignored.
# =====================================================================

saveRDS(
  list(
    m_tte   = m_tte,     m_int   = m_tte_int,
    m_inf   = m_tte_inf, m_uninf = m_tte_uninf,
    m_CON   = m_tte_CON, m_INV   = m_tte_INV,  m_NAT = m_tte_NAT,
    m_slope = m_tte_slope,
    m_rn_inf = m_rn_inf, m_rn_uninf = m_rn_uninf,
    m_pred  = m_pred_tte,
    re_tte  = re_tte,    plast = plast,      ind_wide = ind_wide,
    R_overall = rpt_tte_overall,
    R_inf = rpt_tte_infected, R_uninf = rpt_tte_uninfected,
    R_CON = rpt_tte_CON, R_INV = rpt_tte_INV, R_NAT = rpt_tte_NAT
  ),
  "outputs/models/tte_models.rds"
)
