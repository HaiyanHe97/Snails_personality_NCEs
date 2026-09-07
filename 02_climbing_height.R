# =====================================================================
# 02  CLIMBING HEIGHT — right-censored at 30 cm
# ---------------------------------------------------------------------
# A snail recorded at 30 cm climbed out of the 30 cm cylinder: its true
# height is >= 30 and could not be observed -> right-censored (92/540).
# Zeros are NOT censored and kept as observed. All 540 trials kept.
#
# MAIN:        right-censored Gaussian model (brms).
# SENSITIVITY: old lmer, 30 treated as a real value.
#
# Effect model is the a-priori FULL model with Treatment x Infection —
# no additive/interaction model selection, no LOO. The interaction is
# reported as-is (it is the research question). Repeatability, plasticity
# and predictability use additive fixed effects (adjustment only).
#
# Figures live in 06_figures.R. NOTE: 06 reads models$m_int.
# =====================================================================

source("01_setup.R")

# Ceiling hits by cue (context for the censoring)
table(total_snails$h_cens, total_snails$Treatment)


# =====================================================================
# 1. MAIN MODEL — Treatment x Infection (a-priori full model)
# =====================================================================

# --- SENSITIVITY (OLD): no censoring, 30 treated as a real value -----
model_height_old <- lmer(max_height ~ Treatment * Infection.status +
                           (1 | Snail_ID), data = total_snails)
summary(model_height_old)
anova(model_height_old, type = 2)

# --- MAIN: right-censored at 30, interaction -------------------------
m_height_int <- brm(
  max_height | cens(h_cens) ~ Treatment * Infection.status + (1 | Snail_ID),
  family  = gaussian(),
  data    = total_snails,
  chains  = CHAINS, cores = 4, iter = ITER, warmup = WARMUP,
  control = list(adapt_delta = 0.95), seed = 123
)

# Table 2a: fixed effects, INCLUDING the two Treatment:Infection rows.
summary(m_height_int)

# PP-check will look imperfect: 176/540 heights are exactly 0 and a
# Gaussian cannot spike at zero. Censoring fixes the ceiling, not the
# floor; effects are large and consistent -> note in limitations.
pp_check(m_height_int)

# Table 2b: pairwise contrasts = MARGINAL means (emmeans averages over
# the other factor under the interaction model).
emm_h_trt <- emmeans(m_height_int, ~ Treatment);        pairs(emm_h_trt)
emm_h_inf <- emmeans(m_height_int, ~ Infection.status); pairs(emm_h_inf)


# =====================================================================
# 2. REPEATABILITY  (fixed effects additive — adjustment only)
# =====================================================================

# ---- 2.1 Overall ----------------------------------------------------

# SENSITIVITY (OLD): rptR, no censoring.
rpt_height_old <- rpt(max_height ~ Treatment + Infection.status + (1 | Snail_ID),
                      grname = "Snail_ID", data = total_snails,
                      datatype = "Gaussian", nboot = 1000, npermut = 0)
summary(rpt_height_old)

# MAIN: overall R from the fitted model. R = var(Snail_ID)/(var+sigma^2);
# fixed effects only adjust the mean, so this is numerically identical to
# an additive fit to 3 decimals — Table 1 unaffected.
R_height_overall <- get_R(m_height_int)
print(R_height_overall)




# ---- 2.2 By infection status (Treatment only after subsetting) ------

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


# ---- 2.3 By treatment (Infection only after subsetting) -------------

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
# Compare with uncensored rptR (0.28 / 0.32 / 0.46). NAT has the most
# ceiling hits, so its estimate is the one most likely to shift.


# =====================================================================
# 3. INDIVIDUAL EFFECTS (BLUPs, descriptive; from m_height_int)
# =====================================================================

blups <- ranef(m_height_int)$Snail_ID[, , "Intercept"] %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Snail_ID") %>%
  rename(fit = Estimate, lwr = Q2.5, upr = Q97.5)

blups$infection <- total_snails$Infection.status[
  match(blups$Snail_ID, total_snails$Snail_ID)]


# =====================================================================
# 4. CROSS-CONTEXT CONSISTENCY (raw per-snail means; model-free)
# 92 ceiling values slightly attenuate these; ranking preserved, no refit.
# =====================================================================

snail_mean <- total_snails %>%
  group_by(Snail_ID, Treatment, Infection.status) %>%
  summarise(mean_height = mean(max_height), .groups = "drop")

snail_wide <- snail_mean %>%
  pivot_wider(names_from = Treatment, values_from = mean_height)

cor.test(snail_wide$CON, snail_wide$NAT)
cor.test(snail_wide$CON, snail_wide$INV)


# =====================================================================
# 5. PLASTICITY — random slopes  (fixed effects additive; adjustment)
# Plasticity evidence = sd(TreatmentINV), sd(TreatmentNAT); a CrI clear
# of zero means individuals differ in plasticity. (No LOO / model
# selection, consistent with the a-priori approach.)
# =====================================================================

m_height_slope <- brm(
  max_height | cens(h_cens) ~ Treatment + Infection.status +
    (1 + Treatment | Snail_ID),
  family  = gaussian(),
  data    = total_snails,
  chains  = CHAINS, cores = 4, iter = ITER, warmup = WARMUP,
  control = list(adapt_delta = 0.99, max_treedepth = 12),
  seed    = 123
)

summary(m_height_slope)
rstan::check_hmc_diagnostics(m_height_slope$fit)

# ---- Fit by infection status (for reaction-norm figure / cor) -------
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

# ---- Slope-slope correlation (INV vs NAT): infected vs uninfected ---
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
# 6. PREDICTABILITY — double-hierarchical model  (additive; exploratory)
# Unstable in earlier runs — check convergence before reporting.
# =====================================================================

m_height_pred <- brm(
  bf(max_height | cens(h_cens) ~ Treatment + Infection.status + (1 |q| Snail_ID),
     sigma ~ Treatment + Infection.status + (1 |q| Snail_ID)),
  family = gaussian(),
  data   = total_snails,
  prior = c(
    prior(normal(2.3, 0.5), class = "Intercept", dpar = "sigma"),
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
# 7. SAVE FITTED MODELS  (06 reads models$m_int)
# =====================================================================

saveRDS(
  list(
    m_int   = m_height_int,
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



# =====================================================================
# TEST: one model for infection-split repeatability (climbing height)
# ---------------------------------------------------------------------
# Goal: reproduce the by-infection R (currently from two subset models)
# in a SINGLE model, and get the between-group difference with a CI.
#
# Both the among-individual SD and the residual SD are allowed to differ
# by infection status, so this matches what the two subset models did
# (each subset estimated its own among-individual and residual variance).
#
# Run after 02 has been sourced / m_h_inf & m_h_uninf exist, OR stand-
# alone after source("01_setup.R").
# =====================================================================

# source("01_setup.R")   # uncomment if running standalone

# ---- single model: variances split by infection ---------------------
m_h_byinf <- brm(
  bf(max_height | cens(h_cens) ~ Treatment + Infection.status +
       (1 | gr(Snail_ID, by = Infection.status)),   # among-ind SD per group
     sigma ~ Infection.status),                      # residual SD per group
  family  = gaussian(),
  data    = total_snails,
  chains  = CHAINS, cores = 4, iter = ITER, warmup = WARMUP,
  control = list(adapt_delta = 0.95, max_treedepth = 12),
  seed    = 123
)

summary(m_h_byinf)
rstan::check_hmc_diagnostics(m_h_byinf$fit)   # watch for divergences

# ---- pull the two group-specific R's from the posterior -------------
d <- as_draws_df(m_h_byinf)

# among-individual SDs (brms names them by the 'by' level)
# check names if this errors — print(grep("sd_Snail_ID", names(d), value=TRUE))
sd_inf   <- d[["sd_Snail_ID__Intercept:Infection.statusInfected"]]
sd_uninf <- d[["sd_Snail_ID__Intercept:Infection.statusUninfected"]]

# residual SDs: sigma is on a log link here.
# log(sigma) = b_sigma_Intercept (+ b_sigma_Infection.statusInfected for infected)
sig_uninf <- exp(d[["b_sigma_Intercept"]])                                   # reference = Uninfected
sig_inf   <- exp(d[["b_sigma_Intercept"]] + d[["b_sigma_Infection.statusInfected"]])

R_inf_1   <- sd_inf^2   / (sd_inf^2   + sig_inf^2)
R_uninf_1 <- sd_uninf^2 / (sd_uninf^2 + sig_uninf^2)
R_diff    <- R_inf_1 - R_uninf_1

q <- function(x) c(median = median(x),
                   lwr = quantile(x, .025, names = FALSE),
                   upr = quantile(x, .975, names = FALSE))

cat("\n================ SINGLE MODEL (m_h_byinf) ================\n")
cat("R Infected   :", round(q(R_inf_1),3),   "\n")
cat("R Uninfected :", round(q(R_uninf_1),3), "\n")
cat("Difference (Inf - Uninf):", round(q(R_diff),3),
    " | P(diff>0) =", round(mean(R_diff > 0),3), "\n")

# ---- side-by-side with the two-subset estimates ---------------------
# (only runs if the subset models from 02 are in memory)
if (exists("R_h_inf") && exists("R_h_uninf")) {
  cat("\n---- SUBSET models (from 02) for comparison ----\n")
  cat("R Infected   :", round(R_h_inf,3),   "\n")
  cat("R Uninfected :", round(R_h_uninf,3), "\n")
}
