# =====================================================================
# 02  CLIMBING HEIGHT — height above the waterline, right-censored at 27 cm
# ---------------------------------------------------------------------
# 0 = in water; 0.1 = crossed the waterline, no further climbing;
# 27 = climbed out of the cylinder -> right-censored (92/540).
# Zeros are NOT censored. All 540 trials kept.
#
# Together vs separate (same rule in 02, 03 and 05):
#   Means, and individual stability (within-individual variance,
#   modelled through sigma) -> ONE model, infection as a fixed effect.
#   Among-individual variances and correlations (repeatability,
#   plasticity, syndrome) -> SEPARATE models per group, compared by the
#   posterior of the difference, not by overlap of CrIs.
#
# Effect model = a-priori full model (Treatment x Infection), no model
# selection. Repeatability, plasticity and stability use additive fixed
# effects. Figures live in 06_figures.R (06 reads models$m_int).
# =====================================================================

source("01_setup.R")

# Ceiling hits by cue (context for the censoring)
table(total_snails$h_cens, total_snails$Treatment)


# =====================================================================
# 1. MAIN MODEL — Treatment x Infection   [ONE model: means]
# Do cues and infection relate to mean height, and does the cue effect
# differ between infected and uninfected snails (interaction)?
# =====================================================================

# --- SENSITIVITY (OLD): no censoring, 27 treated as a real value -----
model_height_old <- lmer(max_height_aw ~ Treatment * Infection.status +
                           (1 | Snail_ID), data = total_snails)
summary(model_height_old)
anova(model_height_old, type = 2)

# --- MAIN: right-censored at 27, interaction -------------------------
m_height_int <- brm(
  max_height_aw | cens(h_cens) ~ Treatment * Infection.status + (1 | Snail_ID),
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
# 2. REPEATABILITY (personality)  (fixed effects additive)
# =====================================================================

# ---- 2.1 Overall   [ONE model: do snails differ consistently?] -------

# SENSITIVITY (OLD): rptR, no censoring.
rpt_height_old <- rpt(max_height_aw ~ Treatment + Infection.status + (1 | Snail_ID),
                      grname = "Snail_ID", data = total_snails,
                      datatype = "Gaussian", nboot = 1000, npermut = 0)
summary(rpt_height_old)

# MAIN: overall R from the fitted model. Fixed effects only adjust the
# mean, so this equals an additive fit to 3 decimals (Table 1 unaffected).
R_height_overall <- get_R(m_height_int)
print(R_height_overall)

# ---- 2.2 By infection status   [SEPARATE: R is a variance ratio] ------
# R needs its own among-individual and residual variance in each group.
# Different snails in each group -> independent posteriors, so the
# difference of draws is a valid comparison.

fit_h_inf <- function(status) {
  brm(max_height_aw | cens(h_cens) ~ Treatment + (1 | Snail_ID),
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

R_h_diff <- get_R_draws(m_h_inf) - get_R_draws(m_h_uninf)
summ_diff(R_h_diff)                         # Infected - Uninfected

# ---- 2.3 By treatment   [SEPARATE: R per cue context] ----------------
# Same snails in all treatments: differences between separate fits ignore
# this dependence, so their CrIs are approximate.

fit_h_trt <- function(trt) {
  brm(max_height_aw | cens(h_cens) ~ Infection.status + (1 | Snail_ID),
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

summ_diff(get_R_draws(m_h_INV) - get_R_draws(m_h_CON))   # INV - CON
summ_diff(get_R_draws(m_h_NAT) - get_R_draws(m_h_CON))   # NAT - CON
summ_diff(get_R_draws(m_h_NAT) - get_R_draws(m_h_INV))   # NAT - INV


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
  summarise(mean_height = mean(max_height_aw), .groups = "drop")

snail_wide <- snail_mean %>%
  pivot_wider(names_from = Treatment, values_from = mean_height)

cor.test(snail_wide$CON, snail_wide$NAT)
cor.test(snail_wide$CON, snail_wide$INV)


# =====================================================================
# 5. INDIVIDUAL PLASTICITY — random slopes  (fixed effects additive)
# =====================================================================

# ---- 5.1 Overall   [ONE model: do snails differ in plasticity?] -------
# Evidence = sd(TreatmentINV), sd(TreatmentNAT) clearly above zero;
# also cor(Intercept, slope). Whether infection changes the MEAN cue
# response is the interaction in m_height_int (section 1).

m_height_slope <- brm(
  max_height_aw | cens(h_cens) ~ Treatment + Infection.status +
    (1 + Treatment | Snail_ID),
  family  = gaussian(),
  data    = total_snails,
  chains  = CHAINS, cores = 4, iter = ITER, warmup = WARMUP,
  control = list(adapt_delta = 0.99, max_treedepth = 12),
  seed    = 123
)

summary(m_height_slope)
rstan::check_hmc_diagnostics(m_height_slope$fit)

# ---- 5.2 By infection status   [SEPARATE: compare correlations] -----
# Slope-slope correlation (INV vs NAT) per group, and its difference.

fit_rn_cov <- function(status) {
  brm(bf(max_height_aw | cens(h_cens) ~ Treatment + (1 + Treatment |q| Snail_ID)),
      data = subset(total_snails, Infection.status == status),
      family = gaussian(),
      chains = CHAINS, cores = 4, iter = ITER, warmup = WARMUP,
      control = list(adapt_delta = 0.99), seed = 123)
}
model_reaction_norm_inf   <- fit_rn_cov("Infected")
model_reaction_norm_uninf <- fit_rn_cov("Uninfected")

summary(model_reaction_norm_inf);   rhat(model_reaction_norm_inf)
summary(model_reaction_norm_uninf); rhat(model_reaction_norm_uninf)

posterior_INVNAT_inf   <- as_draws_df(model_reaction_norm_inf)$cor_Snail_ID__TreatmentINV__TreatmentNAT
posterior_INVNAT_uninf <- as_draws_df(model_reaction_norm_uninf)$cor_Snail_ID__TreatmentINV__TreatmentNAT
diff_INVNAT <- posterior_INVNAT_inf - posterior_INVNAT_uninf

summ_draws(posterior_INVNAT_inf)            # Infected
summ_draws(posterior_INVNAT_uninf)          # Uninfected
summ_diff(diff_INVNAT)                      # Infected - Uninfected


# =====================================================================
# 6. INDIVIDUAL STABILITY — double-hierarchical model (additive)
#    [ONE model: stability is a within-individual variance, so infected
#     vs uninfected is the infection effect on sigma]
# Height in dm (h_dm, built in 01): on the cm scale the sampler failed to
# initialise. Rescaling does not change correlations or sigma effects.
# sigma-intercept prior normal(0, 0.5) on dm = normal(2.3, 0.5) on cm.
# MAIN: all 540 trials. SENSITIVITY: emerged trials only (406), as the
# STRICT syndrome model in 05, to check the effect of the 0 cm floor.
# =====================================================================

fit_stab <- function(dat) {
  brm(
    bf(h_dm | cens(h_cens) ~ Treatment + Infection.status + (1 |q| Snail_ID),
       sigma ~ Treatment + Infection.status + (1 |q| Snail_ID)),
    family = gaussian(), data = dat,
    prior = c(prior(normal(0, 0.5), class = "Intercept", dpar = "sigma"),
              prior(normal(0, 0.5), class = "b",         dpar = "sigma"),
              prior(exponential(2), class = "sd",        dpar = "sigma")),
    chains = CHAINS, cores = 4, iter = ITER, warmup = WARMUP,
    init = 0, control = list(adapt_delta = 0.999), seed = 123)
}

# ---- MAIN (540 trials) ----------------------------------------------
m_height_pred <- fit_stab(total_snails)
summary(m_height_pred)
rstan::check_hmc_diagnostics(m_height_pred$fit)
summ_draws(as_draws_df(m_height_pred)$cor_Snail_ID__Intercept__sigma_Intercept)  # mean x stability
summ_draws(as_draws_df(m_height_pred)$b_sigma_Infection.statusInfected)          # Infected vs Uninfected

# ---- SENSITIVITY (emerged trials only, 406) --------------------------
m_height_pred_em <- fit_stab(subset(total_snails, emerged == 1))
summary(m_height_pred_em)
rstan::check_hmc_diagnostics(m_height_pred_em$fit)
summ_draws(as_draws_df(m_height_pred_em)$cor_Snail_ID__Intercept__sigma_Intercept)
summ_draws(as_draws_df(m_height_pred_em)$b_sigma_Infection.statusInfected)


# =====================================================================
# 7. SAVE FITTED MODELS  (06 reads models$m_int)
# =====================================================================

saveRDS(
  list(
    m_int   = m_height_int,
    m_slope = m_height_slope, m_pred = m_height_pred, m_pred_em = m_height_pred_em,
    m_inf   = m_h_inf,       m_uninf = m_h_uninf,
    m_CON   = m_h_CON,       m_INV   = m_h_INV,   m_NAT = m_h_NAT,
    blups   = blups,
    R_overall = R_height_overall,
    R_inf = R_h_inf,  R_uninf = R_h_uninf, R_diff = R_h_diff,
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
# CHECK (not reported): repeatability by infection in ONE model with
# group-specific variances, gr(Snail_ID, by = Infection.status).
# Gives R similar to the subset models in 2.2 (reviewer back-up).
# Run after sections 1-2, or standalone after source("01_setup.R").
# =====================================================================

m_h_byinf <- brm(
  bf(max_height_aw | cens(h_cens) ~ Treatment + Infection.status +
       (1 | gr(Snail_ID, by = Infection.status)),   # among-ind SD per group
     sigma ~ Infection.status),                      # residual SD per group
  family  = gaussian(),
  data    = total_snails,
  chains  = CHAINS, cores = 4, iter = ITER, warmup = WARMUP,
  control = list(adapt_delta = 0.95, max_treedepth = 12),
  seed    = 123
)

summary(m_h_byinf)
rstan::check_hmc_diagnostics(m_h_byinf$fit)

d <- as_draws_df(m_h_byinf)
# if a name errors: print(grep("sd_Snail_ID", names(d), value = TRUE))
sd_inf    <- d[["sd_Snail_ID__Intercept:Infection.statusInfected"]]
sd_uninf  <- d[["sd_Snail_ID__Intercept:Infection.statusUninfected"]]
sig_uninf <- exp(d[["b_sigma_Intercept"]])                                    # reference = Uninfected
sig_inf   <- exp(d[["b_sigma_Intercept"]] + d[["b_sigma_Infection.statusInfected"]])

R_inf_1   <- sd_inf^2   / (sd_inf^2   + sig_inf^2)
R_uninf_1 <- sd_uninf^2 / (sd_uninf^2 + sig_uninf^2)

summ_draws(R_inf_1)                          # Infected
summ_draws(R_uninf_1)                        # Uninfected
summ_diff(R_inf_1 - R_uninf_1)               # Infected - Uninfected
if (exists("R_h_diff")) summ_diff(R_h_diff)  # same, from the subset models
