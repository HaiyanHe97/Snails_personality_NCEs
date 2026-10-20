# =====================================================================
# 03  TIME TO EMERGE — right-censored at 20 min (all 540 trials)
# ---------------------------------------------------------------------
# rptR cannot handle censored data (it would silently drop the 134
# non-emergers: 98 infected, 36 uninfected). So everything is fitted with
# censored lognormal brms models, keeping all 540 trials. Repeatability is
# read with get_R() (from 01), on the log scale:
#     R = sd(Snail_ID)^2 / (sd(Snail_ID)^2 + sigma^2).
#
# Together vs separate: same rule as in 02 (see the header there).
#
# Effect model = a-priori full model (Treatment x Infection), no model
# selection. Plasticity and individual stability use additive fixed
# effects. Figures live in 06_figures.R (06 reads models$m_tte).
# =====================================================================

source("01_setup.R")


# =====================================================================
# 1. MAIN MODEL — Treatment x Infection   [ONE model: means]
# =====================================================================

m_tte <- brm(
  time_to_emerge_cens | cens(censored) ~ Treatment * Infection.status +
    (1 | Snail_ID),
  family  = lognormal(),
  data    = total_snails,
  chains  = CHAINS, cores = 4, iter = ITER, warmup = WARMUP,
  control = list(adapt_delta = 0.95), seed = 1
)

summary(m_tte)

# Table 1a: fixed effects, INCLUDING the two Treatment:Infection rows.
# robust = TRUE gives posterior medians (as in all other tables).
round(fixef(m_tte, robust = TRUE), 2)
pp_check(m_tte)

# Directional evidence on the interaction (for the Results text)
hypothesis(m_tte, "TreatmentINV:Infection.statusInfected < 0")

# Pairwise contrasts = marginal means (emmeans averages over the other factor)
emm_tte_trt <- emmeans(m_tte, ~ Treatment);        pairs(emm_tte_trt)
emm_tte_inf <- emmeans(m_tte, ~ Infection.status); pairs(emm_tte_inf)

# Overall repeatability   [ONE model: do snails differ consistently?]
# Fixed effects only adjust the mean, so numerically identical to an
# additive fit (Table 2 unaffected).
rpt_tte_overall <- get_R(m_tte)
print(rpt_tte_overall)


# =====================================================================
# 2. REPEATABILITY BY INFECTION STATUS   [SEPARATE: R is a variance ratio]
# Infection held constant -> Treatment only. Different snails in each
# group -> independent posteriors, so the difference of draws is valid.
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

R_tte_diff <- get_R_draws(m_tte_inf) - get_R_draws(m_tte_uninf)
summ_diff(R_tte_diff)                       # Infected - Uninfected


# =====================================================================
# 3. REPEATABILITY BY TREATMENT   [SEPARATE: R per cue context]
# Treatment held constant -> Infection only.
# CAUTION: CON is heavily censored (infected emerged in only 47/90), so
# its estimate rests on less information -> wider CrI; say so when reporting.
# Same snails in all treatments: differences between separate fits ignore
# this dependence, so their CrIs are approximate.
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

summ_diff(get_R_draws(m_tte_INV) - get_R_draws(m_tte_CON))   # INV - CON
summ_diff(get_R_draws(m_tte_NAT) - get_R_draws(m_tte_CON))   # NAT - CON
summ_diff(get_R_draws(m_tte_NAT) - get_R_draws(m_tte_INV))   # NAT - INV


# =====================================================================
# 4. INDIVIDUAL EFFECTS (BLUPs, descriptive; from m_tte)
# Sign: on the LOG-TIME scale a POSITIVE value = took LONGER (shyer) —
# opposite to the height BLUPs. Flip one before sharing an axis.
# =====================================================================

re_tte <- ranef(m_tte)$Snail_ID[, , "Intercept"] %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Snail_ID") %>%
  rename(fit = Estimate, lwr = Q2.5, upr = Q97.5)

re_tte$infection <- total_snails$Infection.status[
  match(re_tte$Snail_ID, total_snails$Snail_ID)]





# =====================================================================
# 5. INDIVIDUAL PLASTICITY — random slopes   [ONE model: do snails
#    differ in plasticity?]  (fixed effects additive)
#   Q1 Do individuals differ in plasticity? -> sd(TreatmentINV/NAT) CrI > 0?
#   Q2 Plasticity linked to type?           -> cor(Intercept, slope)
#   Q3 Infection x cue response (mean)      -> the interaction in m_tte
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
rstan::check_hmc_diagnostics(m_tte_slope$fit)

# Individual plasticity estimates (negative slope = faster under the cue)
re <- ranef(m_tte_slope)$Snail_ID
plast <- data.frame(
  Snail_ID       = dimnames(re)[[1]],
  intercept      = re[, "Estimate", "Intercept"],
  plasticity_INV = re[, "Estimate", "TreatmentINV"],
  plasticity_NAT = re[, "Estimate", "TreatmentNAT"]
) %>%
  mutate(Infection.status = total_snails$Infection.status[
    match(Snail_ID, total_snails$Snail_ID)])

# Descriptive (formal infection test = the interaction in m_tte)
plast %>%
  group_by(Infection.status) %>%
  summarise(mean_INV = mean(plasticity_INV), mean_NAT = mean(plasticity_NAT),
            sd_INV = sd(plasticity_INV), sd_NAT = sd(plasticity_NAT),
            .groups = "drop")


# =====================================================================
# 6. REACTION NORMS BY INFECTION STATUS   [SEPARATE: compare correlations]
# Slope-slope correlation (INV vs NAT) per group, and its difference.
# The two models are also used for the reaction-norm figure in 06.
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

rn_cor_inf   <- as_draws_df(m_rn_inf)$cor_Snail_ID__TreatmentINV__TreatmentNAT
rn_cor_uninf <- as_draws_df(m_rn_uninf)$cor_Snail_ID__TreatmentINV__TreatmentNAT

summ_draws(rn_cor_inf)                      # Infected
summ_draws(rn_cor_uninf)                    # Uninfected
summ_diff(rn_cor_inf - rn_cor_uninf)        # Infected - Uninfected


# =====================================================================
# 7. INDIVIDUAL STABILITY — double-hierarchical (sigma) model (additive)
#    [ONE model: stability is a within-individual variance, so infected
#     vs uninfected is the infection effect on sigma]
# sigma has a log link, so a larger sigma = less stable; a positive
# mean-sigma correlation = slower snails are less stable.
# =====================================================================

m_pred_tte <- brm(
  bf(time_to_emerge_cens | cens(censored) ~ Treatment + Infection.status + (1 |q| Snail_ID),
     sigma                                ~ Treatment + Infection.status + (1 |q| Snail_ID)),
  family  = lognormal(),
  data    = total_snails,
  chains  = CHAINS, cores = 4, iter = ITER, warmup = WARMUP,
  control = list(adapt_delta = 0.99), seed = 123
)

summary(m_pred_tte)
rstan::check_hmc_diagnostics(m_pred_tte$fit)
summ_draws(as_draws_df(m_pred_tte)$cor_Snail_ID__Intercept__sigma_Intercept)  # mean x stability
summ_draws(as_draws_df(m_pred_tte)$b_sigma_Infection.statusInfected)          # Infected vs Uninfected


# =====================================================================
# 8. SAVE FITTED MODELS  (06 reads m_tte)
# =====================================================================

saveRDS(
  list(
    m_tte   = m_tte,
    m_inf   = m_tte_inf, m_uninf = m_tte_uninf,
    m_CON   = m_tte_CON, m_INV   = m_tte_INV,  m_NAT = m_tte_NAT,
    m_slope = m_tte_slope,
    m_rn_inf = m_rn_inf, m_rn_uninf = m_rn_uninf,
    m_pred  = m_pred_tte,
    re_tte  = re_tte,    plast = plast,      ind_wide = ind_wide,
    R_overall = rpt_tte_overall,
    R_inf = rpt_tte_infected, R_uninf = rpt_tte_uninfected, R_diff = R_tte_diff,
    R_CON = rpt_tte_CON, R_INV = rpt_tte_INV, R_NAT = rpt_tte_NAT
  ),
  "outputs/models/tte_models.rds"
)


# =====================================================================
# CHECK (not reported): repeatability by infection in ONE model with
# group-specific variances, gr(Snail_ID, by = Infection.status).
# Gives R similar to the subset models in section 2 (reviewer back-up).
# CAUTION: heavily censored (esp. CON); check diagnostics first.
# Run after sections 1-2, or standalone after source("01_setup.R").
# =====================================================================

m_tte_byinf <- brm(
  bf(time_to_emerge_cens | cens(censored) ~ Treatment + Infection.status +
       (1 | gr(Snail_ID, by = Infection.status)),   # among-ind SD per group
     sigma ~ Infection.status),                      # residual SD per group
  family  = lognormal(),
  data    = total_snails,
  chains  = CHAINS, cores = 4, iter = ITER, warmup = WARMUP,
  control = list(adapt_delta = 0.99, max_treedepth = 12),
  seed    = 1
)

summary(m_tte_byinf)
rstan::check_hmc_diagnostics(m_tte_byinf$fit)

d <- as_draws_df(m_tte_byinf)
# if a name errors: print(grep("sd_Snail_ID|b_sigma", names(d), value = TRUE))
sd_inf    <- d[["sd_Snail_ID__Intercept:Infection.statusInfected"]]
sd_uninf  <- d[["sd_Snail_ID__Intercept:Infection.statusUninfected"]]
sig_uninf <- exp(d[["b_sigma_Intercept"]])                                    # reference = Uninfected
sig_inf   <- exp(d[["b_sigma_Intercept"]] + d[["b_sigma_Infection.statusInfected"]])

R_inf_1   <- sd_inf^2   / (sd_inf^2   + sig_inf^2)
R_uninf_1 <- sd_uninf^2 / (sd_uninf^2 + sig_uninf^2)

summ_draws(R_inf_1)                              # Infected
summ_draws(R_uninf_1)                            # Uninfected
summ_diff(R_inf_1 - R_uninf_1)                   # Infected - Uninfected
if (exists("R_tte_diff")) summ_diff(R_tte_diff)  # same, from the subset models
