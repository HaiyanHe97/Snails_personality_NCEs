# =====================================================================
# 03  TIME TO EMERGE — right-censored at 20 min (all 540 trials)
# ---------------------------------------------------------------------
# rptR cannot handle censored data (it would silently drop the 134
# non-emergers: 98 infected, 36 uninfected). So everything is fitted with
# censored lognormal brms models, keeping all 540 trials. Repeatability is
# read with get_R() (from 01), on the log scale:
#     R = sd(Snail_ID)^2 / (sd(Snail_ID)^2 + sigma^2).
#
# Effect model is the a-priori FULL model with Treatment x Infection —
# no additive/interaction model selection, no LOO. m_tte supplies both
# the fixed effects and the overall repeatability. Plasticity and
# predictability use additive fixed effects (adjustment only).
#
# Figures live in 06_figures.R. NOTE: 06 reads models$m_tte.
# =====================================================================

source("01_setup.R")


# =====================================================================
# 1. MAIN MODEL — Treatment x Infection (a-priori full model)
# =====================================================================

m_tte <- brm(
  time_to_emerge_cens | cens(censored) ~ Treatment * Infection.status +
    (1 | Snail_ID),
  family  = lognormal(),
  data    = total_snails,
  chains  = CHAINS, cores = 4, iter = ITER, warmup = WARMUP,
  control = list(adapt_delta = 0.95), seed = 1
)

summary(m_tte)          # incl. the two Treatment:Infection interaction rows
pp_check(m_tte)

# Directional evidence on the interaction (for the Results text)
hypothesis(m_tte, "TreatmentINV:Infection.statusInfected < 0")

# Pairwise contrasts = marginal means (emmeans averages over the other factor)
emm_tte_trt <- emmeans(m_tte, ~ Treatment);        pairs(emm_tte_trt)
emm_tte_inf <- emmeans(m_tte, ~ Infection.status); pairs(emm_tte_inf)

# Overall repeatability (from the fitted model; fixed effects only adjust
# the mean, so numerically identical to an additive fit — Table 1 unaffected)
rpt_tte_overall <- get_R(m_tte)
print(rpt_tte_overall)


# =====================================================================
# 2. REPEATABILITY BY INFECTION STATUS  (treatments pooled)
# Infection held constant -> Treatment only, no interaction to add.
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
# 3. REPEATABILITY BY TREATMENT  (infection pooled)
# Treatment held constant -> Infection only, no interaction to add.
# CAUTION: CON is heavily censored (infected emerged in only 47/90), so
# its estimate rests on less information -> wider CrI; say so when reporting.
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
# 5. CROSS-CONTEXT CONSISTENCY  (CON vs NAT / INV)
# Raw means would mishandle censoring (uneven: CON infected 47/90), so
# correlate per-snail random intercepts from the section-3 models.
# Descriptive: correlating point estimates ignores their uncertainty.
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
# 6. PLASTICITY — random slopes  (fixed effects additive; adjustment)
#   Q1 Do individuals differ in plasticity? -> sd(TreatmentINV/NAT) CrI > 0?
#   Q2 Plasticity linked to type?           -> cor(Intercept, slope)
#   Q3 Infection x cue response             -> the interaction in m_tte
# (No LOO / model selection, consistent with the a-priori approach.)
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
# 7. REACTION-NORM MODELS BY INFECTION STATUS  (for the figure in 06)
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

# Slope-slope (INV vs NAT) correlation: infected vs uninfected
# (MOVED here from section 6, where it ran before these models existed.)
d <- as_draws_df(m_rn_inf)$cor_Snail_ID__TreatmentINV__TreatmentNAT -
     as_draws_df(m_rn_uninf)$cor_Snail_ID__TreatmentINV__TreatmentNAT
cat("slope-slope diff (Infected - Uninfected):", round(mean(d), 3),
    " 95% [", round(quantile(d, .025), 3), ",", round(quantile(d, .975), 3), "]\n")


# =====================================================================
# 8. PREDICTABILITY — double-hierarchical (sigma) model  (additive)
# EXPLORATORY, consistent with the height DHGLM in 02. sigma has a log
# link, so larger sigma = less predictable; a positive level-predictability
# r = slower snails are less predictable. May diverge — check convergence.
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


pred_cor_tte <- as_draws_df(m_pred_tte)$cor_Snail_ID__Intercept__sigma_Intercept
cat("level x predictability (all snails): r =", round(median(pred_cor_tte), 3),
    " 95% [", round(quantile(pred_cor_tte, .025), 3), ",",
    round(quantile(pred_cor_tte, .975), 3), "]",
    " P(r>0) =", round(mean(pred_cor_tte > 0), 3), "\n")

# ---- by infection status (30 snails/group + sigma layer: may not converge)
fit_pred_tte_group <- function(status) {
  brm(bf(time_to_emerge_cens | cens(censored) ~ Treatment + (1 |q| Snail_ID),
         sigma                                ~ Treatment + (1 |q| Snail_ID)),
      family = lognormal(),
      data   = subset(total_snails, Infection.status == status),
      chains = CHAINS, cores = 4, iter = ITER, warmup = WARMUP,
      control = list(adapt_delta = 0.99), seed = 123)
}

m_pred_tte_inf   <- fit_pred_tte_group("Infected")
m_pred_tte_uninf <- fit_pred_tte_group("Uninfected")

rstan::check_hmc_diagnostics(m_pred_tte_inf$fit)
rstan::check_hmc_diagnostics(m_pred_tte_uninf$fit)

pred_cor_tte_inf   <- as_draws_df(m_pred_tte_inf)$cor_Snail_ID__Intercept__sigma_Intercept
pred_cor_tte_uninf <- as_draws_df(m_pred_tte_uninf)$cor_Snail_ID__Intercept__sigma_Intercept

cat("predictability (time to emerge):\n")
cat("  Infected:   r =", round(mean(pred_cor_tte_inf), 3),
    " 95% [", round(quantile(pred_cor_tte_inf, 0.025), 3), ", ",
    round(quantile(pred_cor_tte_inf, 0.975), 3), "]\n", sep = "")
cat("  Uninfected: r =", round(mean(pred_cor_tte_uninf), 3),
    " 95% [", round(quantile(pred_cor_tte_uninf, 0.025), 3), ", ",
    round(quantile(pred_cor_tte_uninf, 0.975), 3), "]\n", sep = "")

d_pred_tte <- pred_cor_tte_inf - pred_cor_tte_uninf
cat("  difference (Infected - Uninfected): ", round(mean(d_pred_tte), 3),
    " 95% [", round(quantile(d_pred_tte, 0.025), 3), ", ",
    round(quantile(d_pred_tte, 0.975), 3), "]\n", sep = "")


# =====================================================================
# 9. SAVE FITTED MODELS  (m_int key removed; 06 reads m_tte)
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
    R_inf = rpt_tte_infected, R_uninf = rpt_tte_uninfected,
    R_CON = rpt_tte_CON, R_INV = rpt_tte_INV, R_NAT = rpt_tte_NAT
  ),
  "outputs/models/tte_models.rds"
)




# =====================================================================
# TEST: one model for infection-split repeatability (TIME TO EMERGE)
# ---------------------------------------------------------------------
# Same idea as the climbing-height test, adapted for time to emerge:
#   - family = lognormal()  (was gaussian)
#   - response is time_to_emerge_cens | cens(censored)  (was max_height | cens(h_cens))
# Both among-individual SD and residual SD are split by infection, to
# match what the two subset models (m_tte_inf / m_tte_uninf) did.
#
# CAUTION: time to emerge is heavily censored (esp. CON), so this model
# may converge less cleanly than the height one. Check diagnostics FIRST.
#
# Run after 03 has been sourced (so m_tte_inf/uninf, R_h_* exist), OR
# standalone after source("01_setup.R").
# =====================================================================

# source("01_setup.R")   # uncomment if running standalone

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
rstan::check_hmc_diagnostics(m_tte_byinf$fit)   # <-- CHECK THIS FIRST

# ---- pull the two group-specific R's from the posterior -------------
d <- as_draws_df(m_tte_byinf)
# If a name below errors, run:
#   grep("sd_Snail_ID|b_sigma", names(d), value = TRUE)
# and paste me the real names.
sd_inf   <- d[["sd_Snail_ID__Intercept:Infection.statusInfected"]]
sd_uninf <- d[["sd_Snail_ID__Intercept:Infection.statusUninfected"]]

sig_uninf <- exp(d[["b_sigma_Intercept"]])                                    # reference = Uninfected
sig_inf   <- exp(d[["b_sigma_Intercept"]] + d[["b_sigma_Infection.statusInfected"]])

R_inf_1   <- sd_inf^2   / (sd_inf^2   + sig_inf^2)
R_uninf_1 <- sd_uninf^2 / (sd_uninf^2 + sig_uninf^2)
R_diff    <- R_inf_1 - R_uninf_1

q <- function(x) c(median = median(x),
                   lwr = quantile(x, .025, names = FALSE),
                   upr = quantile(x, .975, names = FALSE))

cat("\n========== SINGLE MODEL (m_tte_byinf) — time to emerge ==========\n")
cat("R Infected   :", round(q(R_inf_1),3),   "\n")
cat("R Uninfected :", round(q(R_uninf_1),3), "\n")
cat("Difference (Inf - Uninf):", round(q(R_diff),3),
    " | P(diff>0) =", round(mean(R_diff > 0),3), "\n")

# ---- side-by-side with the two-subset estimates from 03 -------------
if (exists("rpt_tte_infected") && exists("rpt_tte_uninfected")) {
  cat("\n---- SUBSET models (from 03) for comparison ----\n")
  cat("R Infected   :", round(rpt_tte_infected,3),   "\n")
  cat("R Uninfected :", round(rpt_tte_uninfected,3), "\n")
}
