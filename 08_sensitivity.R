# =====================================================================
# 08  SENSITIVITY ANALYSES
# ---------------------------------------------------------------------
# All sensitivity analyses in one place, in the order they appear in the
# paper. Each refits a MAIN model from 02 or 05 on a subset of the data
# (or without censoring) and is compared with that main model.
#
#   REPORTED
#   1. Table S2  Climbing height, mean effects: only trials in which
#                snails climbed above the waterline (height > 0; 364).
#   2. Text      Syndrome climbing height x time to emerge: only trials
#                in which snails emerged (406). Reported in the text and
#                the Table 3 legend.
#   3. Table S3  Climbing height, individual stability: only trials in
#                which snails climbed above the waterline (height > 0; 364).
# Both climbing height checks (1 and 3) remove all 0 cm scores; the
# syndrome check (2) removes trials without emergence, because there the
# problem is the censored latency of snails that never left the water.
#   NOT REPORTED (checks of the censoring; set RUN_UNREPORTED <- TRUE)
#   4. Climbing height without censoring: mean effects (lmer),
#      repeatability (rptR), and the syndrome with height not censored.
#
# Run 02 and 05 first: the comparisons read their saved models.
# Each model is fitted once and saved; later runs load the saved file.
# =====================================================================

source("01_setup.R")

h   <- readRDS("outputs/models/height_models.rds")
syn <- readRDS("outputs/models/syndrome_models.rds")

# Fit a model once and save it; on later runs load the saved file.
# (R evaluates `fit` only when it is used, so nothing is refitted when
#  the file already exists.)
fit_or_load <- function(file, fit) {
  if (file.exists(file)) return(readRDS(file))
  saveRDS(fit, file)
  fit
}

dat_climb <- subset(total_snails, max_height_aw > 0)   # climbed above the waterline
dat_em    <- subset(total_snails, emerged == 1)         # emerged at least once
stopifnot(nrow(dat_climb) == 364, nrow(dat_em) == 406)


# =====================================================================
# 1. TABLE S2 — climbing height, mean effects, height > 0 (364 trials)
# Same model as MAIN in 02 (m_height_int). Conditions on climbing, which
# itself depends on treatment and infection, so these are conditional
# rather than total effects. reverse = TRUE gives INV - CON, NAT - CON,
# NAT - INV and Infected - Uninfected (the directions in the tables).
# =====================================================================

m_sens_pos <- fit_or_load("outputs/models/sens_height_pos.rds",
  brm(max_height_aw | cens(h_cens) ~ Treatment * Infection.status + (1 | Snail_ID),
      family  = gaussian(),
      data    = dat_climb,
      chains  = CHAINS, cores = 4, iter = ITER, warmup = WARMUP,
      control = list(adapt_delta = 0.95), seed = 123))

rstan::check_hmc_diagnostics(m_sens_pos$fit)
pairs(emmeans(m_sens_pos, ~ Treatment),        reverse = TRUE)
pairs(emmeans(m_sens_pos, ~ Infection.status), reverse = TRUE)


# =====================================================================
# 2. SYNDROME (STRICT) — climbing height x time to emerge, emerged
#    trials only (406). Removes the shared "did it emerge at all" driver:
#    a snail that never left the water has height 0 AND a censored
#    latency. No censoring is left on latency. Compared with MAIN (05).
# =====================================================================

f_strict <- "outputs/models/sens_syn_strict.rds"
# reuse the earlier fit saved by 05, if present, instead of refitting
if (!file.exists(f_strict) && !is.null(syn$syn_strict)) saveRDS(syn$syn_strict, f_strict)

m_sens_strict <- fit_or_load(f_strict,
  brm(bf(max_height_aw | cens(h_cens) ~ Treatment + Infection.status + (1 |q| Snail_ID)) +
        bf(time_to_emerge_cens ~ Treatment + Infection.status + (1 |q| Snail_ID),
           family = lognormal()) +
        set_rescor(FALSE),
      data    = dat_em,
      chains  = CHAINS, cores = 4, iter = ITER, warmup = WARMUP,
      control = list(adapt_delta = 0.95), seed = 123))

rstan::check_hmc_diagnostics(m_sens_strict$fit)
cor_strict <- get_cor_by(m_sens_strict, "maxheightaw", "timetoemergecens")

summ_draws(syn$cor_full_hcens)   # MAIN, 540 trials (from 05)
summ_draws(cor_strict)           # STRICT, 406 trials


# =====================================================================
# 3. TABLE S3 — climbing height, individual stability, height > 0
#    (364 trials). Same model as MAIN in 02 (fit_stab, defined in 01).
#    Checks whether the 0 cm scores of snails in the water drive the
#    results. Table S3 = MAIN (540) next to this model.
# =====================================================================

m_sens_stab <- fit_or_load("outputs/models/sens_height_stab_pos.rds",
                           fit_stab(dat_climb))
rstan::check_hmc_diagnostics(m_sens_stab$fit)

stab_tab <- function(m) {
  d <- as_draws_df(m)
  rbind("Among-individual SD in individual stability" = summ_draws(d$sd_Snail_ID__sigma_Intercept),
        "Average level x individual stability (r)"    = summ_draws(d$cor_Snail_ID__Intercept__sigma_Intercept),
        "Infected - Uninfected"                       = summ_draws(d$b_sigma_Infection.statusInfected),
        "Invasive - Control"                          = summ_draws(d$b_sigma_TreatmentINV),
        "Native - Control"                            = summ_draws(d$b_sigma_TreatmentNAT))
}

list("All trials (540)"           = stab_tab(h$m_pred),
     "Height above 0 cm (364)"    = stab_tab(m_sens_stab))


# =====================================================================
# 4. NOT REPORTED — climbing height without censoring (27 cm treated as
#    a real value). Checks that the censoring does not drive the results.
# =====================================================================

RUN_UNREPORTED <- FALSE

if (RUN_UNREPORTED) {

  # 4a. Mean effects (lmer)
  model_height_old <- lmer(max_height_aw ~ Treatment * Infection.status +
                             (1 | Snail_ID), data = total_snails)
  print(summary(model_height_old))
  print(anova(model_height_old, type = 2))

  # 4b. Repeatability (rptR)
  rpt_height_old <- rpt(max_height_aw ~ Treatment + Infection.status + (1 | Snail_ID),
                        grname = "Snail_ID", data = total_snails,
                        datatype = "Gaussian", nboot = 1000, npermut = 0)
  print(summary(rpt_height_old))

  # 4c. Syndrome with height not censored
  m_syn_nocens <- fit_or_load("outputs/models/sens_syn_nocens.rds",
    brm(bf(max_height_aw ~ Treatment + Infection.status + (1 |q| Snail_ID)) +
          bf(time_to_emerge_cens | cens(censored) ~ Treatment + Infection.status +
               (1 |q| Snail_ID), family = lognormal()) +
          set_rescor(FALSE),
        data    = total_snails,
        chains  = CHAINS, cores = 4, iter = ITER, warmup = WARMUP,
        control = list(adapt_delta = 0.95), seed = 123))
  print(summ_draws(get_cor_by(m_syn_nocens, "maxheightaw", "timetoemergecens")))
}

