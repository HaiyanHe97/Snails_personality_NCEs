# =====================================================================
# 05  BEHAVIOURAL SYNDROME  (climbing height = height above waterline)
# ---------------------------------------------------------------------
# Do two behaviours covary among individuals? Fitted with brms bivariate
# models. Treatment is a fixed effect (a within-individual manipulation),
# so each correlation is an AMONG-INDIVIDUAL correlation estimated after
# treatment is accounted for. The |q| syntax links the two responses
# into one covariance structure.
#
#   PART A  Syndrome 1: climbing height x TIME TO EMERGE
#   PART B  Syndrome 2: climbing height x EMERGENCE
#
# Climbing height is now height ABOVE the waterline (max_height_aw),
# right-censored at 27 cm (ceiling = 30 cm from base, minus 3 cm water).
#
# Together vs separate: same rule as in 02 (see the header there).
# Does a syndrome exist? -> ONE model (all snails). Does it differ between
# infected and uninfected? -> SEPARATE models, compared by the posterior
# of the difference.
# =====================================================================

source("01_setup.R")

CTRL <- list(adapt_delta = 0.95)

# ---- helper: pretty-print a correlation posterior -------------------
# (define here in case it is not sourced elsewhere)
if (!exists("report_cor")) {
  report_cor <- function(x, label) {
    cat("\n", label, "\n",
        "  r    = ", round(median(x), 3), "\n",
        "  95%  = [", round(quantile(x, 0.025), 3), ", ",
        round(quantile(x, 0.975), 3), "]\n",
        "  P(r<0) = ", round(mean(x < 0), 3),
        " | P(r>0) = ", round(mean(x > 0), 3), "\n", sep = "")
  }
}

# ---- helper: grab a correlation column by matching both responses ---
# Robust to however brms spells the parameter name (handles underscores).
get_cor_by <- function(model, resp_a, resp_b) {
  d  <- as_draws_df(model)
  nm <- grep("^cor_Snail_ID__", names(d), value = TRUE)
  hit <- nm[grepl(resp_a, nm) & grepl(resp_b, nm)]
  if (length(hit) != 1)
    stop("Could not uniquely find cor column. Candidates:\n",
         paste(nm, collapse = "\n"))
  d[[hit]]
}


# =====================================================================
# PART A — SYNDROME 1: climbing height x TIME TO EMERGE
# ---------------------------------------------------------------------
# Do snails that climb high also emerge fast? `emerged` is NOT used
# here: it is the same event as time_to_emerge (emerged == 0 <=> latency
# > 20 min); the censored latency already carries that information.
#
# Sign: NEGATIVE correlation = climbs high AND emerges fast = bold.
#
# TWO VERSIONS, and the comparison is itself a result:
#   FULL   (540 trials, censored) — includes non-emergers.
#   STRICT (406 trials, emerged == 1) — removes the shared driver.
# =====================================================================

# ---- FULL, MAIN: climbing height right-censored at 27 cm   [ONE model]
bf_h2 <- bf(max_height_aw | cens(h_cens) ~ Treatment + Infection.status + (1 |q| Snail_ID))
bf_t2 <- bf(time_to_emerge_cens | cens(censored) ~ Treatment + Infection.status +
              (1 |q| Snail_ID),
            family = lognormal())

m_syn_full_hcens <- brm(
  bf_h2 + bf_t2 + set_rescor(FALSE),
  data    = total_snails,
  chains  = CHAINS, cores = 4, iter = ITER, warmup = WARMUP,
  control = CTRL, seed = 123
)
summary(m_syn_full_hcens)
pp_check(m_syn_full_hcens, resp = "maxheightaw")

cor_full_hcens <- get_cor_by(m_syn_full_hcens, "maxheightaw", "timetoemergecens")
report_cor(cor_full_hcens, "FULL (540 trials, height censored)")

# ---- FULL, SENSITIVITY: height NOT censored ------------------------
bf_h <- bf(max_height_aw ~ Treatment + Infection.status + (1 |q| Snail_ID))
bf_t <- bf(time_to_emerge_cens | cens(censored) ~ Treatment + Infection.status +
             (1 |q| Snail_ID),
           family = lognormal())

m_syn_full <- brm(
  bf_h + bf_t + set_rescor(FALSE),
  data   = total_snails,
  chains = CHAINS, cores = 4, iter = ITER, warmup = WARMUP,
  control = CTRL, seed = 123
)
cor_full <- get_cor_by(m_syn_full, "maxheightaw", "timetoemergecens")
report_cor(cor_full, "FULL (540 trials, height NOT censored — sensitivity)")

# ---- STRICT: emerged == 1 only (no censoring left on latency)   [ONE model]
# Censoring on height uses the ORIGINAL ceiling hit (max_height == 30);
# on the max_height_aw scale this corresponds to 27 cm.
dat_em <- subset(total_snails, emerged == 1)
table(dat_em$h_cens)

bf_h_s <- bf(max_height_aw | cens(h_cens) ~ Treatment + Infection.status + (1 |q| Snail_ID))
bf_t_s <- bf(time_to_emerge_cens ~ Treatment + Infection.status + (1 |q| Snail_ID),
             family = lognormal())

m_syn_strict <- brm(
  bf_h_s + bf_t_s + set_rescor(FALSE),
  data    = dat_em,
  chains  = CHAINS, cores = 4, iter = ITER, warmup = WARMUP,
  control = CTRL, seed = 123
)
summary(m_syn_strict)
pp_check(m_syn_strict, resp = "maxheightaw")

cor_strict <- get_cor_by(m_syn_strict, "maxheightaw", "timetoemergecens")
report_cor(cor_strict, "STRICT (406 emerged trials)")

# How much of the syndrome is the shared "did it emerge at all" driver?
cor_shrink <- cor_full_hcens - cor_strict
cat("FULL - STRICT (median): ", round(median(cor_shrink), 3), "\n")

# ---- By infection status (FULL structure)   [SEPARATE: compare correlations]
bf_h_g <- bf(max_height_aw | cens(h_cens) ~ Treatment + (1 |q| Snail_ID))
bf_t_g <- bf(time_to_emerge_cens | cens(censored) ~ Treatment + (1 |q| Snail_ID),
             family = lognormal())

fit_syn_group <- function(status) {
  brm(bf_h_g + bf_t_g + set_rescor(FALSE),
      data = subset(total_snails, Infection.status == status),
      chains = CHAINS, cores = 4, iter = ITER, warmup = WARMUP,
      control = CTRL, seed = 123)
}
m_syn_inf   <- fit_syn_group("Infected")
m_syn_uninf <- fit_syn_group("Uninfected")
rstan::check_hmc_diagnostics(m_syn_inf$fit)    # 30 snails per group: check
rstan::check_hmc_diagnostics(m_syn_uninf$fit)

cor_ht_inf   <- get_cor_by(m_syn_inf,   "maxheightaw", "timetoemergecens")
cor_ht_uninf <- get_cor_by(m_syn_uninf, "maxheightaw", "timetoemergecens")
report_cor(cor_ht_inf,   "height x tte — Infected (n = 30)")
report_cor(cor_ht_uninf, "height x tte — Uninfected (n = 30)")

cor_ht_diff <- cor_ht_inf - cor_ht_uninf
cat("\nheight x tte, difference (Infected - Uninfected): ",
    round(median(cor_ht_diff), 3),
    " [", round(quantile(cor_ht_diff, 0.025), 3), ", ",
    round(quantile(cor_ht_diff, 0.975), 3), "]\n",
    "P(stronger in infected) = ", round(mean(cor_ht_diff < 0), 3), "\n", sep = "")


# =====================================================================
# PART B — SYNDROME 2: climbing height x EMERGENCE
# ---------------------------------------------------------------------
# Height (gaussian, right-censored at 27 cm) with emergence (bernoulli)
# as the second response. Among-individual correlation of interest.
# =====================================================================

# ---- Overall (all snails)   [ONE model] -----------------------------
bf_height  <- bf(max_height_aw | cens(h_cens) ~ Treatment + Infection.status + (1 |q| Snail_ID))
bf_emerged <- bf(emerged    ~ Treatment + Infection.status + (1 |q| Snail_ID),
                 family = bernoulli())
model_brms_all <- brm(
  bf_height + bf_emerged + set_rescor(FALSE),
  data    = total_snails,
  chains  = CHAINS, cores = 4, iter = ITER, warmup = WARMUP,
  control = CTRL, seed = 123
)
summary(model_brms_all)
rhat(model_brms_all)
pp_check(model_brms_all, resp = "maxheightaw")
pp_check(model_brms_all, resp = "emerged")

cor_he_all <- get_cor_by(model_brms_all, "maxheightaw", "emerged")
cat("height x emerged — overall: r =", round(median(cor_he_all), 3),
    " 95% [", round(quantile(cor_he_all, 0.025), 3), ", ",
    round(quantile(cor_he_all, 0.975), 3), "]\n", sep = "")

# ---- By infection status   [SEPARATE: compare correlations] --------
fit_he_group <- function(status) {
  brm(bf(max_height_aw | cens(h_cens) ~ Treatment + (1 |q| Snail_ID)) +
        bf(emerged  ~ Treatment + (1 |q| Snail_ID), family = bernoulli()) +
        set_rescor(FALSE),
      data = subset(total_snails, Infection.status == status),
      chains = CHAINS, cores = 4, iter = ITER, warmup = WARMUP,
      control = CTRL, seed = 123)
}
model_brms_inf   <- fit_he_group("Infected")
model_brms_uninf <- fit_he_group("Uninfected")
rstan::check_hmc_diagnostics(model_brms_inf$fit)
rstan::check_hmc_diagnostics(model_brms_uninf$fit)

cor_he_inf   <- get_cor_by(model_brms_inf,   "maxheightaw", "emerged")
cor_he_uninf <- get_cor_by(model_brms_uninf, "maxheightaw", "emerged")
cat("height x emerged — Infected:   r =", round(median(cor_he_inf), 3), "\n")
cat("height x emerged — Uninfected: r =", round(median(cor_he_uninf), 3), "\n")

cor_he_diff <- cor_he_inf - cor_he_uninf
cat("Infected   95% [", round(quantile(cor_he_inf,   0.025), 3), ", ",
    round(quantile(cor_he_inf,   0.975), 3), "]\n", sep = "")
cat("Uninfected 95% [", round(quantile(cor_he_uninf, 0.025), 3), ", ",
    round(quantile(cor_he_uninf, 0.975), 3), "]\n", sep = "")
cat("difference (Infected - Uninfected): ", round(median(cor_he_diff), 3),
    " 95% [", round(quantile(cor_he_diff, 0.025), 3), ", ",
    round(quantile(cor_he_diff, 0.975), 3), "]\n",
    "P(infected > uninfected) = ", round(mean(cor_he_diff > 0), 3), "\n", sep = "")


# =====================================================================
# SAVE FITTED MODELS
# =====================================================================
saveRDS(
  list(
    syn_full_hcens = m_syn_full_hcens,
    syn_strict = m_syn_strict,
    syn_inf = m_syn_inf, syn_uninf = m_syn_uninf,
    cor_full_hcens = cor_full_hcens, cor_strict = cor_strict,
    cor_ht_inf = cor_ht_inf, cor_ht_uninf = cor_ht_uninf,
    he_all = model_brms_all, he_inf = model_brms_inf, he_uninf = model_brms_uninf,
    cor_he_inf = cor_he_inf, cor_he_uninf = cor_he_uninf
  ),
  "outputs/models/syndrome_models.rds"
)
