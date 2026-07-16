# =====================================================================
# 05  BEHAVIOURAL SYNDROME
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
# (Reaction-norm covariance and predictability are in 06_covariance.R.)
#
# Figures live in 07_figures.R. This script fits the models, prints the
# numerical results, and saves them to outputs/models/syndrome_models.rds.
#
# WARNING: several multivariate brms models, a few minutes each.
# =====================================================================

source("01_setup.R")

CTRL <- list(adapt_delta = 0.95)


# =====================================================================
# PART A — SYNDROME 1: climbing height x TIME TO EMERGE
# ---------------------------------------------------------------------
# Do snails that climb high also emerge fast? `emerged` is NOT used
# here: it is the same event as time_to_emerge (emerged == 0 <=> latency
# > 20 min); the censored latency already carries that information.
#
# Sign: NEGATIVE correlation = climbs high AND emerges fast = bold
# (latency runs opposite to height, hence the minus sign).
#
# TWO VERSIONS, and the comparison is itself a result:
#   FULL   (540 trials, censored) — includes non-emergers. Part of the
#          correlation is structural (a snail that never left the water
#          has height = 0 AND a censored latency).
#   STRICT (406 trials, emerged == 1) — removes the shared driver and
#          tests whether the two behaviours are genuinely linked.
# =====================================================================

COR_HT <- "cor_Snail_ID__maxheight_Intercept__timetoemergecens_Intercept"
get_cor_ht <- function(model) as_draws_df(model)[[COR_HT]]

# ---- FULL, MAIN: climbing height ALSO right-censored at 30 ----------
# (chosen version: height censoring kept in the syndrome model)
bf_h2 <- bf(max_height | cens(h_cens) ~ Treatment + Infection.status + (1 |q| Snail_ID))
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
pp_check(m_syn_full_hcens, resp = "maxheight")
cor_full_hcens <- get_cor_ht(m_syn_full_hcens)
report_cor(cor_full_hcens, "FULL (540 trials, height censored)")

# ---- FULL, SENSITIVITY (OLD): height NOT censored ------------------
# Kept for comparison; height treated as observed at 30.
bf_h <- bf(max_height ~ Treatment + Infection.status + (1 |q| Snail_ID))
bf_t <- bf(time_to_emerge_cens | cens(censored) ~ Treatment + Infection.status +
             (1 |q| Snail_ID),
           family = lognormal())

m_syn_full <- brm(
  bf_h + bf_t + set_rescor(FALSE),
  data   = total_snails,
  chains = CHAINS, cores = 4, iter = ITER, warmup = WARMUP,
  control = CTRL, seed = 123
)
cor_full <- get_cor_ht(m_syn_full)
report_cor(cor_full, "FULL (540 trials, height NOT censored — sensitivity)")

# ---- STRICT: emerged == 1 only (no censoring left on latency) ------
# max_height still has 42 zeros here (crossed the waterline then came
# back down) — real observations, kept in. Height still censored at 30.
dat_em <- subset(total_snails, emerged == 1) %>%
  mutate(h_cens = ifelse(max_height == 30, "right", "none"))
table(dat_em$h_cens)

bf_h_s <- bf(max_height | cens(h_cens) ~ Treatment + Infection.status + (1 |q| Snail_ID))
bf_t_s <- bf(time_to_emerge_cens ~ Treatment + Infection.status + (1 |q| Snail_ID),
             family = lognormal())

m_syn_strict <- brm(
  bf_h_s + bf_t_s + set_rescor(FALSE),
  data    = dat_em,
  chains  = CHAINS, cores = 4, iter = ITER, warmup = WARMUP,
  control = CTRL, seed = 123
)
summary(m_syn_strict)
pp_check(m_syn_strict, resp = "maxheight")
cor_strict <- get_cor_ht(m_syn_strict)
report_cor(cor_strict, "STRICT (406 emerged trials)")

# How much of the syndrome is the shared "did it emerge at all" driver?
cor_shrink <- cor_full_hcens - cor_strict
cat("FULL - STRICT (median): ", round(median(cor_shrink), 3), "\n")

# ---- By infection status (FULL structure) --------------------------
bf_h_g <- bf(max_height | cens(h_cens) ~ Treatment + (1 |q| Snail_ID))
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

cor_ht_inf   <- get_cor_ht(m_syn_inf)
cor_ht_uninf <- get_cor_ht(m_syn_uninf)
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
# Height (gaussian, right-censored at 30) with emergence (bernoulli) as
# the second response. Treatment as fixed effect; among-individual
# correlation of interest. Height censoring kept consistent with the
# rest of the manuscript.
# =====================================================================
COR_HE <- "cor_Snail_ID__maxheight_Intercept__emerged_Intercept"
get_cor_he <- function(model) as_draws_df(model)[[COR_HE]]
# ---- Overall (all snails) ------------------------------------------
bf_height  <- bf(max_height | cens(h_cens) ~ Treatment + Infection.status + (1 |q| Snail_ID))
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
pp_check(model_brms_all, resp = "maxheight")
pp_check(model_brms_all, resp = "emerged")
cor_he_all <- get_cor_he(model_brms_all)
cat("height x emerged — overall: r =", round(mean(cor_he_all), 3),
    " 95% [", round(quantile(cor_he_all, 0.025), 3), ", ",
    round(quantile(cor_he_all, 0.975), 3), "]\n", sep = "")
# ---- By infection status -------------------------------------------
fit_he_group <- function(status) {
  brm(bf(max_height | cens(h_cens) ~ Treatment + (1 |q| Snail_ID)) +
        bf(emerged  ~ Treatment + (1 |q| Snail_ID), family = bernoulli()) +
        set_rescor(FALSE),
      data = subset(total_snails, Infection.status == status),
      chains = CHAINS, cores = 4, iter = ITER, warmup = WARMUP,
      control = CTRL, seed = 123)
}
model_brms_inf   <- fit_he_group("Infected")
model_brms_uninf <- fit_he_group("Uninfected")
cor_he_inf   <- get_cor_he(model_brms_inf)
cor_he_uninf <- get_cor_he(model_brms_uninf)
cat("height x emerged — Infected:   r =", round(mean(cor_he_inf), 3), "\n")
cat("height x emerged — Uninfected: r =", round(mean(cor_he_uninf), 3), "\n")
cor_he_diff <- cor_he_inf - cor_he_uninf

quantile(cor_he_inf,   c(0.025, 0.975))
quantile(cor_he_uninf, c(0.025, 0.975))

cat("difference (Infected - Uninfected): ", round(mean(cor_he_diff), 3),
    " 95% [", round(quantile(cor_he_diff, 0.025), 3), ", ",
    round(quantile(cor_he_diff, 0.975), 3), "]\n",
    "P(infected > uninfected) = ", round(mean(cor_he_diff > 0), 3), "\n", sep = "")


# =====================================================================
# SAVE FITTED MODELS  (for figures in 07)
# ---------------------------------------------------------------------
# outputs/models/*.rds is git-ignored (large, regenerated by this script).
# =====================================================================

saveRDS(
  list(
    # Part A: height x tte
    syn_full_hcens = m_syn_full_hcens,
    syn_strict = m_syn_strict,
    syn_inf = m_syn_inf, syn_uninf = m_syn_uninf,
    cor_full_hcens = cor_full_hcens, cor_strict = cor_strict,
    cor_ht_inf = cor_ht_inf, cor_ht_uninf = cor_ht_uninf,
    # Part B: height x emerged
    he_all = model_brms_all, he_inf = model_brms_inf, he_uninf = model_brms_uninf,
    cor_he_inf = cor_he_inf, cor_he_uninf = cor_he_uninf
  ),
  "outputs/models/syndrome_models.rds"
)
