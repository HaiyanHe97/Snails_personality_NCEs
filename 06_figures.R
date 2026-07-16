# ============================================================
# 06  FIGURES
# ------------------------------------------------------------
# All manuscript figures, built from the CENSORED models saved
# by 02-05 (so figures match the tables). Reads fitted models
# from outputs/models/*.rds; nothing slow is refitted here.
# Output: outputs/figures/
#
# Data objects available:
#   h   <- readRDS(".../height_models.rds")   # m_add, m_slope, blups, ...
#   tte <- readRDS(".../tte_models.rds")       # m_tte, m_slope, re_tte, ind_wide, ...
#   syn <- readRDS(".../syndrome_models.rds")  # cor_full_hcens, cor_strict, cor_ht_inf/uninf, ...
# ============================================================

source("01_setup.R")

h   <- readRDS("outputs/models/height_models.rds")
tte <- readRDS("outputs/models/tte_models.rds")
syn <- readRDS("outputs/models/syndrome_models.rds")

library(tidybayes)   # for spread_draws / fitted helpers if needed


# ============================================================
# FIG 1a — Climbing height boxplot, RAW per-snail means
#          (treatment x infection). Descriptive version:
#          treats the 92 ceiling values as observed 30 cm.
# ============================================================

snail_mean_height <- total_snails %>%
  group_by(Snail_ID, Treatment, Infection.status) %>%
  summarise(mean_height = mean(max_height), .groups = "drop")

group_means_height <- snail_mean_height %>%
  group_by(Treatment, Infection.status) %>%
  summarise(mean_height = mean(mean_height), .groups = "drop")

dodge <- position_dodge(0.8)

p_box_raw <- ggplot(snail_mean_height,
                    aes(x = Treatment, y = mean_height, fill = Infection.status)) +
  geom_boxplot(position = dodge, width = 0.7, alpha = 0.65, outlier.shape = NA) +
  geom_point(aes(color = Infection.status),
             position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8),
             alpha = 0.5, size = 2, shape = 16) +
  scale_fill_manual(values  = c(Infected = col_infected, Uninfected = col_uninfected)) +
  scale_color_manual(values = c(Infected = col_infected, Uninfected = col_uninfected)) +
  labs(x = NULL, y = "Climbing height (cm)",
       fill = "Infection status", color = "Infection status") +
  theme_classic() +
  theme(axis.text  = element_text(face = "bold", size = 14),
        axis.title = element_text(face = "bold", size = 16),
        legend.title = element_text(face = "bold", size = 13),
        legend.text  = element_text(face = "bold", size = 13))
p_box_raw

ggsave(file.path(fig_dir, "fig1a_height_box_raw.pdf"),
       p_box_raw, width = 8, height = 6)

# ============================================================
# FIG 1b — Climbing height, MODEL-PREDICTED (censored)
#          Predictions from m_slope (censored random-slope),
#          so the 92 ceiling values are handled as >= 30, not
#          pinned at 30. This is the version consistent with
#          the tables.
# ============================================================

newdat_h <- expand.grid(
  Snail_ID  = levels(total_snails$Snail_ID),
  Treatment = levels(total_snails$Treatment)
) %>%
  mutate(Infection.status = total_snails$Infection.status[
    match(Snail_ID, total_snails$Snail_ID)])

newdat_h$pred <- fitted(h$m_slope, newdata = newdat_h,
                        re_formula = ~ (Treatment | Snail_ID),
                        scale = "response")[, "Estimate"]

group_pred_h <- newdat_h %>%
  group_by(Treatment, Infection.status) %>%
  summarise(pred = mean(pred), .groups = "drop")

p_box_pred <- ggplot() +
  geom_line(data = newdat_h,
            aes(x = Treatment, y = pred, group = Snail_ID, color = Infection.status),
            linewidth = 0.45, alpha = 0.18) +
  geom_point(data = newdat_h,
             aes(x = Treatment, y = pred, color = Infection.status),
             size = 1.3, alpha = 0.25) +
  geom_line(data = group_pred_h,
            aes(x = Treatment, y = pred, group = Infection.status, color = Infection.status),
            linewidth = 1.8) +
  geom_point(data = group_pred_h,
             aes(x = Treatment, y = pred, color = Infection.status),
             shape = 18, size = 4) +
  scale_color_manual(values = c(Infected = col_infected, Uninfected = col_uninfected)) +
  labs(x = NULL, y = "Climbing height (cm), model-predicted",
       color = "Infection status") +
  theme_classic() +
  theme(axis.text  = element_text(face = "bold", size = 14),
        axis.title = element_text(face = "bold", size = 16),
        legend.title = element_text(face = "bold", size = 13),
        legend.text  = element_text(face = "bold", size = 13))
p_box_pred
ggsave(file.path(fig_dir, "fig1b_height_predicted.pdf"),
       p_box_pred, width = 8, height = 6)


# ============================================================
# FIG 2 — Emergence: proportion emerged, by treatment x infection
# ============================================================

emerged_prop <- total_snails %>%
  group_by(Treatment, Infection.status, emerged) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(Treatment, Infection.status) %>%
  mutate(prop = n / sum(n),
         emerged = factor(emerged, levels = c(1, 0),
                          labels = c("Emerged", "Not emerged")))

p_emerge <- ggplot(emerged_prop,
                   aes(x = Infection.status, y = prop,
                       fill = interaction(emerged, Infection.status))) +
  geom_bar(stat = "identity", width = 0.7) +
  facet_wrap(~ Treatment) +
  scale_fill_manual(
    values = c("Emerged.Infected"       = col_infected,
               "Not emerged.Infected"   = col_infected_light,
               "Emerged.Uninfected"     = col_uninfected,
               "Not emerged.Uninfected" = col_uninfected_light),
    breaks = c("Emerged.Uninfected", "Not emerged.Uninfected"),
    labels = c("Emerged", "Not emerged")) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = NULL, y = "Proportion", fill = NULL) +
  theme_classic() +
  theme(legend.position = "bottom",
        strip.text = element_text(face = "bold", size = 14),
        axis.text  = element_text(face = "bold", size = 13),
        axis.title = element_text(face = "bold", size = 16),
        legend.text = element_text(face = "bold", size = 13))
p_emerge
ggsave(file.path(fig_dir, "fig2_emergence_proportion.pdf"),
       p_emerge, width = 10, height = 6)


# ============================================================
# FIG 3 — Caterpillar (BLUP) plots, side by side
#   (a) Climbing height  : positive = climbs higher = bolder
#   (b) Time to emerge    : axis REVERSED so up = faster = bolder,
#                           matching panel (a)
# Snails ordered by their climbing-height BLUP in BOTH panels so
# the same individual sits at the same y-position.
# ============================================================

blups_h <- h$blups
re_t    <- tte$re_tte

# common ordering by height BLUP
order_ids <- blups_h %>% arrange(fit) %>% pull(Snail_ID)
blups_h$Snail_ID <- factor(blups_h$Snail_ID, levels = order_ids)
re_t$Snail_ID    <- factor(re_t$Snail_ID,    levels = order_ids)

p_cat_h <- ggplot(blups_h, aes(x = Snail_ID, y = fit, color = infection)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0.3, alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_flip() +
  scale_color_manual(values = c(Infected = col_infected, Uninfected = col_uninfected)) +
  labs(x = "Snail ID", y = "Climbing height (cm)",
       color = "Infection status", title = "a") +
  theme_classic() +
  theme(plot.title = element_text(face = "bold", size = 16),
        axis.text  = element_text(face = "bold", size = 9),
        axis.title = element_text(face = "bold", size = 13),
        legend.position = "none")

p_cat_t <- ggplot(re_t, aes(x = Snail_ID, y = fit, color = infection)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0.3, alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_flip() +
  scale_y_reverse() +                              # up = faster = bolder
  scale_color_manual(values = c(Infected = col_infected, Uninfected = col_uninfected)) +
  labs(x = NULL, y = "Time to emerge(reversed)",
       color = "Infection status", title = "b") +
  theme_classic() +
  theme(plot.title = element_text(face = "bold", size = 16),
        axis.text.y = element_blank(),
        axis.text.x = element_text(face = "bold", size = 9),
        axis.title  = element_text(face = "bold", size = 13),
        legend.title = element_text(face = "bold", size = 12),
        legend.text  = element_text(face = "bold", size = 12))

p_caterpillar <- p_cat_h + p_cat_t + plot_layout(widths = c(1, 1))
p_caterpillar
ggsave(file.path(fig_dir, "fig3_caterpillar_combined.pdf"),
       p_caterpillar, width = 11, height = 10)


# ============================================================
# FIG 4 — Kaplan-Meier "still in the water" curves
#         Faceted by treatment (CON/INV/NAT), two lines per
#         panel (infection status). Right-censoring is visible
#         as curves that do NOT reach 0% by 20 min: that final
#         height is the proportion of trials that never emerged.
#
#         Descriptive (empirical KM), complementing the censored
#         lognormal model in 03.
# ============================================================

library(survival)
library(broom)

# Event = emerged (1); censored (0) = still in the water at 20 min.
# time  = time_to_emerge_cens (emergers: real time; non-emergers: 20).
km_fit <- survfit(
  Surv(time_to_emerge_cens, emerged) ~ Treatment + Infection.status,
  data = total_snails
)

# Tidy into a data frame; "surv" = proportion still in the water.
km_df <- tidy(km_fit) %>%
  # strata looks like "Treatment=CON, Infection.status=Uninfected"
  mutate(
    Treatment = sub(".*Treatment=([^,]+).*", "\\1", strata),
    Infection.status = sub(".*Infection.status=([^,]+).*", "\\1", strata),
    Treatment = factor(Treatment, levels = c("CON", "INV", "NAT")),
    Infection.status = factor(Infection.status, levels = c("Uninfected", "Infected"))
  )

# Add a t = 0 row (100% still in the water) for each group so lines start at the top.
km_start <- km_df %>%
  distinct(Treatment, Infection.status) %>%
  mutate(time = 0, estimate = 1, conf.low = 1, conf.high = 1)

km_df <- data.frame(
  time   = km_fit$time,
  surv   = km_fit$surv,
  strata = rep(names(km_fit$strata), km_fit$strata)
) %>%
  mutate(
    # survfit pads its strata names to equal width, so "Infected  " comes
    # with trailing spaces and never matches scale_colour_manual. trimws()
    # is what fixes the colours.
    Treatment = factor(
      trimws(sub(".*Treatment=([^,]+),.*", "\\1", strata)),
      levels = c("CON", "INV", "NAT")
    ),
    Infection = factor(
      trimws(sub(".*Infection.status=(.*)$", "\\1", strata)),
      levels = c("Uninfected", "Infected")
    )
  )
# every curve should start at 1 (all snails in the water at t = 0)
km_df <- km_df %>%
  group_by(Treatment, Infection) %>%
  group_modify(~ add_row(.x, time = 0, surv = 1, .before = 1)) %>%
  ungroup()

stopifnot(!any(is.na(km_df$Treatment)), !any(is.na(km_df$Infection)))

p_km_tte <- ggplot(km_df, aes(x = time, y = surv, colour = Infection)) +
  geom_step(linewidth = 1.1) +
  facet_wrap(~ Treatment) +
  scale_colour_manual(values = c("Infected"   = col_infected,
                                 "Uninfected" = col_uninfected)) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  scale_x_continuous(limits = c(0, 20), breaks = seq(0, 20, 5)) +
  labs(x = "Time (min)",
       y = "Proportion still in the water",
       colour = "Infection status") +
  theme_classic() +
  theme(strip.text   = element_text(face = "bold", size = 16),
        axis.text    = element_text(face = "bold", size = 14),
        axis.title   = element_text(face = "bold", size = 18),
        legend.title = element_text(face = "bold", size = 14),
        legend.text  = element_text(face = "bold", size = 16))

p_km_tte


ggsave(file.path(fig_dir, "fig4_time_to_emerge_km.pdf"),
       p_km_tte, width = 10, height = 5)

# ============================================================
# FIG 5 — Reaction norms: climbing height (a) + time to emerge (b)
#         Both panels use MODEL-PREDICTED values (censored models),
#         not raw per-snail means, so the figure is consistent with
#         the tables and doesn't distort by the ceiling/censoring
#         issues described below.
# ============================================================

# --------------------------------------------------------------
# Panel a — Climbing height (from the censored random-slope model)
#
# Predictions come from h$m_slope, not from raw per-snail means.
# Raw means treat the 92 ceiling trials as exactly 30 cm, which is
# the very thing the censored model corrects: those snails climbed
# OUT of the cylinder, so their true height is >= 30 and unobserved.
# Using raw means here would show a flattened figure that
# contradicts the model.
# --------------------------------------------------------------

newdat_h_rn <- expand.grid(
  Snail_ID  = levels(total_snails$Snail_ID),
  Treatment = levels(total_snails$Treatment)
) %>%
  mutate(Infection.status = total_snails$Infection.status[
    match(Snail_ID, total_snails$Snail_ID)])

newdat_h_rn$pred <- fitted(h$m_slope, newdata = newdat_h_rn,
                           re_formula = ~ (Treatment | Snail_ID),
                           scale = "response")[, "Estimate"]

group_pred_h_rn <- newdat_h_rn %>%
  group_by(Treatment, Infection.status) %>%
  summarise(pred = mean(pred), .groups = "drop")

p_height_rn <- ggplot() +
  geom_line(data = newdat_h_rn,
            aes(x = Treatment, y = pred, group = Snail_ID,
                colour = Infection.status),
            linewidth = 0.45, alpha = 0.18) +
  geom_point(data = newdat_h_rn,
             aes(x = Treatment, y = pred, colour = Infection.status),
             size = 1.3, alpha = 0.25) +
  geom_line(data = group_pred_h_rn,
            aes(x = Treatment, y = pred, group = Infection.status,
                colour = Infection.status),
            linewidth = 1.8) +
  geom_point(data = group_pred_h_rn,
             aes(x = Treatment, y = pred, colour = Infection.status),
             shape = 18, size = 4) +
  facet_wrap(~ Infection.status) +
  scale_colour_manual(values = c("Infected"   = col_infected,
                                 "Uninfected" = col_uninfected)) +
  geom_hline(yintercept = 30, linetype = "dashed", colour = "gray40") +
  labs(x = NULL, y = "Climbing height (cm)", colour = "Infection status") +
  theme_classic() +
  theme(strip.text   = element_text(face = "bold", size = 15),
        axis.text.x  = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y  = element_text(face = "bold", size = 13),
        axis.title   = element_text(face = "bold", size = 15),
        legend.title = element_text(face = "bold", size = 13),
        legend.text  = element_text(face = "bold", size = 13))

# --------------------------------------------------------------
# Panel b — Time to emerge
#
# Predictions come from tte$m_slope, not from raw per-snail means.
# Raw means would silently drop the 134 trials where the snail
# never emerged, and those are not random (98 infected vs 36
# uninfected).
#
# We predict the MEDIAN, i.e. exp(mu). Do NOT use scale =
# "response": for a lognormal, that returns the MEAN,
# exp(mu + sigma^2/2), which with sigma = 0.81 inflates every value
# by ~39% and pushes predictions past the 20-min limit of the trial.
# --------------------------------------------------------------

newdat_tte_rn <- expand.grid(
  Snail_ID  = levels(total_snails$Snail_ID),
  Treatment = levels(total_snails$Treatment)
) %>%
  mutate(Infection.status = total_snails$Infection.status[
    match(Snail_ID, total_snails$Snail_ID)])

newdat_tte_rn$pred <- exp(
  fitted(tte$m_slope, newdata = newdat_tte_rn,
         re_formula = ~ (Treatment | Snail_ID),
         scale = "linear")[, "Estimate"]
)

group_pred_tte_rn <- newdat_tte_rn %>%
  group_by(Treatment, Infection.status) %>%
  summarise(pred = median(pred), .groups = "drop")

p_tte_rn <- ggplot() +
  geom_line(data = newdat_tte_rn,
            aes(x = Treatment, y = pred, group = Snail_ID,
                colour = Infection.status),
            linewidth = 0.45, alpha = 0.2) +
  geom_point(data = newdat_tte_rn,
             aes(x = Treatment, y = pred, colour = Infection.status),
             size = 1.3, alpha = 0.25) +
  geom_line(data = group_pred_tte_rn,
            aes(x = Treatment, y = pred, group = Infection.status,
                colour = Infection.status),
            linewidth = 1.8) +
  geom_point(data = group_pred_tte_rn,
             aes(x = Treatment, y = pred, colour = Infection.status),
             shape = 18, size = 4) +
  facet_wrap(~ Infection.status) +
  scale_colour_manual(values = c("Infected"   = col_infected,
                                 "Uninfected" = col_uninfected)) +
  # reversed and capped at the 20-min trial length, so that UP =
  # faster = bolder, matching panel a
  scale_y_reverse(breaks = seq(0, 45, 10)) +
  geom_hline(yintercept = 20, linetype = "dashed", colour = "gray40") +
  labs(x = NULL, y = "Time to emerge (min)", colour = "Infection status") +
  theme_classic() +
  theme(strip.text   = element_blank(),      # facet labels already on panel a
        axis.text    = element_text(face = "bold", size = 13),
        axis.title   = element_text(face = "bold", size = 15),
        legend.title = element_text(face = "bold", size = 13),
        legend.text  = element_text(face = "bold", size = 13))

# --------------------------------------------------------------
# Combine panels a + b (needs `patchwork`, already used above for
# FIG 3's caterpillar plots)
# --------------------------------------------------------------

p_a_rn <- p_height_rn +
  labs(title = "a") +
  theme(plot.title = element_text(face = "bold", size = 18, hjust = 0),
        legend.position = "none")

p_b_rn <- p_tte_rn +
  labs(title = "b") +
  theme(plot.title = element_text(face = "bold", size = 18, hjust = 0))

p_reaction_norms <- (p_a_rn / p_b_rn) + plot_layout(guides = "collect")
p_reaction_norms

ggsave(file.path(fig_dir, "fig5_reaction_norms_combined.pdf"),
       p_reaction_norms, width = 8, height = 9, units = "in")
