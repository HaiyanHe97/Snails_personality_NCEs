# ============================================================
# 06  FIGURES
# ------------------------------------------------------------
# All manuscript figures, built from the CENSORED models saved
# by 02-05. Climbing height = height above the waterline
# (max_height_aw), right-censored at 27 cm.
# Numbering follows the manuscript (Figure 1 = design, not here):
#   Figure 2  Emergence: proportion emerged
#   Figure 3  Time to emerge: Kaplan-Meier curves
#   Figure 4  Climbing height: per-snail means (boxplot)
#   Figure 5  Individual differences: BLUPs (a height, b time)
#   Figure 6  Reaction norms (a height, b time)
#   Figure S  Appendix: climbing height of snails that climbed
# Global font size = 15 (except Figure 5 ID labels, kept small
# to avoid overlap).
# Output: outputs/figures/
# ============================================================

source("01_setup.R")

h   <- readRDS("outputs/models/height_models.rds")
tte <- readRDS("outputs/models/tte_models.rds")

library(survival)   # Figure 3

# ---- GLOBAL THEME: font size 15 everywhere -----------------
theme_set(
  theme_classic(base_size = 15) +
    theme(
      axis.text    = element_text(face = "bold", size = 15),
      axis.title   = element_text(face = "bold", size = 15),
      strip.text   = element_text(face = "bold", size = 15),
      legend.title = element_text(face = "bold", size = 15),
      legend.text  = element_text(face = "bold", size = 15),
      plot.title   = element_text(face = "bold", size = 15, hjust = 0)
    )
)



# ============================================================
# FIGURE 2 — Emergence: proportion emerged, treatment x infection
# colour = infection status (named on the x axis);
# dark = emerged, light = not emerged (one shared legend)
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
                       fill = Infection.status, alpha = emerged)) +
  geom_bar(stat = "identity", width = 0.7) +
  facet_wrap(~ Treatment) +
  scale_fill_manual(values = c(Infected = col_infected, Uninfected = col_uninfected),
                    guide = "none") +
  scale_alpha_manual(values = c("Emerged" = 1, "Not emerged" = 0.5)) +
  guides(alpha = guide_legend(override.aes = list(fill = "gray30"))) +
  scale_y_continuous(labels = function(x) x * 100) +
  labs(x = NULL, y = "Proportion (%)", alpha = NULL) +
  theme(legend.position = "bottom")
p_emerge
ggsave(file.path(fig_dir, "fig2_emergence_proportion.pdf"),
       p_emerge, width = 10, height = 6)
# ============================================================
# FIGURE 3 — Kaplan-Meier "still in the water" curves
# Descriptive: trials treated as independent; inference is from
# the censored lognormal model in 03.
# ============================================================

km_fit <- survfit(
  Surv(time_to_emerge_cens, emerged) ~ Treatment + Infection.status,
  data = total_snails
)

km_df <- data.frame(
  time   = km_fit$time,
  surv   = km_fit$surv,
  strata = rep(names(km_fit$strata), km_fit$strata)
) %>%
  mutate(
    Treatment = factor(
      trimws(sub(".*Treatment=([^,]+),.*", "\\1", strata)),
      levels = c("CON", "INV", "NAT")
    ),
    Infection = factor(
      trimws(sub(".*Infection.status=(.*)$", "\\1", strata)),
      levels = c("Uninfected", "Infected")
    )
  )

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
  scale_y_continuous(labels = function(x) x * 100, limits = c(0, 1)) +
  scale_x_continuous(limits = c(0, 20), breaks = seq(0, 20, 5)) +
  labs(x = "Time (min)",
       y = "Proportion still in the water (%)",
       colour = "Infection status")
p_km_tte
ggsave(file.path(fig_dir, "fig3_time_to_emerge_km.pdf"),
       p_km_tte, width = 10, height = 5)

# ============================================================
# FIGURE 3 (combined) — Kaplan-Meier curves, all predator cues
# and both infection groups in one panel
# colour = infection status, linetype = predator cue
# ============================================================



km_fit <- survfit(
  Surv(time_to_emerge_cens, emerged) ~ Treatment + Infection.status,
  data = total_snails
)

km_df <- data.frame(
  time   = km_fit$time,
  surv   = km_fit$surv,
  strata = rep(names(km_fit$strata), km_fit$strata)
) %>%
  mutate(
    Treatment = factor(trimws(sub(".*Treatment=([^,]+),.*", "\\1", strata)),
                       levels = c("CON", "INV", "NAT")),
    Infection = factor(trimws(sub(".*Infection.status=(.*)$", "\\1", strata)),
                       levels = c("Uninfected", "Infected"))
  ) %>%
  group_by(Treatment, Infection) %>%
  group_modify(~ add_row(.x, time = 0, surv = 1, .before = 1)) %>%
  ungroup()

p_km_tte_all <- ggplot(km_df,
                       aes(x = time, y = surv,
                           colour = Infection, linetype = Treatment,
                           group = interaction(Infection, Treatment))) +
  geom_step(linewidth = 1.1) +
  scale_colour_manual(values = c("Infected"   = col_infected,
                                 "Uninfected" = col_uninfected)) +
  scale_linetype_manual(values = c(CON = "dotted", INV = "dashed", NAT = "solid")) +
  scale_y_continuous(labels = function(x) x * 100, limits = c(0, 1)) +
  scale_x_continuous(limits = c(0, 20), breaks = seq(0, 20, 5)) +
  labs(x = "Time (min)",
       y = "Proportion still in the water (%)",
       colour = "Infection status", linetype = "Predator cue") +
  theme(legend.key.width = grid::unit(1.5, "cm"))
p_km_tte_all
ggsave(file.path(fig_dir, "fig3b_time_to_emerge_km_combined.pdf"),
       p_km_tte_all, width = 8, height = 6)

# ============================================================
# FIGURE 4 — Climbing height boxplot, RAW per-snail means
#            (height above the waterline; ceiling values shown
#            at their observed 27 cm for this descriptive plot)
# ============================================================

snail_mean_height <- total_snails %>%
  group_by(Snail_ID, Treatment, Infection.status) %>%
  summarise(mean_height = mean(max_height_aw), .groups = "drop")

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
       fill = "Infection status", color = "Infection status")
p_box_raw

ggsave(file.path(fig_dir, "fig4_height_box_raw.pdf"),
       p_box_raw, width = 8, height = 6)


# ============================================================
# FIGURE 5 — Caterpillar (BLUP) plots  [ID labels kept small]
#   (a) Climbing height  (b) Time to emerge (axis reversed)
# BLUPs from m_height_int / m_tte: deviations from the value
# predicted by treatment and infection status, so 0 (dashed)
# = the average snail after accounting for both.
# ============================================================

blups_h <- h$blups
re_t    <- tte$re_tte

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
  theme_classic(base_size = 15) +
  theme(plot.title = element_text(face = "bold", size = 15),
        axis.text.y = element_text(face = "bold", size = 9),   # kept small
        axis.text.x = element_text(face = "bold", size = 15),
        axis.title  = element_text(face = "bold", size = 15),
        legend.position = "none")

p_cat_t <- ggplot(re_t, aes(x = Snail_ID, y = fit, color = infection)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0.3, alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_flip() +
  scale_y_reverse() +
  scale_color_manual(values = c(Infected = col_infected, Uninfected = col_uninfected)) +
  labs(x = NULL, y = "Time to emerge\n(log-scale BLUP, reversed)",
       color = "Infection status", title = "b") +
  theme_classic(base_size = 15) +
  theme(plot.title = element_text(face = "bold", size = 15),
        axis.text.y = element_blank(),
        axis.text.x = element_text(face = "bold", size = 15),
        axis.title  = element_text(face = "bold", size = 15),
        legend.title = element_text(face = "bold", size = 15),
        legend.text  = element_text(face = "bold", size = 15))

p_caterpillar <- p_cat_h + p_cat_t + plot_layout(widths = c(1, 1))
p_caterpillar
ggsave(file.path(fig_dir, "fig5_caterpillar_combined.pdf"),
       p_caterpillar, width = 11, height = 10)




# ============================================================
# FIGURE 6 — Reaction norms: climbing height (a) + time to emerge (b)
#            Model-predicted from the censored random-slope models
#            (h$m_slope, tte$m_slope). Group line = median (a) and
#            median (b) of the individual predictions.
# ============================================================

# --- Panel a: climbing height (censored random-slope model) ---
# Predictions from h$m_slope; the 92 ceiling trials are handled
# as >= 27 (climbed out of the cylinder), not pinned at 27.
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
  summarise(pred = median(pred), .groups = "drop")


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
  geom_hline(yintercept = 27, linetype = "dashed", colour = "gray40") +
  labs(x = NULL, y = "Climbing height (cm)", colour = "Infection status") +
  theme(axis.text.x  = element_blank(),
        axis.ticks.x = element_blank())

# --- Panel b: time to emerge (predict median = exp(mu)) ---
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
  scale_y_reverse(breaks = seq(0, 45, 10)) +
  geom_hline(yintercept = 20, linetype = "dashed", colour = "gray40") +
  labs(x = NULL, y = "Time to emerge (min)", colour = "Infection status") +
  theme(strip.text = element_blank())

p_a_rn <- p_height_rn + labs(title = "a") + theme(legend.position = "none")
p_b_rn <- p_tte_rn + labs(title = "b")

p_reaction_norms <- (p_a_rn / p_b_rn) + plot_layout(guides = "collect")
p_reaction_norms
ggsave(file.path(fig_dir, "fig6_reaction_norms_combined.pdf"),
       p_reaction_norms, width = 8, height = 9, units = "in")


