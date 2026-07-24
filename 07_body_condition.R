# =====================================================================
# 07  BODY CONDITION — is the infection effect confounded by size?
# ---------------------------------------------------------------------
# Infected and uninfected snails were collected from the same site on the
# same day, so any difference in body condition would be a confound for
# the behavioural comparisons rather than a result in its own right.
#
# Shell dimensions (height, aperture, width) are strongly intercorrelated,
# so they are reduced to a single size axis with a PCA (PC1). Body
# condition is then weight at a given shell size: the Infection.status
# coefficient in
#
#     log(weight) ~ shell_PC1 + Infection.status
#
# tests whether infected snails are lighter or heavier than uninfected
# snails of the same shell size.
# =====================================================================

source("01_setup.R")


# ---- 1. One row per snail -------------------------------------------
# Shell measurements and weight are already in total_snails: one value per
# snail, repeated across that snail's nine trials. Decimal separators and
# factor levels were handled in 00_merge_clean.R and 01_setup.R, so nothing
# is re-read or re-coerced here.

morph <- total_snails %>%
  distinct(Snail_ID, Weight_g, Height_cm, Aperture_cm, Width_cm, Infection.status)

stopifnot(
  nrow(morph) == 60,
  n_distinct(morph$Snail_ID) == 60,
  !any(is.na(morph[, c("Weight_g", "Height_cm", "Aperture_cm", "Width_cm")])),
  all(levels(morph$Infection.status) == c("Uninfected", "Infected")),
  all(table(morph$Infection.status) == 30)
)


# ---- 2. Shell size axis (PC1) ---------------------------------------
# The sign of a principal component is arbitrary. PC1 is oriented here so
# that larger values always mean a larger shell, otherwise the sign of the
# regression coefficient would depend on the run.

shell_pca <- prcomp(morph %>% select(Height_cm, Aperture_cm, Width_cm),
                    center = TRUE, scale. = TRUE)

pc1_sign <- sign(shell_pca$rotation["Height_cm", 1])
shell_pca$rotation[, 1] <- shell_pca$rotation[, 1] * pc1_sign
shell_pca$x[, 1]        <- shell_pca$x[, 1]        * pc1_sign

morph <- morph %>%
  mutate(shell_PC1  = shell_pca$x[, 1],
         log_weight = log(Weight_g))

pc1_var <- summary(shell_pca)$importance[2, 1] * 100

cat("PC1 loadings:\n"); print(round(shell_pca$rotation[, 1], 3))
cat(sprintf("PC1 explains %.1f%% of shell shape variance\n", pc1_var))
stopifnot(all(shell_pca$rotation[, 1] > 0))


# ---- 3. Condition model ---------------------------------------------
# MAIN: the Infection.status coefficient is the condition test.
# The shell_PC1 coefficient and the R2 quoted in the text both come from
# this model, so that a single model is reported throughout.

m_condition <- lm(log_weight ~ shell_PC1 + Infection.status, data = morph)

summary(m_condition)
confint(m_condition)

cat("\ncondition check (log weight ~ shell size + infection status)\n")
print(round(cbind(Estimate = coef(m_condition),
                  confint(m_condition)), 4))
cat(sprintf("  R2 = %.3f (adjusted %.3f)\n",
            summary(m_condition)$r.squared,
            summary(m_condition)$adj.r.squared))

# Size-only model, for comparison with the R2 above.
m_size_only <- lm(log_weight ~ shell_PC1, data = morph)
cat(sprintf("  R2, shell size alone = %.3f\n", summary(m_size_only)$r.squared))

# Model assumptions.
par(mfrow = c(2, 2)); plot(m_condition); par(mfrow = c(1, 1))


# ---- 4. Figure S1a: what the size axis represents --------------------

loadings_df <- data.frame(dimension = rownames(shell_pca$rotation),
                          loading   = shell_pca$rotation[, 1])

plot_pca_loadings <- ggplot(loadings_df, aes(x = dimension, y = loading)) +
  geom_col(fill = col_uninfected) +
  labs(subtitle = sprintf("PC1 explains %.1f%% of variance", pc1_var),
       x = NULL, y = "Loading on PC1") +
  theme_classic() +
  theme(axis.text  = element_text(face = "bold", size = 14),
        axis.title = element_text(face = "bold", size = 18))

print(plot_pca_loadings)

ggsave(file.path(fig_dir, "FigS1a_pca_loadings.pdf"),
       plot = plot_pca_loadings, width = 6, height = 5, units = "in")


# ---- 5. Figure S1b: weight against shell size ------------------------
# The fitted line is the size-only model, shown for display. The two
# infection groups overlap, consistent with the condition test above.

pred_con <- data.frame(
  shell_PC1 = seq(min(morph$shell_PC1), max(morph$shell_PC1), length.out = 100)
)
pred_con <- cbind(pred_con,
                  predict(m_size_only, newdata = pred_con, interval = "confidence"))

plot_condition <- ggplot(morph, aes(x = shell_PC1, y = log_weight)) +
  geom_ribbon(data = pred_con, aes(x = shell_PC1, ymin = lwr, ymax = upr),
              inherit.aes = FALSE, fill = col_con, alpha = 0.15) +
  geom_line(data = pred_con, aes(x = shell_PC1, y = fit),
            inherit.aes = FALSE, colour = col_con, linewidth = 1.3) +
  geom_point(aes(colour = Infection.status), size = 2, alpha = 0.8) +
  scale_colour_manual(values = c("Infected"   = col_infected,
                                 "Uninfected" = col_uninfected)) +
  labs(x = "Shell size (PC1)", y = "log(Weight)", colour = "Infection status") +
  theme_classic() +
  theme(axis.text    = element_text(face = "bold", size = 14),
        axis.title   = element_text(face = "bold", size = 18),
        legend.title = element_text(face = "bold", size = 14),
        legend.text  = element_text(face = "bold", size = 16))

print(plot_condition)

ggsave(file.path(fig_dir, "FigS1b_condition_check.pdf"),
       plot = plot_condition, width = 8, height = 6, units = "in")


# ---- 6. Save ---------------------------------------------------------

saveRDS(
  list(morph = morph, pca = shell_pca, m_condition = m_condition,
       m_size_only = m_size_only),
  "outputs/models/body_condition.rds"
)
