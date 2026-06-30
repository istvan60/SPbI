# ============================================================
# d-excess diagnostic analyses  [KZ4.1 – KZ7.1 + seasonal]
# Requires: all_imputed in memory (from Stage 0)
# ============================================================

library(dplyr)
library(lubridate)
library(ggplot2)
library(patchwork)

stopifnot(exists("all_imputed"))

# ── Base dataset: non-fallback rows, error columns ───────────
base <- all_imputed %>%
  filter(!fallback) %>%
  mutate(
    month    = month(Date),
    mon_lab  = factor(month.abb[month], levels = month.abb),
    X_pct    = factor(paste0(round(X * 100), "%"),
                      levels = paste0(round(sort(unique(X)) * 100), "%")),
    err_O18  = O18_imp - O18_orig,
    err_H2   = H2_imp  - H2_orig,
    err_dex  = d_ex_imp - d_ex_orig
  )

method_order <- c("LOCF","Linear","Spline","Stine",
                  "Kalman","Moving-average","Sinusoidal","SPbI")
base$method <- factor(base$method, levels = method_order)

# ── [KZ4.1]  Error propagation O18 & H2 → d-excess ─────────
# d-excess error = err_H2 - 8*err_O18  (analytical identity)
# Show how much each component contributes per method

prop <- base %>%
  group_by(method) %>%
  summarise(
    sd_O18_component = sd(8 * err_O18,  na.rm = TRUE),   # 8 * σ(ε_O18)
    sd_H2_component  = sd(err_H2,       na.rm = TRUE),   # σ(ε_H2)
    sd_dex_actual    = sd(err_dex,      na.rm = TRUE),   # empirical
    r_O18_H2         = cor(err_O18, err_H2, use = "complete.obs"),
    .groups = "drop"
  ) %>%
  mutate(
    sd_dex_if_uncorr = sqrt(sd_H2_component^2 + sd_O18_component^2),
    dominance = ifelse(sd_O18_component > sd_H2_component,
                       "8×O18 dominates", "H2 dominates")
  )

message("\n── [KZ4.1] Error propagation ──")
print(as.data.frame(round(prop[,-1], 4)))

# ── [KZ5.1]  Correlation of δ¹⁸O and δ²H errors ────────────
err_cor <- base %>%
  group_by(method, X_pct) %>%
  summarise(r = cor(err_O18, err_H2, use = "complete.obs"), .groups = "drop")

p_cor <- ggplot(err_cor, aes(x = X_pct, y = r, colour = method, group = method)) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 2) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  scale_colour_brewer(palette = "Dark2") +
  labs(
    title   = "[KZ5.1] Pearson r between δ¹⁸O and δ²H imputation errors",
    x       = "Masking fraction",
    y       = "Pearson r (err_O18 vs err_H2)",
    colour  = NULL
  ) +
  theme_minimal(base_size = 11)

print(p_cor)
ggsave(file.path(out_dir, "KZ5.1_error_correlation.pdf"),
       p_cor, width = 7, height = 4)

# ── [KZ6.1]  Sinusoidal amplitude damping in d-excess ───────
sin_dat <- base %>%
  filter(method == "Sinusoidal") %>%
  group_by(mon_lab) %>%
  summarise(
    d_ex_obs = mean(d_ex_orig, na.rm = TRUE),
    d_ex_imp = mean(d_ex_imp,  na.rm = TRUE),
    .groups = "drop"
  ) %>%
  tidyr::pivot_longer(c(d_ex_obs, d_ex_imp),
                      names_to = "source", values_to = "d_excess") %>%
  mutate(source = recode(source,
                         d_ex_obs = "Observed",
                         d_ex_imp = "Sinusoidal imputed"))

amp <- base %>%
  filter(method == "Sinusoidal") %>%
  summarise(
    amp_obs = diff(range(tapply(d_ex_orig, month(Date), mean, na.rm = TRUE))),
    amp_imp = diff(range(tapply(d_ex_imp,  month(Date), mean, na.rm = TRUE)))
  )
amp_label <- sprintf("Amplitude: observed = %.1f‰, imputed = %.1f‰  (retention = %.0f%%)",
                     amp$amp_obs, amp$amp_imp, 100 * amp$amp_imp / amp$amp_obs)

p_sin <- ggplot(sin_dat, aes(x = mon_lab, y = d_excess,
                              colour = source, group = source)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_colour_manual(values = c("Observed" = "black",
                                 "Sinusoidal imputed" = "#E41A1C")) +
  labs(
    title    = "[KZ6.1] Sinusoidal imputation: seasonal d-excess cycle",
    subtitle = amp_label,
    x        = "Month", y = "Mean d-excess (‰)", colour = NULL
  ) +
  theme_minimal(base_size = 11)

print(p_sin)
ggsave(file.path(out_dir, "KZ6.1_sinusoidal_dexcess_seasonality.pdf"),
       p_sin, width = 7, height = 4)

# ── [KZ7.1]  SPbI preservation of anomalous d-excess events ─
q   <- quantile(base$d_ex_orig, c(0.25, 0.75), na.rm = TRUE)
iqr <- diff(q)
base <- base %>%
  mutate(anomaly = d_ex_orig < q[1] - 1.5 * iqr |
                   d_ex_orig > q[2] + 1.5 * iqr)

anom_perf <- base %>%
  group_by(method, anomaly) %>%
  summarise(
    MAD = mean(abs(err_dex), na.rm = TRUE),
    n   = n(),
    .groups = "drop"
  ) %>%
  mutate(label = ifelse(anomaly, "Anomalous", "Normal"))

p_anom <- ggplot(anom_perf,
                 aes(x = method, y = MAD, fill = label)) +
  geom_col(position = "dodge", colour = "white", linewidth = 0.3) +
  scale_fill_manual(values = c("Normal" = "#74C476", "Anomalous" = "#E41A1C")) +
  labs(
    title = "[KZ7.1] d-excess MAD: normal vs anomalous observed events",
    x     = NULL, y = "MAD d-excess (‰)", fill = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))

print(p_anom)
ggsave(file.path(out_dir, "KZ7.1_SPbI_anomalous_dexcess.pdf"),
       p_anom, width = 7, height = 4)

# ── [KZ seasonal]  Seasonal variation of d-excess errors ────
seasonal <- base %>%
  group_by(method, mon_lab) %>%
  summarise(MAD = mean(abs(err_dex), na.rm = TRUE), .groups = "drop")

p_seas <- ggplot(seasonal,
                 aes(x = mon_lab, y = MAD,
                     colour = method, group = method)) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.8) +
  scale_colour_brewer(palette = "Dark2") +
  labs(
    title  = "[KZ seasonal] Seasonal variation of d-excess imputation error",
    x      = "Month", y = "MAD d-excess (‰)", colour = NULL
  ) +
  theme_minimal(base_size = 11)

print(p_seas)
ggsave(file.path(out_dir, "KZ_seasonal_dexcess_error.pdf"),
       p_seas, width = 8, height = 4)

message("\nAll d-excess analyses complete. PDFs written to: ", out_dir)
