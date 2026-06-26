# ============================================================
# Stage 2 — CI tests, threshold detection, and significance plots
# Hatvani I.G. & Kern Z. — Mind the gap (HESS, 2026)
# https://github.com/istvan60/SPbI
#
# Input : spbi_nonoverlapping_band_errors.rds
#         paired_SPbI_vs_baselines_nonoverlapping_bands.rds
#         baseline_errors_site_X_boot_method_isotope.rds
#         (all from Stage 1)
# Output: paired_difference_CI_and_tests.csv
#         threshold_summary_SPbI_better_than_both.csv
#         plot_SPbI_bandwise_box_significance_with_baselines_*.png/pdf
#
# Run in a fresh RStudio session.
# ============================================================

# install.packages(c("tidyverse", "data.table", "scales"))
library(tidyverse)
library(data.table)
library(scales)
library(grid)


# -----------------------------
# 0) User settings
# -----------------------------

extract_dir <- "SPbI_threshold_extracted"

out_dir   <- "D:/HIG/teszt, hogy megy-e/SPbI_threshold_extracted"

file.exists(file.path(out_dir, "spbi_nonoverlapping_band_errors_d18O.rds"))
file.exists(file.path(out_dir, "spbi_nonoverlapping_band_errors_d2H.rds"))
file.exists(file.path(out_dir, "baseline_errors_site_X_boot_method_isotope.rds"))
file.exists(file.path(out_dir, "paired_SPbI_vs_baselines_nonoverlapping_bands.rds"))

setwd(out_dir)

chunk_dir <- "~/Downloads/tst/SPbI_threshold_extracted/stage1_parallel_chunks"

# Confidence interval method:
# "t" is fast.
# "bootstrap" is more robust but slower.
CI_METHOD <- "bootstrap"
B_BOOT <- 2000
CONF_LEVEL <- 0.95

# Final significance criterion:
# "CI_only"       = CI of paired difference entirely above zero
# "CI_and_Holm"   = CI above zero AND Holm-adjusted Wilcoxon p < 0.05
FINAL_CRITERION <- "CI_and_Holm"

set.seed(123)

# -----------------------------
# 1) Load compact data
# -----------------------------

spbi <- readRDS(file.path(out_dir, "spbi_nonoverlapping_band_errors.rds"))

baseline <- readRDS(file.path(out_dir, "baseline_errors_site_X_boot_method_isotope.rds"))

paired <- readRDS(file.path(out_dir, "paired_SPbI_vs_baselines_nonoverlapping_bands.rds"))

spbi <- as_tibble(spbi)
baseline <- as_tibble(baseline)
paired <- as_tibble(paired)

# Formatting
x_levels <- percent(sort(unique(paired$X)), accuracy = 1)

band_levels <- paired %>%
  distinct(band_low, band_high, band_label) %>%
  arrange(band_low, band_high) %>%
  pull(band_label)

format_common <- function(dat) {
  dat %>%
    mutate(
      X_pct = factor(percent(X, accuracy = 1), levels = x_levels),
      band_label = factor(band_label, levels = band_levels),
      metric = factor(metric, levels = c("MAD", "RMSE")),
      isotope = factor(isotope, levels = c("d18O", "d2H")),
      method = tolower(as.character(method))
    )
}

paired <- paired %>% format_common()

spbi <- spbi %>%
  mutate(
    X_pct = factor(percent(X, accuracy = 1), levels = x_levels),
    band_label = factor(band_label, levels = band_levels),
    isotope = factor(isotope, levels = c("d18O", "d2H"))
  )

baseline <- baseline %>%
  mutate(
    X_pct = factor(percent(X, accuracy = 1), levels = x_levels),
    isotope = factor(isotope, levels = c("d18O", "d2H")),
    method = tolower(as.character(method))
  )

# -----------------------------
# 2) Helper functions
# -----------------------------

mean_ci <- function(x, conf = 0.95, method = c("t", "bootstrap"), B = 2000) {
  
  method <- match.arg(method)
  x <- x[is.finite(x)]
  n <- length(x)
  
  if (n == 0) {
    return(tibble(
      n = 0L,
      mean = NA_real_,
      median = NA_real_,
      sd = NA_real_,
      ci_low = NA_real_,
      ci_high = NA_real_
    ))
  }
  
  x_mean <- mean(x)
  x_median <- median(x)
  x_sd <- sd(x)
  
  if (n < 3 || !is.finite(x_sd) || x_sd == 0) {
    return(tibble(
      n = n,
      mean = x_mean,
      median = x_median,
      sd = x_sd,
      ci_low = x_mean,
      ci_high = x_mean
    ))
  }
  
  if (method == "t") {
    
    alpha <- 1 - conf
    se <- x_sd / sqrt(n)
    crit <- qt(1 - alpha / 2, df = n - 1)
    
    ci_low <- x_mean - crit * se
    ci_high <- x_mean + crit * se
    
  } else {
    
    alpha <- 1 - conf
    
    boot_means <- replicate(
      B,
      mean(sample(x, size = n, replace = TRUE))
    )
    
    ci_low <- unname(quantile(boot_means, alpha / 2, na.rm = TRUE))
    ci_high <- unname(quantile(boot_means, 1 - alpha / 2, na.rm = TRUE))
  }
  
  tibble(
    n = n,
    mean = x_mean,
    median = x_median,
    sd = x_sd,
    ci_low = ci_low,
    ci_high = ci_high
  )
}

wilcox_greater_p <- function(x) {
  
  x <- x[is.finite(x)]
  
  if (length(x) < 3 || length(unique(x)) < 2) {
    return(NA_real_)
  }
  
  suppressWarnings(
    wilcox.test(
      x,
      mu = 0,
      alternative = "greater",
      exact = FALSE
    )$p.value
  )
}

# -----------------------------
# 3) CI for SPbI MAD/RMSE in each band
# -----------------------------

spbi_long <- spbi %>%
  pivot_longer(
    cols = c(MAD, RMSE),
    names_to = "metric",
    values_to = "error"
  ) %>%
  mutate(metric = factor(metric, levels = c("MAD", "RMSE")))

spbi_ci <- spbi_long %>%
  group_by(isotope, X, X_pct, metric, band_low, band_high, band_label) %>%
  group_modify(~ mean_ci(
    .x$error,
    conf = CONF_LEVEL,
    method = CI_METHOD,
    B = B_BOOT
  )) %>%
  ungroup()

write_csv(
  spbi_ci, "SPbI_bandwise_MAD_RMSE_CI.csv")

# -----------------------------
# 4) CI for baseline MAD/RMSE
# -----------------------------

baseline_long <- baseline %>%
  pivot_longer(
    cols = c(MAD, RMSE),
    names_to = "metric",
    values_to = "error"
  ) %>%
  mutate(metric = factor(metric, levels = c("MAD", "RMSE")))

baseline_ci <- baseline_long %>%
  group_by(isotope, X, X_pct, metric, method) %>%
  group_modify(~ mean_ci(
    .x$error,
    conf = CONF_LEVEL,
    method = CI_METHOD,
    B = B_BOOT
  )) %>%
  ungroup()

write_csv(
  baseline_ci, "baseline_MAD_RMSE_CI.csv")

# -----------------------------
# 5) Paired CI of baseline - SPbI difference
# -----------------------------
# Positive difference means:
# baseline error > SPbI error
# Therefore SPbI is better.
#
# SPbI significantly better by CI if:
# ci_low > 0

diff_ci <- paired %>%
  dplyr::group_by(
    isotope, X, X_pct, metric,
    band_low, band_high, band_label, method
  ) %>%
  dplyr::group_modify(~ mean_ci(
    .x$error_difference,
    conf = CONF_LEVEL,
    method = CI_METHOD,
    B = B_BOOT
  )) %>%
  dplyr::ungroup() %>%
  dplyr::rename(
    mean_difference   = mean,
    median_difference = median,
    sd_difference     = sd,
    diff_ci_low       = ci_low,
    diff_ci_high      = ci_high,
    n_paired          = n
  )

# -----------------------------
# 6) Paired one-sided Wilcoxon tests
# -----------------------------
# Positive error_difference means:
# baseline error > SPbI error
# Therefore SPbI is better.

group_vars <- c(
  "isotope", "X", "X_pct", "metric",
  "band_low", "band_high", "band_label", "method"
)

# Optional sanity check
missing_group_vars <- setdiff(group_vars, names(paired))

if (length(missing_group_vars) > 0) {
  stop(
    "paired is missing these required columns: ",
    paste(missing_group_vars, collapse = ", ")
  )
}

diff_tests <- paired %>%
  dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) %>%
  dplyr::summarise(
    p_value_wilcox_one_sided = wilcox_greater_p(.data[["error_difference"]]),
    .groups = "drop"
  ) %>%
  dplyr::group_by(
    .data[["isotope"]],
    .data[["X"]],
    .data[["metric"]]
  ) %>%
  dplyr::mutate(
    p_adj_holm = stats::p.adjust(
      .data[["p_value_wilcox_one_sided"]],
      method = "holm"
    )
  ) %>%
  dplyr::ungroup()


diff_results <- diff_ci %>%
  left_join(
    diff_tests,
    by = c(
      "isotope", "X", "X_pct", "metric",
      "band_low", "band_high", "band_label", "method"
    )
  ) %>%
  mutate(
    spbi_better_by_CI = diff_ci_low > 0,
    spbi_better_by_Holm_test =
      mean_difference > 0 &
      !is.na(p_adj_holm) &
      p_adj_holm < 0.05,
    spbi_better_final = case_when(
      FINAL_CRITERION == "CI_only" ~ spbi_better_by_CI,
      FINAL_CRITERION == "CI_and_Holm" ~ spbi_better_by_CI & spbi_better_by_Holm_test,
      TRUE ~ spbi_better_by_CI
    )
  )

write_csv(diff_results, "paired_difference_CI_and_tests.csv")



# -----------------------------
# 7) Threshold detection
# -----------------------------



decision_by_band <- diff_results %>%
  dplyr::group_by(
    .data[["isotope"]],
    .data[["X"]],
    .data[["X_pct"]],
    .data[["metric"]],
    .data[["band_low"]],
    .data[["band_high"]],
    .data[["band_label"]]
  ) %>%
  dplyr::summarise(
    n_methods_tested = dplyr::n_distinct(.data[["method"]]),
    
    better_than_linear = any(
      .data[["method"]] == "linear" &
        .data[["spbi_better_final"]],
      na.rm = TRUE
    ),
    
    better_than_sinusoid = any(
      .data[["method"]] == "sinusoidal" &
        .data[["spbi_better_final"]],
      na.rm = TRUE
    ),
    
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    better_than_both = .data[["better_than_linear"]] &
      .data[["better_than_sinusoid"]]
  )

threshold_summary <- decision_by_band %>%
  dplyr::group_by(
    .data[["isotope"]],
    .data[["X"]],
    .data[["X_pct"]],
    .data[["metric"]]
  ) %>%
  dplyr::summarise(
    farthest_supported_band_high_km = if (any(.data[["better_than_both"]], na.rm = TRUE)) {
      max(.data[["band_high"]][.data[["better_than_both"]]], na.rm = TRUE)
    } else {
      NA_real_
    },
    
    supported_bands = if (any(.data[["better_than_both"]], na.rm = TRUE)) {
      paste(
        as.character(.data[["band_label"]][.data[["better_than_both"]]]),
        collapse = "; "
      )
    } else {
      NA_character_
    },
    
    .groups = "drop"
  )

write_csv(
  decision_by_band, "decision_by_band_SPbI_better_than_both.csv")

write_csv(
  threshold_summary, "threshold_summary_SPbI_better_than_both.csv")

print(threshold_summary)

# -----------------------------
# 8) Final plot: SPbI boxes with baseline means and significance colouring
# -----------------------------
# Box = Q1-Q3 of SPbI errors
# Middle line = mean SPbI error
# Baseline lines = mean Linear and Sinusoidal errors
#
# Box fill shows whether SPbI is significantly better than:
# - both baselines
# - Linear only
# - Sinusoidal only
# - neither baseline
#
# Significance is taken from diff_results$spbi_better_final,
# therefore it follows FINAL_CRITERION defined above.

# ------------------------------------------------------------
# 8.1) Prepare SPbI box statistics from paired data
# ------------------------------------------------------------
# Use paired, not spbi_ci, because paired contains the exact
# SPbI-vs-baseline comparison structure and usually both isotopes.

spbi_box_source <- paired %>%
  distinct(
    Site, X, X_pct, boot,
    isotope, metric,
    band_low, band_high, band_label,
    spbi_error
  ) %>%
  transmute(
    Site,
    X,
    X_pct,
    boot,
    isotope = as.character(isotope),
    metric = as.character(metric),
    band_low,
    band_high,
    band_label = as.character(band_label),
    error = spbi_error
  ) %>%
  filter(is.finite(error))

spbi_box_stats <- spbi_box_source %>%
  dplyr::group_by(
    .data[["isotope"]],
    .data[["X"]],
    .data[["X_pct"]],
    .data[["metric"]],
    .data[["band_low"]],
    .data[["band_high"]],
    .data[["band_label"]]
  ) %>%
  dplyr::summarise(
    n = dplyr::n(),
    q_low = stats::quantile(.data[["error"]], 0.25, na.rm = TRUE),
    mean_error = mean(.data[["error"]], na.rm = TRUE),
    q_high = stats::quantile(.data[["error"]], 0.75, na.rm = TRUE),
    .groups = "drop"
  )
# ------------------------------------------------------------
# 8.2) Prepare significance classes for box fill
# ------------------------------------------------------------


sig_for_boxes <- diff_results %>%
  dplyr::mutate(
    isotope = as.character(.data[["isotope"]]),
    metric = as.character(.data[["metric"]]),
    band_label = as.character(.data[["band_label"]]),
    
    method_clean = dplyr::case_when(
      tolower(as.character(.data[["method"]])) == "linear" ~ "Linear",
      tolower(as.character(.data[["method"]])) %in% c("sinusoid", "sinusoidal") ~ "Sinusoidal",
      TRUE ~ as.character(.data[["method"]])
    ),
    
    sig = .data[["spbi_better_final"]] %in% TRUE
  ) %>%
  dplyr::group_by(
    .data[["isotope"]],
    .data[["X"]],
    .data[["X_pct"]],
    .data[["metric"]],
    .data[["band_low"]],
    .data[["band_high"]],
    .data[["band_label"]]
  ) %>%
  dplyr::summarise(
    sig_vs_linear = any(
      .data[["method_clean"]] == "Linear" & .data[["sig"]],
      na.rm = TRUE
    ),
    
    sig_vs_sinusoidal = any(
      .data[["method_clean"]] == "Sinusoidal" & .data[["sig"]],
      na.rm = TRUE
    ),
    
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    sig_class = dplyr::case_when(
      .data[["sig_vs_linear"]] & .data[["sig_vs_sinusoidal"]] ~
        "SPbI significantly better than both",
      
      .data[["sig_vs_linear"]] & !.data[["sig_vs_sinusoidal"]] ~
        "SPbI better than Linear only",
      
      !.data[["sig_vs_linear"]] & .data[["sig_vs_sinusoidal"]] ~
        "SPbI better than Sinusoidal only",
      
      TRUE ~
        "Not significantly better than either"
    )
  )


spbi_box_stats <- spbi_box_stats %>%
  left_join(
    sig_for_boxes,
    by = c(
      "isotope", "X", "X_pct", "metric",
      "band_low", "band_high", "band_label"
    )
  ) %>%
  mutate(
    sig_class = replace_na(sig_class, "No test available"),
    isotope = factor(isotope, levels = c("d18O", "d2H")),
    metric = factor(metric, levels = c("MAD", "RMSE")),
    band_label = factor(band_label, levels = band_levels),
    sig_class = factor(
      sig_class,
      levels = c(
        "SPbI significantly better than both",
        "SPbI better than Linear only",
        "SPbI better than Sinusoidal only",
        "Not significantly better than either",
        "No test available"
      )
    )
  )

# ------------------------------------------------------------
# 8.3) Prepare baseline mean lines
# ------------------------------------------------------------

baseline_lines <- baseline_long %>%
  dplyr::mutate(
    isotope = as.character(.data[["isotope"]]),
    metric  = as.character(.data[["metric"]]),
    method  = dplyr::case_when(
      tolower(as.character(.data[["method"]])) == "linear" ~ "Linear",
      tolower(as.character(.data[["method"]])) %in% c("sinusoid", "sinusoidal") ~ "Sinusoidal",
      TRUE ~ as.character(.data[["method"]])
    )
  ) %>%
  dplyr::filter(
    .data[["method"]] %in% c("Linear", "Sinusoidal"),
    is.finite(.data[["error"]])
  ) %>%
  dplyr::group_by(
    .data[["isotope"]],
    .data[["X"]],
    .data[["X_pct"]],
    .data[["metric"]],
    .data[["method"]]
  ) %>%
  dplyr::summarise(
    mean_error = mean(.data[["error"]], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    isotope = factor(.data[["isotope"]], levels = c("d18O", "d2H")),
    metric  = factor(.data[["metric"]], levels = c("MAD", "RMSE")),
    method  = factor(.data[["method"]], levels = c("Linear", "Sinusoidal"))
  )

# ------------------------------------------------------------
# 8.4) Plot function
# Uniform y-axis min/max per X_pct row
# Baseline hlines retained
# ------------------------------------------------------------

plot_error_with_baselines <- function(iso_name) {
  
  spbi_plot <- spbi_box_stats %>%
    dplyr::filter(as.character(.data[["isotope"]]) == iso_name)
  
  base_plot <- baseline_lines %>%
    dplyr::filter(as.character(.data[["isotope"]]) == iso_name)
  
  message("Plotting ", iso_name)
  message("SPbI box rows: ", nrow(spbi_plot))
  message("Baseline rows: ", nrow(base_plot))
  
  # Diagnostic: check that baselines are really present
  print(
    base_plot %>%
      dplyr::count(.data[["X_pct"]], .data[["metric"]], .data[["method"]])
  )
  
  # ----------------------------------------------------------
  # Build invisible row-wise y-limits
  # These limits are calculated from:
  #   1) SPbI IQR boxes: q_low to q_high
  #   2) baseline mean lines: mean_error
  # grouped only by X_pct, so both metric columns in a row
  # get the same y-axis range.
  # ----------------------------------------------------------
  
  y_from_spbi <- spbi_plot %>%
    dplyr::transmute(
      X_pct = .data[["X_pct"]],
      y_min = .data[["q_low"]],
      y_max = .data[["q_high"]]
    )
  
  y_from_base <- base_plot %>%
    dplyr::transmute(
      X_pct = .data[["X_pct"]],
      y_min = .data[["mean_error"]],
      y_max = .data[["mean_error"]]
    )
  
  y_limits_row <- dplyr::bind_rows(y_from_spbi, y_from_base) %>%
    dplyr::filter(
      is.finite(.data[["y_min"]]),
      is.finite(.data[["y_max"]])
    ) %>%
    dplyr::group_by(.data[["X_pct"]]) %>%
    dplyr::summarise(
      y_min = min(.data[["y_min"]], na.rm = TRUE),
      y_max = max(.data[["y_max"]], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      y_pad = 0.06 * (.data[["y_max"]] - .data[["y_min"]]),
      y_pad = dplyr::if_else(
        is.finite(.data[["y_pad"]]) & .data[["y_pad"]] > 0,
        .data[["y_pad"]],
        0.1
      ),
      y_min = .data[["y_min"]] - .data[["y_pad"]],
      y_max = .data[["y_max"]] + .data[["y_pad"]]
    )
  
  # Use the existing facet levels from the actual data
  metric_levels <- levels(spbi_box_stats$metric)
  if (is.null(metric_levels)) {
    metric_levels <- unique(as.character(spbi_box_stats$metric))
  }
  
  first_band <- levels(spbi_box_stats$band_label)[1]
  if (is.null(first_band) || is.na(first_band)) {
    first_band <- unique(as.character(spbi_plot$band_label))[1]
  }
  
  y_blank <- y_limits_row %>%
    tidyr::crossing(
      metric = factor(metric_levels, levels = metric_levels),
      band_label = factor(
        first_band,
        levels = levels(spbi_box_stats$band_label)
      )
    ) %>%
    tidyr::pivot_longer(
      cols = c("y_min", "y_max"),
      names_to = "limit_type",
      values_to = "y"
    )
  
  ggplot() +
    
    # Invisible layer forcing row-wise y-limits
    geom_blank(
      data = y_blank,
      aes(
        x = band_label,
        y = y
      )
    ) +
    
    # SPbI IQR box with mean as central line
    geom_crossbar(
      data = spbi_plot,
      aes(
        x = band_label,
        y = mean_error,
        ymin = q_low,
        ymax = q_high,
        fill = sig_class
      ),
      width = 0.55,
      linewidth = 0.45,
      colour = "black",
      alpha = 0.75
    ) +
    
    # Baseline mean lines — unchanged from original
    geom_hline(
      data = base_plot,
      aes(
        yintercept = mean_error,
        colour = method,
        linetype = method
      ),
      linewidth = 1.1,
      alpha = 0.95
    ) +
    
    facet_grid(
      rows = vars(X_pct),
      cols = vars(metric),
      scales = "free_y",
      switch = "y"
    ) +
    
    scale_fill_manual(
      name = "SPbI significance",
      values = c(
        "SPbI significantly better than both" = "#2ca25f",
        "SPbI better than Linear only" = "#a1d99b",
        "SPbI better than Sinusoidal only" = "#fdae6b",
        "Not significantly better than either" = "grey80",
        "No test available" = "white"
      ),
      drop = FALSE
    ) +
    
    scale_linetype_manual(
      name = "Baseline",
      values = c(
        "Linear" = "dashed",
        "Sinusoidal" = "dotted"
      ),
      drop = FALSE
    ) +
    
    scale_colour_manual(
      name = "Baseline",
      values = c(
        "Linear" = "brown",
        "Sinusoidal" = "#40E0D0"
      ),
      drop = FALSE
    ) +
    
    labs(
      title = NULL,
      x = "Non-overlapping distance band",
      y = "SPbI error (‰)"
    ) +
    
    theme_minimal(base_size = 16) +
    theme(
      legend.position = "top",
      legend.key.width = unit(1.6, "cm"),
      strip.placement = "outside",
      strip.text = element_text(face = "bold", size = 16),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 13),
      axis.text.y = element_text(size = 13),
      axis.title = element_text(face = "bold", size = 18),
      panel.spacing.x = unit(0.35, "lines"),
      panel.spacing.y = unit(0.35, "lines"),
      plot.margin = margin(6, 8, 6, 8)
    ) +
    
    guides(
      fill = guide_legend(nrow = 2),
      colour = guide_legend(override.aes = list(linewidth = 1.1)),
      linetype = guide_legend(override.aes = list(linewidth = 1.1))
    )
}
# ------------------------------------------------------------
# 8.5) Make and save plots
# ------------------------------------------------------------

p_error_d18O <- plot_error_with_baselines("d18O")
p_error_d2H  <- plot_error_with_baselines("d2H")

print(p_error_d18O)
print(p_error_d2H)

ggsave(
  "plot_SPbI_bandwise_box_significance_with_baselines_d18O.png",
  p_error_d18O,
  width = 9,
  height = 14,
  dpi = 300
)

ggsave(
  "plot_SPbI_bandwise_box_significance_with_baselines_d2H.png",
  p_error_d2H,
  width = 9,
  height = 14,
  dpi = 300
)

ggsave(
  "plot_SPbI_bandwise_box_significance_with_baselines_d18O.pdf",
  p_error_d18O,
  width = 9,
  height = 14
)

ggsave(
  "plot_SPbI_bandwise_box_significance_with_baselines_d2H.pdf",
  p_error_d2H,
  width = 9,
  height = 14
)



# -----------------------------
# 9) Optional manuscript-ready compact table
# -----------------------------

# ------------------------------------------------------------
# Manuscript table from current diff_results
# Works whether Wilcoxon/Holm columns exist or not
# ------------------------------------------------------------


manuscript_table <- diff_results %>%
  dplyr::select(
    .data[["isotope"]],
    .data[["X_pct"]],
    .data[["metric"]],
    .data[["band_label"]],
    .data[["method"]],
    .data[["n_paired"]],
    .data[["mean_difference"]],
    .data[["diff_ci_low"]],
    .data[["diff_ci_high"]],
    .data[["p_value_wilcox_one_sided"]],
    .data[["p_adj_holm"]],
    .data[["spbi_better_by_CI"]],
    .data[["spbi_better_by_Holm_test"]],
    .data[["spbi_better_final"]]
  ) %>%
  dplyr::arrange(
    .data[["isotope"]],
    .data[["X_pct"]],
    .data[["metric"]],
    .data[["band_label"]],
    .data[["method"]]
  )

data.table::fwrite(
    manuscript_table,
    file.path(out_dir, "manuscript_table_paired_difference_tests.csv")
)

message("Done.")
message("Important outputs:")
message("- paired_difference_CI_and_tests.csv")
message("- threshold_summary_SPbI_better_than_both.csv")
message("- plot_paired_difference_baseline_minus_SPbI_d18O.png")
message("- plot_paired_difference_baseline_minus_SPbI_d2H.png")
message("- plot_decision_heatmap_SPbI_better_than_both.png")




















# ============================================================
# Reviewer response analysis
# Question 13: fallback rates + gap structure
#
# Version with explicit package references.
# This is intended for sessions where kgc keeps plyr attached.
# Do NOT detach plyr. Instead, all dplyr/tidyr/lubridate/ggplot2
# calls are explicitly namespaced.
#
# Run in the R session where Full_futas.RData is already loaded,
# OR uncomment load() below.
#
# Required objects in the session:
#   all_imputed
#   df_O18_trim
#   removed_dates
# ============================================================

# ============================================================
# PART 1: FALLBACK RATES
# ============================================================

message("Part 1: fallback rates ...")

fallback_raw <- all_imputed %>%
  dplyr::group_by(.data[["method"]], .data[["X"]], .data[["X_pct"]]) %>%
  dplyr::summarise(
    n_total    = dplyr::n(),
    n_fail_O18 = sum(!is.finite(.data[["O18_imp"]])),
    n_fail_H2  = sum(!is.finite(.data[["H2_imp"]])),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    pct_fail_O18 = 100 * .data[["n_fail_O18"]] / .data[["n_total"]],
    pct_fail_H2  = 100 * .data[["n_fail_H2"]]  / .data[["n_total"]]
  ) %>%
  dplyr::arrange(.data[["method"]], .data[["X"]])

utils::write.csv(
  fallback_raw,
  file.path(out_dir, "fallback_rates_table.csv"),
  row.names = FALSE
)

fallback_wide <- dplyr::bind_rows(
  fallback_raw %>%
    dplyr::select(.data[["method"]], .data[["X_pct"]], .data[["pct_fail_O18"]]) %>%
    tidyr::pivot_wider(
      names_from = .data[["X_pct"]],
      values_from = .data[["pct_fail_O18"]]
    ) %>%
    dplyr::mutate(isotope = "d18O") %>%
    dplyr::select(.data[["isotope"]], dplyr::everything()),
  
  fallback_raw %>%
    dplyr::select(.data[["method"]], .data[["X_pct"]], .data[["pct_fail_H2"]]) %>%
    tidyr::pivot_wider(
      names_from = .data[["X_pct"]],
      values_from = .data[["pct_fail_H2"]]
    ) %>%
    dplyr::mutate(isotope = "d2H") %>%
    dplyr::select(.data[["isotope"]], dplyr::everything())
)

utils::write.csv(
  fallback_wide,
  file.path(out_dir, "fallback_rates_wide.csv"),
  row.names = FALSE
)

p_fallback <- fallback_raw %>%
  ggplot2::ggplot(
    ggplot2::aes(
      x = .data[["X_pct"]],
      y = .data[["method"]],
      fill = .data[["pct_fail_O18"]]
    )
  ) +
  ggplot2::geom_tile(color = "white") +
  ggplot2::geom_text(
    ggplot2::aes(label = sprintf("%.1f%%", .data[["pct_fail_O18"]])),
    size = 3.5
  ) +
  ggplot2::scale_fill_gradient(
    low = "white",
    high = "#d73027",
    name = "Fallback %"
  ) +
  ggplot2::labs(
    title    = "Fallback rates per method and masking fraction (δ18O)",
    subtitle = "Fallback = imputation returned NA",
    x = "Masking fraction",
    y = NULL
  ) +
  ggplot2::theme_minimal(base_size = 13) +
  ggplot2::theme(legend.position = "right")

print(p_fallback)

ggplot2::ggsave(
  file.path(out_dir, "fallback_heatmap_d18O.png"),
  p_fallback,
  width = 9,
  height = 5,
  dpi = 300
)

# ============================================================
# PART 2: SEASONAL PATTERN OF REAL MISSING DATA
# ============================================================

message("Part 2: seasonal missingness ...")

gap_by_site <- df_O18_trim %>%
  dplyr::mutate(
    month_floor = lubridate::floor_date(.data[["Date"]], "month")
  ) %>%
  dplyr::filter(!is.na(.data[["month_floor"]])) %>%
  dplyr::group_by(.data[["Site"]]) %>%
  dplyr::summarise(
    obs_months = list(sort(unique(.data[["month_floor"]]))),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    full_seq = purrr::map(
      .data[["obs_months"]],
      ~ seq.Date(
        from = min(.x, na.rm = TRUE),
        to   = max(.x, na.rm = TRUE),
        by   = "month"
      )
    ),
    missing_months = purrr::map2(
      .data[["full_seq"]],
      .data[["obs_months"]],
      ~ as.Date(setdiff(.x, .y), origin = "1970-01-01")
    )
  )

seasonal_missing <- gap_by_site %>%
  dplyr::select(.data[["Site"]], .data[["missing_months"]]) %>%
  tidyr::unnest(
    cols = c("missing_months"),
    keep_empty = FALSE
  ) %>%
  dplyr::mutate(
    calendar_month = factor(
      lubridate::month(.data[["missing_months"]]),
      levels = 1:12,
      labels = month.abb
    )
  ) %>%
  dplyr::count(
    .data[["calendar_month"]],
    name = "n_missing_months",
    .drop = FALSE
  ) %>%
  tidyr::complete(
    calendar_month = factor(month.abb, levels = month.abb),
    fill = list(n_missing_months = 0L)
  )

expected_per_month <- sum(seasonal_missing$n_missing_months, na.rm = TRUE) / 12

seasonal_missing <- seasonal_missing %>%
  dplyr::mutate(
    expected = expected_per_month,
    ratio_obs_exp = dplyr::if_else(
      .data[["expected"]] > 0,
      .data[["n_missing_months"]] / .data[["expected"]],
      NA_real_
    )
  )

utils::write.csv(
  seasonal_missing,
  file.path(out_dir, "seasonal_missingness.csv"),
  row.names = FALSE
)

p_seasonal <- seasonal_missing %>%
  ggplot2::ggplot(
    ggplot2::aes(
      x = .data[["calendar_month"]],
      y = .data[["n_missing_months"]]
    )
  ) +
  ggplot2::geom_col(fill = "#4393c3") +
  ggplot2::geom_hline(
    yintercept = expected_per_month,
    linetype = "dashed",
    color = "red",
    linewidth = 1
  ) +
  ggplot2::annotate(
    "text",
    x = "Jan",
    y = expected_per_month * 1.05,
    label = "Expected if MCAR",
    hjust = 0,
    color = "red",
    size = 3.5
  ) +
  ggplot2::labs(
    title    = "Seasonal distribution of real missing months",
    subtitle = "Dashed line = expected count under MCAR (uniform across months)",
    x = "Calendar month",
    y = "Number of missing months (all sites combined)"
  ) +
  ggplot2::theme_minimal(base_size = 13)

print(p_seasonal)

ggplot2::ggsave(
  file.path(out_dir, "seasonal_missingness.png"),
  p_seasonal,
  width = 8,
  height = 5,
  dpi = 300
)

# ============================================================
# PART 3: GAP LENGTH DISTRIBUTION IN THE BOOTSTRAP MASKING
# ============================================================

message("Part 3: bootstrap gap structure by rarefaction ...")

calc_run_lengths <- function(dates) {
  dates <- sort(unique(as.Date(dates)))
  
  if (length(dates) == 0) return(integer(0))
  if (length(dates) == 1) return(1L)
  
  diffs <- as.numeric(diff(dates))
  break_after <- which(diffs > 32)
  breaks <- c(0L, break_after, length(dates))
  
  as.integer(diff(breaks))
}

boot_gap_runs <- removed_dates %>%
  dplyr::mutate(Date = as.Date(.data[["Date"]])) %>%
  dplyr::arrange(.data[["Site"]], .data[["X"]], .data[["boot"]], .data[["Date"]]) %>%
  dplyr::group_by(.data[["Site"]], .data[["X"]], .data[["boot"]]) %>%
  dplyr::summarise(
    run_lengths = list(calc_run_lengths(.data[["Date"]])),
    .groups = "drop"
  ) %>%
  tidyr::unnest_longer(
    col = .data[["run_lengths"]],
    values_to = "run_length"
  ) %>%
  dplyr::mutate(run_length = as.integer(.data[["run_length"]])) %>%
  dplyr::filter(is.finite(.data[["run_length"]]))

boot_gap_summary <- boot_gap_runs %>%
  dplyr::mutate(
    X_pct = factor(
      paste0(round(as.numeric(.data[["X"]]) * 100), "%"),
      levels = x_levels
    )
  ) %>%
  dplyr::group_by(.data[["X_pct"]]) %>%
  dplyr::summarise(
    n_single_month   = sum(.data[["run_length"]] == 1, na.rm = TRUE),
    n_2month         = sum(.data[["run_length"]] == 2, na.rm = TRUE),
    n_3month         = sum(.data[["run_length"]] == 3, na.rm = TRUE),
    n_4to6month      = sum(.data[["run_length"]] >= 4 & .data[["run_length"]] <= 6, na.rm = TRUE),
    n_gt6month       = sum(.data[["run_length"]] > 6, na.rm = TRUE),
    n_runs           = dplyr::n(),
    pct_single_month = 100 * .data[["n_single_month"]] / .data[["n_runs"]],
    pct_2plus_month  = 100 * (.data[["n_runs"]] - .data[["n_single_month"]]) / .data[["n_runs"]],
    median_run       = stats::median(.data[["run_length"]], na.rm = TRUE),
    max_run          = max(.data[["run_length"]], na.rm = TRUE),
    .groups = "drop"
  )

utils::write.csv(
  boot_gap_summary,
  file.path(out_dir, "bootstrap_consecutive_gaps.csv"),
  row.names = FALSE
)

message("  Bootstrap gap structure:")
print(boot_gap_summary)

boot_run_dist <- boot_gap_runs %>%
  dplyr::mutate(
    X_pct = factor(
      paste0(round(as.numeric(.data[["X"]]) * 100), "%"),
      levels = x_levels
    ),
    run_capped = pmin(.data[["run_length"]], 12L)
  ) %>%
  dplyr::count(.data[["X_pct"]], .data[["run_capped"]]) %>%
  dplyr::group_by(.data[["X_pct"]]) %>%
  dplyr::mutate(
    pct = 100 * .data[["n"]] / sum(.data[["n"]], na.rm = TRUE)
  ) %>%
  dplyr::ungroup() %>%
  tidyr::complete(
    X_pct = factor(x_levels, levels = x_levels),
    run_capped = 1:12,
    fill = list(n = 0L, pct = 0)
  )

p_boot_runs <- boot_run_dist %>%
  ggplot2::ggplot(
    ggplot2::aes(
      x = .data[["run_capped"]],
      y = .data[["pct"]]
    )
  ) +
  ggplot2::geom_col(fill = "#2166ac", width = 0.7) +
  ggplot2::scale_x_continuous(breaks = 1:12) +
  ggplot2::facet_wrap(~ X_pct, ncol = 3) +
  ggplot2::scale_y_continuous(limits = c(0, 100)) +
  ggplot2::labs(
    title    = "Gap length distribution in bootstrap masking by rarefaction level",
    subtitle = "Run lengths > 12 months capped at 12",
    x        = "Gap length (months)",
    y        = "% of masked runs"
  ) +
  ggplot2::theme_minimal(base_size = 13) +
  ggplot2::theme(
    strip.text    = ggplot2::element_text(face = "bold"),
    panel.spacing = grid::unit(0.6, "lines")
  )

print(p_boot_runs)

ggplot2::ggsave(
  file.path(out_dir, "bootstrap_run_length_distribution.png"),
  p_boot_runs,
  width = 10,
  height = 7,
  dpi = 300
)

# ============================================================
# PART 4: DO MASKED VALUES OVER-REPRESENT ISOTOPIC EXTREMES?
# ============================================================

message("Part 4: masked vs. full value distribution ...")

all_O18 <- df_O18_trim %>%
  dplyr::filter(is.finite(.data[["O18"]]))

masked_O18 <- all_imputed %>%
  dplyr::filter(
    .data[["method"]] == "Linear",
    is.finite(.data[["O18_orig"]])
  ) %>%
  dplyr::select(
    .data[["Site"]],
    .data[["Date"]],
    .data[["X"]],
    .data[["X_pct"]],
    .data[["O18_orig"]]
  )

extremes_comparison <- masked_O18 %>%
  dplyr::group_by(.data[["X_pct"]]) %>%
  dplyr::summarise(
    n           = dplyr::n(),
    mean_masked = mean(.data[["O18_orig"]], na.rm = TRUE),
    sd_masked   = stats::sd(.data[["O18_orig"]], na.rm = TRUE),
    q05_masked  = stats::quantile(.data[["O18_orig"]], 0.05, na.rm = TRUE),
    q50_masked  = stats::quantile(.data[["O18_orig"]], 0.50, na.rm = TRUE),
    q95_masked  = stats::quantile(.data[["O18_orig"]], 0.95, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    mean_all = mean(all_O18$O18, na.rm = TRUE),
    sd_all   = stats::sd(all_O18$O18, na.rm = TRUE),
    q05_all  = stats::quantile(all_O18$O18, 0.05, na.rm = TRUE),
    q50_all  = stats::quantile(all_O18$O18, 0.50, na.rm = TRUE),
    q95_all  = stats::quantile(all_O18$O18, 0.95, na.rm = TRUE)
  )

utils::write.csv(
  extremes_comparison,
  file.path(out_dir, "extremes_masked_vs_all.csv"),
  row.names = FALSE
)

p_extremes <- masked_O18 %>%
  ggplot2::ggplot(ggplot2::aes(x = .data[["O18_orig"]])) +
  ggplot2::geom_density(
    ggplot2::aes(
      color = .data[["X_pct"]],
      group = .data[["X_pct"]]
    ),
    linewidth = 0.7
  ) +
  ggplot2::geom_density(
    data = all_O18,
    ggplot2::aes(x = .data[["O18"]]),
    color = "black",
    linewidth = 1.2,
    linetype = "dashed",
    inherit.aes = FALSE
  ) +
  ggplot2::annotate(
    "text",
    x = min(all_O18$O18, na.rm = TRUE),
    y = 0,
    label = "All observed",
    hjust = 0,
    size = 3.5,
    color = "black"
  ) +
  ggplot2::scale_color_brewer(
    palette = "Reds",
    name = "Masking %"
  ) +
  ggplot2::labs(
    title    = "δ18O value distribution: masked vs. all observed",
    subtitle = "Dashed = full dataset; colours = randomly masked subset by masking fraction",
    x = "δ18O (‰)",
    y = "Density"
  ) +
  ggplot2::theme_minimal(base_size = 13) +
  ggplot2::theme(legend.position = "top")

print(p_extremes)

ggplot2::ggsave(
  file.path(out_dir, "masked_vs_all_distribution_d18O.png"),
  p_extremes,
  width = 8,
  height = 5,
  dpi = 300
)

# ============================================================
# SUMMARY
# ============================================================

message("\n========== SUMMARY ==========\n")

message("Fallback rates (d18O):")
print(
  fallback_wide %>%
    dplyr::filter(.data[["isotope"]] == "d18O")
)

message("\nBootstrap gap structure:")
print(boot_gap_summary)

message("\nMCAR check — mean d18O:")
print(
  extremes_comparison %>%
    dplyr::select(
      .data[["X_pct"]],
      .data[["mean_masked"]],
      .data[["mean_all"]],
      .data[["sd_masked"]],
      .data[["sd_all"]]
    )
)

message("\nAll files written to: ", out_dir)
















