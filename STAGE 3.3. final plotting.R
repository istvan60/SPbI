# ============================================================
# Stage 3.3 — Combined Bland-Altman, observed-vs-predicted, and performance plots
# Hatvani I.G. & Kern Z. — Mind the gap (HESS, 2026)
# https://github.com/istvan60/SPbI
#
# Input : combined_point_eval_with_ICM.rds  (from Stage 3.2)
# Output: BA_with_ICM_*.pdf/png
#         observed_vs_predicted_with_ICM_*.pdf/png
#         MAD_RMSE_with_ICM_*.pdf/png
#         combined_stats_*.csv
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(grid)
})

if (!requireNamespace("hexbin", quietly = TRUE)) {
  stop("Package 'hexbin' is required for geom_hex(). Install it with install.packages('hexbin').", call. = FALSE)
}

# -----------------------------
# USER SETTINGS
# -----------------------------
# Folder containing combined_point_eval_with_ICM.rds
base_dir <- "/Users/hatvaniistvan/Downloads/tst/SPbI_threshold_extracted/extracted_for_fig2"

# If you run the script from inside the output folder, this fallback keeps it portable.
if (!dir.exists(base_dir) && file.exists("combined_point_eval_with_ICM.rds")) {
  base_dir <- "."
}

input_file <- file.path(base_dir, "combined_point_eval_with_ICM.rds")
out_dir <- file.path(base_dir, "combined_BA_performance_LinCCC_with_ICM")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Keep all X levels by default. To shorten figures, use e.g.:
# keep_X_pct <- c("2%", "4%", "8%", "16%", "32%")
keep_X_pct <- NULL

preferred_method_order <- c(
  "LOCF", "Linear", "Spline", "Stine", "Kalman", "Moving-average",
  "Sinusoidal", "SPbI", "ICM-grid"
)

bins_ba  <- 60
bins_hex <- 70

# Bland-Altman axis limits. Set individual entries to NULL if you want free limits.
ba_axis_limits <- list(
  "d18O" = list(x = c(-40, 0),    y = c(-30, 30)),
  "d2H"  = list(x = c(-250, 50),  y = c(-200, 400)),
  "d-excess" = list(x = c(-40, 40), y = c(-40, 40))
)

# Observed-vs-predicted axis limits for visualisation.
obs_pred_axis_limits <- list(
  "d18O" = c(-40, 0),
  "d2H"  = c(-250, 50),
  "d-excess" = c(-40, 40)
)

# -----------------------------
# Helpers
# -----------------------------
stop_if_missing <- function(path) {
  if (!file.exists(path)) stop("Cannot find file: ", path, call. = FALSE)
}

safe_filename <- function(x) {
  x <- gsub("δ", "d", x, fixed = TRUE)
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("_+", "_", x)
  gsub("^_|_$", "", x)
}

label_var <- function(which_var) {
  switch(
    which_var,
    "d18O"     = "δ¹⁸O",
    "d2H"      = "δ²H",
    "d-excess" = "d-excess",
    which_var
  )
}

x_axis_label_ba <- function(which_var) {
  switch(
    which_var,
    "d18O"     = expression(Observed~delta^18*O~("‰")),
    "d2H"      = expression(Observed~delta^2*H~("‰")),
    "d-excess" = expression(Observed~d-excess~("‰")),
    "Observed value"
  )
}

y_axis_label_ba <- function(which_var) {
  switch(
    which_var,
    "d18O"     = expression(Predicted - observed~delta^18*O~("‰")),
    "d2H"      = expression(Predicted - observed~delta^2*H~("‰")),
    "d-excess" = expression(Predicted - observed~d-excess~("‰")),
    "Predicted - observed"
  )
}

x_axis_label_obs_pred <- function(which_var) {
  switch(
    which_var,
    "d18O"     = expression(Observed~delta^18*O~("‰")),
    "d2H"      = expression(Observed~delta^2*H~("‰")),
    "d-excess" = expression(Observed~d-excess~("‰")),
    "Observed value"
  )
}

y_axis_label_obs_pred <- function(which_var) {
  switch(
    which_var,
    "d18O"     = expression(Predicted~delta^18*O~("‰")),
    "d2H"      = expression(Predicted~delta^2*H~("‰")),
    "d-excess" = expression(Predicted~d-excess~("‰")),
    "Predicted value"
  )
}

fmt2 <- function(x) {
  ifelse(is.finite(x), sprintf("%.2f", x), "NA")
}

calc_stats_one <- function(orig, pred, diff) {
  ok <- is.finite(orig) & is.finite(pred) & is.finite(diff)
  x <- as.numeric(orig[ok])
  y <- as.numeric(pred[ok])
  e <- as.numeric(diff[ok])
  n <- length(e)
  
  if (n == 0) {
    return(data.table(
      n_points = 0L,
      MAD = NA_real_, RMSE = NA_real_, bias = NA_real_,
      mean_orig = NA_real_, mean_pred = NA_real_,
      sd_orig = NA_real_, sd_pred = NA_real_,
      Pearson_r = NA_real_, slope = NA_real_, intercept = NA_real_, Lin_CCC = NA_real_
    ))
  }
  
  mad_i <- mean(abs(e), na.rm = TRUE)
  rmse_i <- sqrt(mean(e^2, na.rm = TRUE))
  bias_i <- mean(e, na.rm = TRUE)
  mx <- mean(x, na.rm = TRUE)
  my <- mean(y, na.rm = TRUE)
  sx <- stats::sd(x, na.rm = TRUE)
  sy <- stats::sd(y, na.rm = TRUE)
  
  if (n < 3 || !is.finite(sx) || !is.finite(sy) || sx == 0 || sy == 0) {
    pearson_i <- NA_real_
    slope_i <- NA_real_
    intercept_i <- NA_real_
    ccc_i <- NA_real_
  } else {
    vx <- stats::var(x, na.rm = TRUE)
    vy <- stats::var(y, na.rm = TRUE)
    cov_xy <- stats::cov(x, y, use = "complete.obs")
    pearson_i <- stats::cor(x, y, use = "complete.obs")
    slope_i <- cov_xy / vx
    intercept_i <- my - slope_i * mx
    ccc_i <- (2 * cov_xy) / (vx + vy + (mx - my)^2)
  }
  
  data.table(
    n_points = as.integer(n),
    MAD = mad_i,
    RMSE = rmse_i,
    bias = bias_i,
    mean_orig = mx,
    mean_pred = my,
    sd_orig = sx,
    sd_pred = sy,
    Pearson_r = pearson_i,
    slope = slope_i,
    intercept = intercept_i,
    Lin_CCC = ccc_i
  )
}

# -----------------------------
# Load combined current-run interpolation + ICM data
# -----------------------------
stop_if_missing(input_file)
combined_eval <- readRDS(input_file)
data.table::setDT(combined_eval)

required_cols <- c("Site", "date", "X", "X_pct", "boot", "method", "var", "orig", "pred")
missing_cols <- setdiff(required_cols, names(combined_eval))
if (length(missing_cols) > 0) {
  stop(
    "combined_point_eval_with_ICM.rds is missing required column(s): ",
    paste(missing_cols, collapse = ", "),
    "\nAvailable columns: ", paste(names(combined_eval), collapse = ", "),
    call. = FALSE
  )
}

if (!"diff" %in% names(combined_eval)) {
  combined_eval[, diff := as.numeric(pred) - as.numeric(orig)]
}
if (!"source" %in% names(combined_eval)) {
  combined_eval[, source := fifelse(as.character(method) == "ICM-grid", "ICM_database_grid", "masked_interpolation")]
}

combined_eval[, `:=`(
  Site = as.character(Site),
  date = as.IDate(date),
  X = as.numeric(X),
  X_pct = as.character(X_pct),
  method = as.character(method),
  var = as.character(var),
  orig = as.numeric(orig),
  pred = as.numeric(pred),
  diff = as.numeric(diff),
  source = as.character(source)
)]

combined_eval <- combined_eval[is.finite(orig) & is.finite(pred) & is.finite(diff)]

if (!is.null(keep_X_pct)) {
  combined_eval <- combined_eval[X_pct %in% keep_X_pct]
}

if (nrow(combined_eval) == 0) {
  stop("No finite rows remain in combined_eval after filtering.", call. = FALSE)
}

# Factor order for methods and rarefaction levels.
all_methods <- unique(as.character(combined_eval$method))
method_levels <- c(
  preferred_method_order[preferred_method_order %in% all_methods],
  setdiff(sort(all_methods), preferred_method_order)
)

x_order <- unique(combined_eval[order(X)]$X_pct)
if (!is.null(keep_X_pct)) x_order <- x_order[x_order %in% keep_X_pct]

combined_eval[, method := factor(method, levels = method_levels)]
combined_eval[, X_pct := factor(X_pct, levels = x_order)]
combined_eval <- combined_eval[!is.na(X_pct)]

message("Loaded combined current-run interpolation + ICM data: ", format(nrow(combined_eval), big.mark = ","), " rows")
message("Variables: ", paste(sort(unique(as.character(combined_eval$var))), collapse = ", "))
message("Methods: ", paste(levels(combined_eval$method), collapse = ", "))
message("X levels: ", paste(levels(combined_eval$X_pct), collapse = ", "))

# -----------------------------
# Statistics tables
# -----------------------------
stats_by_X_method <- combined_eval[
  ,
  calc_stats_one(orig, pred, diff),
  by = .(source, var, X, X_pct, method)
][order(var, X, method)]

stats_overall_by_method <- combined_eval[
  ,
  calc_stats_one(orig, pred, diff),
  by = .(source, var, method)
][order(var, method)]

stats_by_boot <- combined_eval[
  ,
  calc_stats_one(orig, pred, diff),
  by = .(source, var, X, X_pct, boot, method)
][order(var, X, method, boot)]

stats_by_boot_summary <- stats_by_boot[
  ,
  .(
    n_boot = .N,
    n_points_total = sum(n_points, na.rm = TRUE),
    MAD_median = median(MAD, na.rm = TRUE),
    MAD_mean = mean(MAD, na.rm = TRUE),
    RMSE_median = median(RMSE, na.rm = TRUE),
    RMSE_mean = mean(RMSE, na.rm = TRUE),
    bias_median = median(bias, na.rm = TRUE),
    bias_mean = mean(bias, na.rm = TRUE),
    slope_median = median(slope, na.rm = TRUE),
    slope_mean = mean(slope, na.rm = TRUE),
    Pearson_r_median = median(Pearson_r, na.rm = TRUE),
    Pearson_r_mean = mean(Pearson_r, na.rm = TRUE),
    Lin_CCC_median = median(Lin_CCC, na.rm = TRUE),
    Lin_CCC_mean = mean(Lin_CCC, na.rm = TRUE)
  ),
  by = .(source, var, X, X_pct, method)
][order(var, X, method)]

ba_stats <- combined_eval[
  ,
  .(
    bias = mean(diff, na.rm = TRUE),
    n_points = .N
  ),
  by = .(source, var, X, X_pct, method)
][order(var, X, method)]

saveRDS(stats_by_X_method, file.path(out_dir, "combined_stats_by_X_method_with_ICM.rds"), compress = FALSE)
saveRDS(stats_overall_by_method, file.path(out_dir, "combined_stats_overall_by_method_with_ICM.rds"), compress = FALSE)
saveRDS(stats_by_boot, file.path(out_dir, "combined_stats_by_boot_with_ICM.rds"), compress = FALSE)
saveRDS(stats_by_boot_summary, file.path(out_dir, "combined_stats_by_boot_summary_with_ICM.rds"), compress = FALSE)
saveRDS(ba_stats, file.path(out_dir, "combined_ba_stats_recalculated_with_ICM.rds"), compress = FALSE)

data.table::fwrite(stats_by_X_method, file.path(out_dir, "combined_stats_by_X_method_with_ICM.csv"))
data.table::fwrite(stats_overall_by_method, file.path(out_dir, "combined_stats_overall_by_method_with_ICM.csv"))
data.table::fwrite(stats_by_boot, file.path(out_dir, "combined_stats_by_boot_with_ICM.csv"))
data.table::fwrite(stats_by_boot_summary, file.path(out_dir, "combined_stats_by_boot_summary_with_ICM.csv"))
data.table::fwrite(ba_stats, file.path(out_dir, "combined_ba_stats_recalculated_with_ICM.csv"))

# -----------------------------
# Plot functions
# -----------------------------
make_ba_plot <- function(which_var) {
  if (!which_var %in% combined_eval$var) {
    message("Skipping BA plot for ", which_var, ": variable not present.")
    return(NULL)
  }
  
  d <- combined_eval[var == which_var]
  b <- ba_stats[var == which_var]
  s <- stats_by_boot_summary[var == which_var]
  s[, label := sprintf("(%.2f, %.2f)", MAD_median, RMSE_median)]  
  lim <- ba_axis_limits[[which_var]]
  xlim_i <- if (!is.null(lim)) lim$x else NULL
  ylim_i <- if (!is.null(lim)) lim$y else NULL
  
  p <- ggplot(d, aes(x = orig, y = diff)) +
    geom_hline(yintercept = 0, linetype = "dotted", colour = "grey50", linewidth = 0.35) +
    geom_bin2d(bins = bins_ba, alpha = 0.90) +
    scale_fill_viridis_c(name = "count", option = "C", trans = "sqrt") +
    geom_hline(
      data = b,
      aes(yintercept = bias),
      inherit.aes = FALSE,
      colour = "steelblue",
      linewidth = 0.45
    ) +
    geom_label(
      data = s,
      aes(x = Inf, y = Inf, label = label),
      inherit.aes = FALSE,
      hjust = 1.08,
      vjust = 1.20,
      size = 3.4,
      label.size = NA,
      fill = "white",
      alpha = 0.75
    ) +
    facet_grid(rows = vars(X_pct), cols = vars(method), scales = "free_x") +
    coord_cartesian(xlim = xlim_i, ylim = ylim_i) +
    labs(
      title = paste0("Bland-Altman comparison with ICM-grid — ", label_var(which_var)),
      subtitle = "Labels show median MAD and median RMSE across bootstrap replicates; blue horizontal line = mean bias.",
      x = x_axis_label_ba(which_var),
      y = y_axis_label_ba(which_var)
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
      plot.subtitle = element_text(hjust = 0.5, size = 12),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 10),
      strip.background = element_rect(fill = "#F0F0F0", colour = NA),
      strip.text = element_text(face = "bold", size = 11),
      panel.grid.minor = element_blank(),
      legend.position = "right"
    )
  
  width  <- max(12, length(method_levels) * 1.35)
  height <- max(7, length(x_order) * 0.75 + 2.5)
  
  outfile_pdf <- file.path(out_dir, paste0("BA_with_ICM_", safe_filename(which_var), ".pdf"))
  outfile_png <- file.path(out_dir, paste0("BA_with_ICM_", safe_filename(which_var), ".png"))
  
  ggsave(outfile_pdf, p, width = width, height = height, units = "in", limitsize = FALSE)
  ggsave(outfile_png, p, width = width, height = height, units = "in", dpi = 300, limitsize = FALSE)
  p
}

make_obs_pred_plot <- function(which_var) {
  if (!which_var %in% combined_eval$var) {
    message("Skipping observed-vs-predicted plot for ", which_var, ": variable not present.")
    return(NULL)
  }
  
  d <- combined_eval[var == which_var]
  panel_stats <- stats_by_X_method[var == which_var]
  panel_stats[, `:=`(
    slope_label = paste0("slope = ", fmt2(slope)),
    ccc_label   = paste0("Lin CCC = ", fmt2(Lin_CCC))
  )]
  
  lim <- obs_pred_axis_limits[[which_var]]
  
  p <- ggplot(d, aes(x = orig, y = pred)) +
    geom_hex(bins = bins_hex, alpha = 0.95) +
    scale_fill_viridis_c(name = "count", option = "C", trans = "sqrt") +
    geom_abline(
      slope = 1,
      intercept = 0,
      linetype = "dashed",
      colour = "red",
      linewidth = 0.45,
      alpha = 0.75
    ) +
    geom_abline(
      data = panel_stats[is.finite(slope) & is.finite(intercept)],
      aes(slope = slope, intercept = intercept),
      inherit.aes = FALSE,
      colour = "steelblue",
      linewidth = 0.60
    ) +
    geom_label(
      data = panel_stats,
      aes(x = -Inf, y = Inf, label = slope_label),
      inherit.aes = FALSE,
      hjust = -0.05,
      vjust = 1.10,
      size = 3.7,
      label.size = NA,
      fill = "white",
      alpha = 0.85,
      label.padding = grid::unit(0.08, "lines")
    ) +
    geom_label(
      data = panel_stats,
      aes(x = Inf, y = -Inf, label = ccc_label),
      inherit.aes = FALSE,
      hjust = 1.05,
      vjust = -0.10,
      size = 3.7,
      label.size = NA,
      fill = "white",
      alpha = 0.85,
      label.padding = grid::unit(0.08, "lines")
    ) +
    facet_grid(rows = vars(X_pct), cols = vars(method)) +
    coord_cartesian(xlim = lim, ylim = lim) +
    labs(
      title = paste0("Observed vs. predicted values with ICM-grid — ", label_var(which_var)),
      subtitle = "Colour scale indicates point density; dashed red line = 1:1 relationship; blue line = fitted regression.",
      x = x_axis_label_obs_pred(which_var),
      y = y_axis_label_obs_pred(which_var)
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
      plot.subtitle = element_text(hjust = 0.5, size = 12),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 10),
      strip.background = element_rect(fill = "#F0F0F0", colour = NA),
      strip.text = element_text(face = "bold", size = 11),
      panel.grid.minor = element_blank(),
      legend.position = "right"
    )
  
  width  <- max(14, length(method_levels) * 1.55)
  height <- max(7.5, length(x_order) * 0.80 + 2.8)
  
  outfile_pdf <- file.path(out_dir, paste0("observed_vs_predicted_with_ICM_", safe_filename(which_var), ".pdf"))
  outfile_png <- file.path(out_dir, paste0("observed_vs_predicted_with_ICM_", safe_filename(which_var), ".png"))
  
  ggsave(outfile_pdf, p, width = width, height = height, units = "in", limitsize = FALSE)
  ggsave(outfile_png, p, width = width, height = height, units = "in", dpi = 300, limitsize = FALSE)
  p
}

make_perf_plot <- function(which_var) {
  if (!which_var %in% stats_by_boot_summary$var) {
    message("Skipping performance plot for ", which_var, ": variable not present.")
    return(NULL)
  }
  
  s <- copy(stats_by_boot_summary[var == which_var])
  s[, X_num := as.numeric(sub("%", "", as.character(X_pct)))]
  
  p_mad <- ggplot(s, aes(x = X_num, y = MAD_median, group = method, colour = method, linetype = method)) +
    geom_line(linewidth = 0.60) +
    geom_point(size = 1.8) +
    labs(
      title = paste0("MAD — ", label_var(which_var)),
      x = "Removed observations (%)",
      y = "Median MAD (‰)",
      colour = "Method",
      linetype = "Method"
    ) +
    theme_bw(base_size = 14) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5))
  
  p_rmse <- ggplot(s, aes(x = X_num, y = RMSE_median, group = method, colour = method, linetype = method)) +
    geom_line(linewidth = 0.60) +
    geom_point(size = 1.8) +
    labs(
      title = paste0("RMSE — ", label_var(which_var)),
      x = "Removed observations (%)",
      y = "Median RMSE (‰)",
      colour = "Method",
      linetype = "Method"
    ) +
    theme_bw(base_size = 14) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5))
  
  p <- (p_mad + p_rmse) + patchwork::plot_layout(guides = "collect") &
    theme(legend.position = "bottom")
  
  outfile_pdf <- file.path(out_dir, paste0("MAD_RMSE_with_ICM_", safe_filename(which_var), ".pdf"))
  outfile_png <- file.path(out_dir, paste0("MAD_RMSE_with_ICM_", safe_filename(which_var), ".png"))
  
  ggsave(outfile_pdf, p, width = 12, height = 5.8, units = "in", limitsize = FALSE)
  ggsave(outfile_png, p, width = 12, height = 5.8, units = "in", dpi = 300, limitsize = FALSE)
  p
}

# -----------------------------
# Create outputs
# -----------------------------
vars_to_plot <- intersect(c("d18O", "d2H", "d-excess"), unique(as.character(combined_eval$var)))

ba_plots <- lapply(vars_to_plot, make_ba_plot)
obs_pred_plots <- lapply(vars_to_plot, make_obs_pred_plot)
perf_plots <- lapply(vars_to_plot, make_perf_plot)

# Print the most manuscript-relevant plots in RStudio.
if (length(ba_plots) > 0 && !is.null(ba_plots[[1]])) print(ba_plots[[1]])
if ("d-excess" %in% vars_to_plot) {
  i_dex <- match("d-excess", vars_to_plot)
  if (!is.null(ba_plots[[i_dex]])) print(ba_plots[[i_dex]])
  if (!is.null(obs_pred_plots[[i_dex]])) print(obs_pred_plots[[i_dex]])
}

message("Done. Combined BA/performance/Lin CCC outputs saved in: ", normalizePath(out_dir))
message("Key input used: ", normalizePath(input_file))


