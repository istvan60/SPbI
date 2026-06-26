# SPbI — Spatial Proximity-Based Imputation

Code repository for:

> Hatvani I.G. & Kern Z. — *Mind the gap: benchmarking imputation methods for stable isotope time series in precipitation* (HESS, 2026)

Benchmarking of eight imputation methods for monthly δ¹⁸O and δ²H time series from stations across Austria, Slovenia, and Hungary (1973–2024). Six common methods (LOCF, linear interpolation, spline, Stineman, Kalman filter, moving average) are compared against a sinusoidal seasonal fit and a novel Spatial Proximity-Based Imputation (SPbI) approach. Performance is evaluated using MAD, RMSE, and Bland–Altman analysis across multiple masking fractions.

---

## Pipeline

```mermaid
flowchart TD
    A([inputation_test6.xlsx\nraw GNIP station data]) --> S0

    S0["**Stage 0** — Data preparation, imputation bootstrap, bundle export\n`Stage_0_prepare_and_bootstrap.R`\n─────────────────────────────────────────\n• Reads δ¹⁸O, δ²H, and station coordinates\n• Filters sites with ≥ 84 consecutive months\n• Applies lapse-rate altitude corrections\n• Runs bootstrap masking across 6 gap fractions\n  (1–32%) with 8 imputation methods in parallel\n• Saves Bland–Altman diagnostic plots\n• Exports reproducibility bundle"]
    S0 --> B([minimal_SPbI_input_bundle.rds])

    B --> S1
    S1["**Stage 1** — SPbI error extraction\n`STAGE 1 PARALLEL.R`\n─────────────────────────────────────────\n• Computes SPbI errors in non-overlapping\n  distance bands around each target station\n• Pairs SPbI errors against Linear and\n  Sinusoidal baseline errors per bootstrap\n• Fully parallelised with future.apply"]
    S1 --> C([spbi_nonoverlapping_band_errors.rds\npaired_SPbI_vs_baselines.rds\nbaseline_errors_site_X_boot_method_isotope.rds])

    C --> S2
    S2["**Stage 2** — CI tests, threshold detection, significance plots\n`STAGE 2 CI detection.R`\n─────────────────────────────────────────\n• Computes paired differences SPbI vs. baselines\n• Derives 95% bootstrap confidence intervals\n• Detects the spatial radius threshold at which\n  SPbI outperforms both baseline methods\n• Produces bandwise boxplots with significance annotations"]
    S2 --> D([paired_difference_CI_and_tests.csv\nthreshold_summary.csv\nsignificance plots .png/.pdf])

    B --> S3
    S3["**Stage 3** — Point-level evaluation export\n`STAGE 3 plotting.R`\n─────────────────────────────────────────\n• Run in the same session as Stage 0\n• Converts the all_imputed object in memory\n  into a structured point-level evaluation table\n• Computes per-method MAD and RMSE summary"]
    S3 --> E([interp_point_eval.rds\ninterp_perf_summary.csv])

    E --> S32
    F([Monthly_194709_202403.nc\nIAEA GNIP ICM NetCDF]) --> S32
    S32["**Stage 3.2** — ICM grid extraction and join\n`STAGE 3.2 ICM extraction.R`\n─────────────────────────────────────────\n• Reads the IAEA GNIP Interpolated Climatology\n  Model NetCDF file\n• Assigns each station to its nearest ICM grid cell\n• Extracts monthly δ¹⁸O and δ²H ICM predictions\n• Joins ICM values to the withheld observations\n  as an additional benchmark method\n• Reports match rate and flags temporal gaps"]
    S32 --> G([combined_point_eval_with_ICM.rds])

    G --> S33
    S33["**Stage 3.3** — Final figures\n`STAGE 3.3. final plotting.R`\n─────────────────────────────────────────\n• Bland–Altman plots for all methods incl. ICM\n• Observed-vs-predicted scatterplots\n• MAD / RMSE performance summary plots\n• All figures exported as .pdf and .png"]
    S33 --> H([BA_with_ICM_*.pdf/png\nobserved_vs_predicted_*.pdf/png\nMAD_RMSE_*.pdf/png])
```

---

## How to run

Run the scripts in order from a fresh R session. Each stage saves its output to disk so subsequent stages can be run independently.

1. **Stage 0** — edit the `input_xlsx` and `out_dir` paths at the top of the script, then source it. Runtime: 30–120 min depending on CPU cores.
2. **Stage 1** — reads `minimal_SPbI_input_bundle.rds`; parallelised with `furrr`.
3. **Stage 2** — reads `spbi_nonoverlapping_band_errors.rds`; produces CI tests and threshold plots.
4. **Stage 3** — must be run in the same R session as Stage 0 (uses the `all_imputed` object in memory); exports the point-level evaluation table.
5. **Stage 3.2** — reads `interp_point_eval.rds` and the IAEA GNIP ICM NetCDF file.
6. **Stage 3.3** — reads `combined_point_eval_with_ICM.rds`; produces all final figures.

## Dependencies

```r
install.packages(c(
  "readxl", "dplyr", "lubridate", "tidyr", "purrr",
  "geosphere", "imputeTS", "minpack.lm",
  "ggplot2", "patchwork", "scales",
  "furrr", "parallelly", "data.table", "ncdf4"
))
```

## License

CC0 1.0 — see [LICENSE](LICENSE).
