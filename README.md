# Public R code archive

This directory contains the curated R code needed to reproduce the reported target-trial emulation analyses, diagnostics, tables, and figures. It is intended for public archiving with the manuscript.

## Scope

The archive contains 40 active scripts in three groups:

- `analysis_pipeline/`: data validation, multiple imputation, propensity-score weighting, nested bootstrap, primary outcome estimation, sensitivity analyses, supplementary analyses, table generation, and validation.
- `final_figure_scripts/`: the scripts that generated or finalized the four main figures and Supplementary Figures 1-8.
- `R_script_manifest.csv`: a plain-language inventory of each public script, its internal source, its purpose, and whether local file paths were adapted for public use.

One-time diagnostic probes, convergence-investigation scripts, obsolete analysis versions, superseded plots, screenshot/contact-sheet utilities, synchronization-only scripts, and duplicate manuscript copies have been excluded. They remain preserved in the internal project archive and are not required to reproduce the submitted results.

Two files retain `repair` in their historical names because they are part of the final audited analysis chain rather than disposable diagnostics:

- `10k_repair_ato_mortality_sensitivity_2026-08-26.R` generated the final ATO 30-day mortality sensitivity estimates.
- `22_repair_table1_ato_microbiology_2026-08-28.R` added the final microbiology-status rows to the ATO-weighted Table 1 without refitting the propensity-score or outcome models.

## Analysis sequence

1. Configure paths and analysis constants in `analysis_pipeline/02_配置与变量字典/analysis_config.R`.
2. Validate the trial datasets with `analysis_pipeline/00_运行说明/01_validate_input_data.R`.
3. Run multiple imputation and MICE diagnostics with the scripts in `analysis_pipeline/03_多重插补/`.
4. Fit primary ATO overlap weights and complementary CBPS-ATE weights with `analysis_pipeline/04_倾向评分与权重/03_fit_cbps_ate_weights.R`, then generate weighting diagnostics.
5. Run the nested patient-level bootstrap using the three execution scripts in `analysis_pipeline/05_Bootstrap/`; pool the locked results with `06b_pool_bootstrap_results_controlled.R`.
6. Generate the primary outcome tables with `analysis_pipeline/06_主要结局分析/07b_primary_outcome_tables_controlled.R`.
7. Run and finalize sensitivity analyses with the scripts in `analysis_pipeline/07_敏感性分析/`.
8. Generate supplementary definitions, diagnostics, publication tables, and validation outputs using the retained scripts in `analysis_pipeline/08_未完成图表/`. The historical folder name is retained to preserve the audited source layout; the retained files are the finalized reporting scripts.
9. Run the shared-participant 2-by-2 exploratory analysis in `analysis_pipeline/09_探索性分析/`.
10. Generate the final figures in this order: `plot_bmc_compliance_figures.R`, `plot_english_dual_dag.R`, `plot_figure2_aki_curve_ato.R`, `plot_figure3_death_curve.R`, `plot_main_findings.R`, `plot_primary_ato_diagnostics_2026-08-20.R`, and `plot_refined_supplementary_figures_2_4_8_3.R`. Later scripts intentionally apply the final refinements to Supplementary Figures 2, 3, 4, 5-8.

## Reproducibility notes

- R version: 4.4.2.
- Multiple imputation: 20 completed datasets, 50 iterations, 5 predictive-mean-matching donors where applicable.
- Bootstrap: 1,000 patient-level resamples per trial; the anti-MRSA analysis retained 1,000 successful resamples and the antipseudomonal analysis retained 999 successful resamples, with failed replicate `b = 934` not replaced.
- Primary estimand: ATO using overlap weights from a bias-reduced logistic propensity-score model.
- Complementary estimand: CBPS-ATE.
- Random seeds and deterministic seed derivation are defined in the configuration and bootstrap scripts.

Before running the scripts, define two environment variables:

```r
Sys.setenv(
  TTE_PROJECT_ROOT = "/path/to/research-project",
  TTE_RUNTIME_ROOT = "/path/to/independent-R-runtime"
)
```

`TTE_PROJECT_ROOT` identifies the research package containing `data/`, `manuscript/`, and the reporting-output folders. `TTE_RUNTIME_ROOT` identifies the independent analysis runtime containing the locked input datasets, configuration files, package library, and derived-output directories. The public copies use these variables instead of the authors' local absolute paths. Preserve the documented directory structure when running the scripts.

MIMIC-IV access remains subject to PhysioNet credentialing and its data-use agreement; this code archive does not grant access to the source database.

## Code checks

All included R files must pass `parse()` under R 4.4.2. The manifest explains what each script does and identifies scripts for which local absolute paths were replaced by `TTE_PROJECT_ROOT` or `TTE_RUNTIME_ROOT`. These portability changes did not alter the analysis methods or statistical settings.
