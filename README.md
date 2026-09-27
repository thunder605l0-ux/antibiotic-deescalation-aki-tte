# Public R code archive

This directory contains the curated R code needed to reproduce the reported target-trial emulation analyses, diagnostics, tables, and figures. It is intended for public archiving with the manuscript.

## Scope

The archive contains the active scripts listed in `R_script_manifest.csv` (44 R scripts), together with MICE configuration CSVs, aggregate-result CSVs, and package-level audit files:

- `analysis_pipeline/`: data validation, multiple imputation, propensity-score weighting, nested bootstrap, primary outcome estimation, sensitivity analyses, supplementary analyses, table generation, and validation.
- `final_figure_scripts/`: the scripts that generated or finalized the four main figures and Figures S1-S8.
- `analysis_results/`: aggregate ATO AKI sensitivity CSVs that allow `10o --mode=archived` to recompute the submitted confidence intervals without patient-level data.
- `R_script_manifest.csv`: a plain-language inventory of each public script, its internal source, its purpose, and whether local file paths were adapted for public use.
- `analysis_code_index.md`: the manuscript-level mapping from reported analyses to their source scripts and audited SHA-256 hashes.
- `file_manifest_sha256.csv`: the size and SHA-256 digest of every other file in this public archive (the manifest itself is excluded to avoid a recursive self-hash).

One-time diagnostic probes, convergence-investigation scripts, obsolete analysis versions, superseded plots, screenshot/contact-sheet utilities, synchronization-only scripts, and duplicate manuscript copies have been excluded. They remain preserved in the internal project archive and are not required to reproduce the submitted results.

Several files retain `repair` in their historical names because they are part of the final audited analysis chain rather than disposable diagnostics:

- `10k_repair_ato_mortality_sensitivity_2026-08-26.R` generated the final ATO 30-day mortality sensitivity estimates.
- `22_repair_table1_ato_microbiology_2026-08-28.R` added the final microbiology-status rows to the ATO-weighted Table 1 without refitting the propensity-score or outcome models.
- `10n_repair_ato_aki_sensitivity_2026-09-07.R` recalculates two ATO AKI phenotype specifications from the original bootstrap and MICE checkpoints.
- `10o_finalize_ato_aki_sensitivity_2026-09-07.R` validates the saved replicate estimates and recomputes percentile intervals, or performs full checkpoint reconciliation when those files are available.

## Analysis sequence

1. Configure paths and analysis constants in `analysis_pipeline/02_配置与变量字典/analysis_config.R`.
2. Validate the trial datasets with `analysis_pipeline/00_运行说明/01_validate_input_data.R`.
3. Run multiple imputation and MICE diagnostics with the scripts in `analysis_pipeline/03_多重插补/`.
4. Fit primary ATO overlap weights and complementary CBPS-ATE weights with `analysis_pipeline/04_倾向评分与权重/03_fit_cbps_ate_weights.R`, then generate weighting diagnostics.
5. Run the nested patient-level bootstrap using the three execution scripts in `analysis_pipeline/05_Bootstrap/`; pool the locked results with `06b_pool_bootstrap_results_controlled.R`.
6. Generate the primary outcome tables with `analysis_pipeline/06_主要结局分析/07b_primary_outcome_tables_controlled.R`.
7. Run and finalize sensitivity analyses with the scripts in `analysis_pipeline/07_敏感性分析/`. The locked-checkpoint scripts are `10l_run_ato_discharge_censoring_bootstrap_2026-08-30.R`, `10m_run_ato_drug_subgroup_bootstrap_2026-08-30.R`, and `10n_repair_ato_aki_sensitivity_2026-09-07.R`. Script `10o_finalize_ato_aki_sensitivity_2026-09-07.R` assembles or verifies the final ATO AKI results; its default archived mode is described below.
8. Generate supplementary definitions, diagnostics, publication tables, and validation outputs using the retained scripts in `analysis_pipeline/08_supplementary_outputs/`. The historical folder name is retained to preserve the audited source layout; the retained files are the finalized reporting scripts.
9. Run the shared-participant 2-by-2 exploratory analysis in `analysis_pipeline/09_探索性分析/`.
10. Generate the final figures in this order: `plot_bmc_compliance_figures.R`, `plot_english_dual_dag.R`, `plot_figure2_aki_curve_ato.R`, `plot_figure3_death_curve.R`, `plot_main_findings.R`, `plot_primary_ato_diagnostics_2026-08-20.R`, and `plot_refined_supplementary_figures_2_4_8_3.R`. Later scripts intentionally apply the final refinements to Figures S2, S3, S4, and S5-S8.

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

`TTE_PROJECT_ROOT` identifies the research package containing `data/`, `manuscript/`, and the reporting-output folders. `TTE_RUNTIME_ROOT` identifies the analysis runtime. Copy the contents of `analysis_pipeline/` directly into this runtime, preserving its numbered subfolders. The three required MICE configuration files are included in `analysis_pipeline/02_配置与变量字典/`: `mice_working_variable_dictionary.csv`, `mice_method_specification.csv`, and `mice_predictor_matrix_directed_uo.csv`. They are byte-identical to the locked analysis configurations.

By default, input datasets are read from `TTE_RUNTIME_ROOT/01_原始分析数据/` and derived objects from `TTE_RUNTIME_ROOT/09_输出结果/`. Set `TTE_DATA_DIR` and `TTE_OUTPUT_ROOT` to use existing directories without copying patient-level data. Install the locked package versions in the configuration's `R_library/` locations, or set `TTE_R_LIBRARY_ROOT` to an existing library root containing the `mice-github-61083e667fa41cc4bc655492591b100bf8a7e443/` and `phase47-weightit-cobalt-github/` libraries. Other dependencies must be available on `.libPaths()`. The archive supplies code and configuration, not an installed R library or the complete private runtime.

### Recompute ATO AKI confidence intervals from the public archive

This mode requires only base R and the aggregate CSVs already included in this archive. It checks replicate identities, counts, risk ranges, RD/RR arithmetic, and all type-7 percentile confidence limits against the submitted results. It also checks the two phenotype point estimates against their saved original-sample outputs. Run from the extracted `R/` directory:

```text
Rscript --vanilla "analysis_pipeline/07_敏感性分析/10o_finalize_ato_aki_sensitivity_2026-09-07.R" --mode=archived --output=ATO_AKI_sensitivity_recomputed.csv
```

The input CSVs remain unchanged. This is an aggregate-result verification, not a patient-level model refit.

### Recalculate the two ATO AKI definitions from original checkpoints

Script `10n` additionally requires the following authorized, patient-level inputs and derived objects:

- In `TTE_DATA_DIR`: `analysis_dataset_anti_mrsa.csv`, `analysis_dataset_anti_psa.csv`, and `aki_three_path_patient_summary.csv` (with `trial`, `stay_id`, `first_creatinine_aki_time`, `renal_followup_end`, `creatinine_aki`, and `followup_end_type`).
- In `TTE_SUBMISSION_ROOT`: `final.data_anti-MRSA.csv` and `final.data_anti-PSA.csv`, used to verify that the runtime datasets match the submission data version. This directory defaults to `TTE_PROJECT_ROOT/final_data`.
- In `TTE_OUTPUT_ROOT/02_MICE/`: `mids_anti_mrsa_m20_maxit50.rds` and `mids_anti_psa_m20_maxit50.rds`.
- In `TTE_OUTPUT_ROOT/04_Bootstrap/checkpoints/`: the original `<trial>_bNNNN_checkpoint.rds` files containing the sampled-row keys, MICE objects, and trial/replicate identity. The antipseudomonal failed replicate 934 is excluded, not replaced.
- In `TTE_OUTPUT_ROOT/06_敏感性分析/ato_mortality_repair_20260826/full/`: `ST-08_ato_mortality_sensitivity_bootstrap_replicates.csv` and `ST-08_ato_mortality_sensitivity_point_estimates.csv`, generated by `10k`. These outputs also contain the four reusable ATO AKI specifications and the primary AKI estimates under the `death_calendar_boundary` label.

Set `TTE_ATO_AKI_OUTPUT` to an empty writable directory for recalculated outputs. For a one-checkpoint validation in each trial, run:

```text
Rscript --vanilla "analysis_pipeline/07_敏感性分析/10n_repair_ato_aki_sensitivity_2026-09-07.R" --trial=anti_mrsa --start=1 --end=1 --preflight=true
Rscript --vanilla "analysis_pipeline/07_敏感性分析/10n_repair_ato_aki_sensitivity_2026-09-07.R" --trial=anti_psa --start=1 --end=1 --preflight=true
```

For all saved successful checkpoints, use `--start=1 --end=1000` for each trial, then run `10o` with `--mode=checkpoints` and the same environment variables. No new bootstrap samples or MICE datasets are generated by `10n`. Full checkpoint mode checks the recalculated outputs against the original checkpoint hashes; it cannot run from aggregate CSVs alone.

MIMIC-IV access remains subject to PhysioNet credentialing and its data-use agreement; this code archive does not grant access to the source database.

## Submission additional-file mapping

The submission package uses descriptive filenames beginning with `Additional file`, followed by the manuscript label (`Table S1-S11` or `Figure S1-S8`) and a short content description. Supplementary methods and the estimand amendment are in `Additional file Supplementary methods.docx`. Script-internal ST/SF identifiers remain unchanged because they are part of the locked analysis and audit trail.

| Submission filename | Manuscript label |
|---|---|
| `Additional file Supplementary methods.docx` | Supplementary methods, Sections 1-4 |
| `Additional file Table S1 Antibacterial exposure dictionary.docx` | Table S1 |
| `Additional file Table S2 Target trial emulation mapping.docx` | Table S2 |
| `Additional file Table S3 AKI phenotype components, KDIGO stages, and urine-output missingness.docx` | Table S3 |
| `Additional file Figure S1 Full causal and measured operational DAGs.pdf` | Figure S1 |
| `Additional file Table S4 Variable definitions and coding.docx` | Table S4 |
| `Additional file Figure S2 Missing data by trial and treatment strategy.pdf` | Figure S2 |
| `Additional file Figure S3 Spearman correlations among continuous variables.pdf` | Figure S3 |
| `Additional file Figure S4 MICE chain means across 20 imputations.pdf` | Figure S4 |
| `Additional file Table S5 Primary ATO weighting diagnostics.docx` | Table S5 |
| `Additional file Figure S5 Covariate balance before and after primary ATO weighting.pdf` | Figure S5 |
| `Additional file Figure S6 Propensity-score overlap for primary ATO analysis.pdf` | Figure S6 |
| `Additional file Figure S7 Distribution of primary ATO overlap weights.pdf` | Figure S7 |
| `Additional file Table S6 Complementary CBPS-ATE estimand results.docx` | Table S6 |
| `Additional file Figure S8 AKI analyses under CBPS-ATE and ATO estimands.pdf` | Figure S8 |
| `Additional file Table S7 Primary ATO mortality sensitivity analyses.docx` | Table S7 |
| `Additional file Table S8 E-values for unmeasured confounding.docx` | Table S8 |
| `Additional file Table S9 Unweighted baseline characteristics.docx` | Table S9 |
| `Additional file Table S10 Exploratory drug-subgroup AKI analyses.docx` | Table S10 |
| `Additional file Table S11 Joint strategy cross-classification and AKI risk.docx` | Table S11 |

## Code checks

All included R files must pass `parse()` under R 4.4.2. `R_script_manifest.csv` explains what each script does and identifies scripts for which local absolute paths were replaced by `TTE_PROJECT_ROOT` or `TTE_RUNTIME_ROOT`. `file_manifest_sha256.csv` provides the package-wide integrity check. These portability changes did not alter the analysis methods or statistical settings.
