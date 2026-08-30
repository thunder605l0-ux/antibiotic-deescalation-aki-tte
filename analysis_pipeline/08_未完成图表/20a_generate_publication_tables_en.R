# Analysis: Generate English publication-ready tables from locked controlled sources
# Date: 2026-08-08
# Random seed: not applicable (presentation-only transformation)
# R: 4.4.2
#
# Method: deterministic column projection, reader-label mapping, and three-rule DOCX rendering
# Source: locked MT-01 to MT-03 and ST-01 to ST-10 CSV outputs
# Citation/version: RD-4.6-ANALYSIS-OUTPUT-PLAN-V1.0; publication display dictionary v1.0
# Access date: 2026-08-08
# Adaptation: no model is fitted and no estimate is recomputed; machine fields remain in source/audit files.
# Verification command: D:/software/R4.4.2/R-4.4.2/bin/x64/Rscript.exe 08_未完成图表/20a_generate_publication_tables_en.R

options(stringsAsFactors = FALSE)

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg[[1L]]), winslash = "/", mustWork = TRUE)
project_root <- dirname(dirname(script_path))
source(file.path(project_root, "02_配置与变量字典", "publication_style_helpers.R"), encoding = "UTF-8")

if (!requireNamespace("flextable", quietly = TRUE)) stop("Package 'flextable' is required.")
if (!requireNamespace("officer", quietly = TRUE)) stop("Package 'officer' is required.")
if (!requireNamespace("digest", quietly = TRUE)) stop("Package 'digest' is required.")

output_root <- file.path(project_root, "09_输出结果", "11_publication_ready_en")
main_dir <- file.path(output_root, "02_main_tables")
supp_dir <- file.path(output_root, "04_supplementary_tables")
audit_dir <- file.path(output_root, "_audit")
dir.create(main_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(supp_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)

dictionary <- load_publication_dictionary(project_root)
display_by_field <- setNames(dictionary$display_label, dictionary$clean_field)
unit_by_field <- setNames(dictionary$unit, dictionary$clean_field)

read_locked_csv <- function(relative_path) {
  path <- file.path(project_root, relative_path)
  if (!file.exists(path)) stop("Missing locked table source: ", path)
  utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE, fileEncoding = "UTF-8")
}

trial_label <- function(x) {
  map <- c(anti_mrsa = "Anti-MRSA trial", anti_psa = "Anti-PSA trial")
  out <- unname(map[x])
  if (anyNA(out)) stop("Unmapped trial code: ", paste(unique(x[is.na(out)]), collapse = ", "))
  out
}

model_label <- function(x) {
  map <- c(
    cbps_ate = "CBPS-ATE",
    ato = "Overlap weighting (ATO)",
    ato_brlogit = "Overlap weighting (ATO; bias-reduced logistic regression)"
  )
  out <- unname(map[x])
  if (anyNA(out)) stop("Unmapped model code: ", paste(unique(x[is.na(out)]), collapse = ", "))
  out
}

specification_labels <- c(
  main_cbps_ate = "Main analysis: CBPS-ATE",
  overlap_weight_ato_brlogit = "Overlap weighting with bias-reduced logistic regression",
  aki_cross_t0_uo_boundary = "Urine-output window allowed to cross index time",
  death_calendar_boundary = "Alternative death-date boundary",
  weight_truncation_1_99 = "ATE weights truncated at the 1st and 99th percentiles",
  weight_truncation_5_95 = "ATE weights truncated at the 5th and 95th percentiles",
  creatinine_only_phenotype = "Creatinine-only AKI phenotype",
  NEW_MS03 = "Alternative DAG adjustment set including Race",
  ABX_TIME_proxy = "Adjustment including the antibiotic-timing proxy",
  MICRO_STATE_proxy = "Microbiology-status proxy perturbation",
  complete_case = "Complete-case analysis"
)

estimand_labels <- c(
  risk_aki_0 = "Continuation: 7-day AKI risk",
  risk_aki_1 = "De-escalation: 7-day AKI risk",
  rd_aki = "Risk difference for 7-day AKI",
  rr_aki = "Risk ratio for 7-day AKI",
  risk_death30_0 = "Continuation: 30-day all-cause mortality risk",
  risk_death30_1 = "De-escalation: 30-day all-cause mortality risk",
  rd_death30 = "Risk difference for 30-day all-cause mortality",
  rr_death30 = "Risk ratio for 30-day all-cause mortality"
)

title_map <- c(
  "MT-01" = "Weighted baseline characteristics after CBPS-ATE weighting",
  "MT-02" = "Outcomes in the anti-MRSA trial",
  "MT-03" = "Outcomes in the anti-PSA trial",
  "ST-01" = "Target trial specification and MIMIC-IV emulation",
  "ST-02" = "Antibacterial classification and route dictionary",
  "ST-03" = "Variable definitions, coding, and analysis handling",
  "ST-04" = "Unweighted baseline characteristics",
  "ST-05" = "Multiple-imputation specification and convergence diagnostics",
  "ST-06" = "Propensity-score and weight diagnostics",
  "ST-07" = "Prespecified sensitivity analyses",
  "ST-08" = "E-values for unmeasured confounding",
  "ST-09" = "Components of the 7-day AKI phenotype"
)

write_csv <- function(data, path) {
  utils::write.csv(data, path, row.names = FALSE, na = "", fileEncoding = "UTF-8")
}

add_tnr_paragraph <- function(doc, text, bold = FALSE, size = 10, align = "left") {
  officer::body_add_fpar(
    doc,
    officer::fpar(
      officer::ftext(text, officer::fp_text(font.family = publication_font_family, font.size = size, bold = bold)),
      fp_p = officer::fp_par(text.align = align)
    )
  )
}

make_docx <- function(data, path, title, note = NULL, subtables = NULL, font_size = 8) {
  doc <- officer::read_docx()
  doc <- add_tnr_paragraph(doc, title, bold = TRUE, size = 11, align = "center")
  if (is.null(subtables)) subtables <- list(data)
  for (i in seq_along(subtables)) {
    current <- subtables[[i]]
    subtitle <- attr(current, "subtitle", exact = TRUE)
    if (!is.null(subtitle) && nzchar(subtitle)) doc <- add_tnr_paragraph(doc, subtitle, bold = TRUE, size = 9)
    ft <- publication_three_rule_table(current, font_size = font_size)
    ft <- flextable::set_table_properties(ft, layout = "autofit", opts_word = list(repeat_headers = TRUE, split = FALSE))
    doc <- flextable::body_add_flextable(doc, value = ft)
    if (i < length(subtables)) doc <- officer::body_add_break(doc)
  }
  if (!is.null(note) && nzchar(note)) {
    doc <- add_tnr_paragraph(doc, paste0("Note: ", note), size = 8)
  }
  section <- officer::prop_section(
    page_size = officer::page_size(orient = "landscape"),
    page_margins = officer::page_mar(top = 0.45, bottom = 0.45, left = 0.45, right = 0.45)
  )
  doc <- officer::body_end_block_section(doc, value = officer::block_section(section))
  print(doc, target = path)
  invisible(path)
}

assert_reader_ready <- function(data, id) {
  character_values <- unlist(data[vapply(data, is.character, logical(1))], use.names = FALSE)
  character_values <- character_values[nzchar(character_values)]
  cjk <- grepl("[\u3400-\u9FFF]", character_values, perl = TRUE)
  forbidden <- grepl("race_white|White race|factor\\(|anti_mrsa|anti_psa|cbps_ate|ato_brlogit|NEW_MS0[1-9]", character_values) |
    grepl("[A-Za-z0-9]+_[A-Za-z0-9_]+", character_values) |
    grepl("MICE:|PS formula|target-antibiotic", character_values, ignore.case = TRUE)
  if (any(cjk)) stop(id, " contains Chinese text in a formal publication table: ", character_values[which(cjk)[1L]])
  if (any(forbidden)) stop(id, " exposes a forbidden machine/readability label: ", character_values[which(forbidden)[1L]])
  invisible(TRUE)
}

sanitize_reader_text <- function(data) {
  replacements <- c(
    "MICE: pmm" = "Multiple imputation using predictive mean matching",
    "MICE: logreg" = "Multiple imputation using logistic regression",
    "MICE: polr" = "Multiple imputation using proportional-odds logistic regression",
    "MICE auxiliary variable" = "Multiple-imputation auxiliary variable",
    "PS formula" = "propensity-score formula",
    "Target-antibiotic" = "Target antibiotic",
    "target-antibiotic" = "target antibiotic",
    "NEW_MS02/NEW_MS04" = "the secondary sensitivity analyses",
    "NEW_MS01" = "the main prespecified adjustment set",
    "NEW_MS02" = "a secondary sensitivity analysis",
    "NEW_MS03" = "the alternative DAG adjustment set",
    "NEW_MS04" = "a secondary sensitivity analysis",
    "ABX_TIME_proxy" = "the antibiotic-timing proxy adjustment",
    "ABX_TIME" = "the antibiotic-timing proxy",
    "MICRO_STATE_proxy" = "the microbiology-status proxy perturbation",
    "diagnoses_icd" = "diagnosis code table",
    "cbps_ate" = "CBPS-ATE", "ato_brlogit" = "overlap weighting with bias-reduced logistic regression",
    "anti_mrsa" = "Anti-MRSA", "anti_psa" = "Anti-PSA", "race_white" = "Race", "White race" = "Race"
  )
  for (column in names(data)) {
    if (!is.character(data[[column]])) next
    for (from in names(replacements)) data[[column]] <- gsub(from, replacements[[from]], data[[column]], fixed = TRUE)
    data[[column]] <- gsub("T0", "index time", data[[column]], fixed = TRUE)
  }
  data
}

format_publication_number <- function(x, digits = 3L) {
  out <- rep("", length(x))
  finite <- is.finite(x)
  out[finite] <- formatC(x[finite], format = "f", digits = digits)
  out[is.infinite(x) & x > 0] <- "Inf"
  out[is.infinite(x) & x < 0] <- "-Inf"
  out
}

compact_outcome_table <- function(data) {
  ci_text <- function(point, lower, upper) {
    paste0(
      format_publication_number(point), " (",
      format_publication_number(lower), " to ",
      format_publication_number(upper), ")"
    )
  }
  data.frame(
    `Outcome` = data$Outcome,
    `Model` = data$Model,
    `Continuation events/N` = paste0(data$`Continuation events`, "/", data$`Continuation N`),
    `Continuation risk (95% CI)` = ci_text(data$`Continuation risk`, data$`Continuation 95% CI lower`, data$`Continuation 95% CI upper`),
    `De-escalation events/N` = paste0(data$`De-escalation events`, "/", data$`De-escalation N`),
    `De-escalation risk (95% CI)` = ci_text(data$`De-escalation risk`, data$`De-escalation 95% CI lower`, data$`De-escalation 95% CI upper`),
    `Risk difference (95% CI)` = ci_text(data$`Risk difference`, data$`RD 95% CI lower`, data$`RD 95% CI upper`),
    `Risk ratio (95% CI)` = ci_text(data$`Risk ratio`, data$`RR 95% CI lower`, data$`RR 95% CI upper`),
    `Successful replicates` = data$`Successful bootstrap replicates`,
    check.names = FALSE
  )
}

outputs <- list()
sources <- list()

# MT-01
src <- read_locked_csv("09_输出结果/09_tables/MT-01_weighted_baseline_characteristics.csv")
characteristic <- src$characteristic
map_rows <- nzchar(src$field) & src$row_type %in% c("continuous", "binary", "categorical_header")
mt01_field <- src$field
mt01_field[mt01_field == "target_class_grouped"] <- "initial_target_antibiotic"
characteristic[map_rows] <- unname(display_by_field[mt01_field[map_rows]])
if (anyNA(characteristic[map_rows])) stop("MT-01 contains an unmapped field.")
out <- data.frame(
  `Trial` = trial_label(src$trial), `Section` = src$section,
  `Characteristic` = characteristic, `De-escalation` = src$deescalated,
  `Continuation` = src$continued, `Maximum absolute SMD` = src$maximum_absolute_smd,
  `Variance-ratio range` = src$variance_ratio_range, check.names = FALSE
)
outputs[["MT-01"]] <- out; sources[["MT-01"]] <- src

# MT-02 and MT-03
for (id in c("MT-02", "MT-03")) {
  src <- read_locked_csv(paste0(
    "09_输出结果/05_主要结局/controlled_reporting/", id,
    ifelse(id == "MT-02", "_anti_mrsa", "_anti_psa"), "_outcomes_controlled.csv"
  ))
  out <- data.frame(
    `Outcome` = src$outcome, `Model` = model_label(src$model),
    `Continuation events` = src$continuation_events, `Continuation N` = src$continuation_n,
    `Continuation risk` = src$continuation_risk, `Continuation 95% CI lower` = src$continuation_ci_lower,
    `Continuation 95% CI upper` = src$continuation_ci_upper,
    `De-escalation events` = src$deescalation_events, `De-escalation N` = src$deescalation_n,
    `De-escalation risk` = src$deescalation_risk, `De-escalation 95% CI lower` = src$deescalation_ci_lower,
    `De-escalation 95% CI upper` = src$deescalation_ci_upper,
    `Risk difference` = src$risk_difference, `RD 95% CI lower` = src$rd_ci_lower,
    `RD 95% CI upper` = src$rd_ci_upper, `Risk ratio` = src$risk_ratio,
    `RR 95% CI lower` = src$rr_ci_lower, `RR 95% CI upper` = src$rr_ci_upper,
    `Successful bootstrap replicates` = src$bootstrap_success_n,
    `Confidence level` = src$confidence_level, `Quantile type` = src$quantile_type,
    check.names = FALSE
  )
  outputs[[id]] <- out; sources[[id]] <- src
}

# ST-01
src <- read_locked_csv("09_输出结果/09_tables/ST-01_target_trial_mapping.csv")
out <- src
names(out) <- c("Element", "Hypothetical target trial", "MIMIC-IV emulation")
outputs[["ST-01"]] <- out; sources[["ST-01"]] <- src

# ST-02: publication projection excludes internal Chinese audit columns.
src <- read_locked_csv("09_输出结果/09_tables/ST-02_antibacterial_dictionary.csv")
classification_map <- c(
  eligible_systemic_antibacterial = "Eligible systemic antibacterial",
  local_or_non_systemic_product = "Local or non-systemic product",
  non_bacterial_antimicrobial = "Non-bacterial antimicrobial",
  non_systemic_route = "Non-systemic route",
  non_treatment_or_placeholder = "Non-treatment or placeholder",
  nonabsorbed_enteral_antibacterial = "Nonabsorbed enteral antibacterial",
  not_an_antibacterial = "Not an antibacterial",
  oral_enteral_vancomycin = "Oral or enteral vancomycin"
)
route_map <- c(
  catheter_local = "Catheter-local", cns_local = "Central nervous system-local",
  enteral_systemic = "Enteral systemic", inhaled_local = "Inhaled local",
  irrigation_local = "Irrigation local", ophthalmic_local = "Ophthalmic local",
  order_route_proxy = "Order-route proxy", otic_local = "Otic local",
  parenteral_systemic = "Parenteral systemic", peritoneal_local = "Peritoneal local",
  rectal_local = "Rectal local", topical_local = "Topical local"
)
classification_reason <- ifelse(src$classification_reason %in% names(classification_map), unname(classification_map[src$classification_reason]), src$classification_reason)
route_group <- ifelse(src$route_group %in% names(route_map), unname(route_map[src$route_group]), src$route_group)
out <- data.frame(
  `Section` = src$table_section,
  `Trial class` = ifelse(src$trial_class == "anti_mrsa", "Anti-MRSA", ifelse(src$trial_class == "anti_psa", "Anti-PSA", src$trial_class)),
  `Canonical generic name` = src$canonical_generic,
  `Name in source publication` = src$source_name_in_paper_a,
  `Alias pattern` = src$alias_pattern,
  `Detected in MIMIC-IV eMAR` = src$mimic_emar_raw_detected,
  `Raw drug name` = src$raw_drug_name,
  `Raw route` = src$raw_route,
  `Included in main analysis` = src$include_main,
  `Classification rationale` = classification_reason,
  `Source rows` = src$source_rows,
  `Database` = src$database,
  `Route group` = route_group,
  `Eligible systemic coverage` = src$eligible_for_systemic_sepsis_coverage,
  check.names = FALSE
)
outputs[["ST-02"]] <- out; sources[["ST-02"]] <- src

# ST-03
src <- read_locked_csv("09_输出结果/09_tables/ST-03_variable_definitions.csv")
labels <- unname(display_by_field[src$clean_field])
if (anyNA(labels)) stop("ST-03 contains an unmapped variable.")
main_status <- gsub("NEW_MS03", "the alternative DAG sensitivity analysis", src$main_model_status, fixed = TRUE)
main_status <- gsub("SEX", "Male sex", main_status, fixed = TRUE)
missing_handling <- gsub("NEW_MS01", "the main analysis", src$missing_data_handling, fixed = TRUE)
missing_handling <- gsub("NEW_MS03", "the alternative DAG sensitivity analysis", missing_handling, fixed = TRUE)
measurement_window <- gsub("T0", "index time", src$measurement_window, fixed = TRUE)
out <- data.frame(
  `No.` = src$number, `Variable` = labels, `Variable type` = src$variable_type,
  `Unit` = ifelse(src$unit == "None", "", src$unit), `Coding or definition` = src$numeric_code_definition,
  `Measurement window` = measurement_window,
  `Main-model status` = main_status, `Missing-data handling` = missing_handling,
  `Data source` = src$data_source, check.names = FALSE
)
outputs[["ST-03"]] <- out; sources[["ST-03"]] <- src

# ST-04
src <- read_locked_csv("09_输出结果/09_tables/controlled_reporting/ST-04_unweighted_baseline_characteristics_controlled.csv")
characteristic <- src$characteristic
characteristic[characteristic == "White race"] <- "Race"
out <- data.frame(
  `Trial` = trial_label(src$trial), `Section` = src$section, `Characteristic` = characteristic,
  `Overall` = src$overall, `De-escalation` = src$deescalated, `Continuation` = src$continued,
  `Absolute SMD` = src$absolute_smd, check.names = FALSE
)
outputs[["ST-04"]] <- out; sources[["ST-04"]] <- src

# ST-05
src <- read_locked_csv("09_输出结果/09_tables/ST-06_multiple_imputation_specification_QA.csv")
out <- data.frame(
  `Trial` = trial_label(src$trial), `N` = src$n, `Imputations` = src$m,
  `Maximum iterations` = src$maxit, `Predictive mean-matching donors` = src$donors,
  `Seed` = src$seed, `Directed urine-output order` = src$directed_urine_order,
  `Active imputation targets` = src$active_imputation_targets,
  `Logged events` = src$logged_events, `Maximum final PSRF` = src$maximum_final_psrf,
  `Variables with final PSRF >1.10` = src$variables_with_final_psrf_above_1_10,
  `Residual missing cells in main propensity-score model` = src$residual_missing_cells_in_main_ps,
  check.names = FALSE
)
outputs[["ST-05"]] <- out; sources[["ST-05"]] <- src

# ST-06
src <- read_locked_csv("09_输出结果/09_tables/ST-07_weight_diagnostics.csv")
out <- src
out$trial <- trial_label(out$trial)
out$model <- model_label(out$model)
names(out) <- c(
  "Trial", "Model", "Imputations", "Minimum propensity score", "Maximum propensity score",
  "Maximum N with propensity score <0.01", "Maximum N with propensity score >0.99",
  "Maximum de-escalation N with propensity score <0.01", "Maximum continuation N with propensity score >0.99",
  "Maximum fraction outside 0.01–0.99", "Minimum empirical overlap width",
  "Empirical overlap absent in any imputation", "Minimum weight", "Minimum 1st percentile of weights",
  "Minimum 5th percentile of weights", "Minimum median weight", "Maximum median weight",
  "Maximum 95th percentile of weights", "Maximum 99th percentile of weights", "Maximum weight",
  "Maximum N with weight >10", "Maximum fraction with weight >10", "Minimum continuation ESS",
  "Minimum de-escalation ESS", "Minimum continuation ESS fraction", "Minimum de-escalation ESS fraction",
  "Maximum absolute SMD"
)
outputs[["ST-06"]] <- out; sources[["ST-06"]] <- src

# ST-07
src <- read_locked_csv("09_输出结果/06_敏感性分析/controlled_reporting/ST-08_prespecified_sensitivity_results.csv")
spec <- unname(specification_labels[src$specification])
estimand <- unname(estimand_labels[src$estimand_field])
if (anyNA(spec) || anyNA(estimand)) stop("ST-07 contains an unmapped specification or estimand.")
out <- data.frame(
  `Trial` = trial_label(src$trial), `Sensitivity analysis` = spec, `Estimand` = estimand,
  `Successful bootstrap replicates` = src$bootstrap_B, `95% CI lower` = src$ci_lower,
  `Bootstrap median` = src$bootstrap_median, `95% CI upper` = src$ci_upper,
  `Point estimate` = src$point_estimate, `Analysis N` = src$n_analysis,
  `Truncated fraction` = src$truncated_fraction,
  `Minimum continuation ESS` = src$ess_continuation_min,
  `Minimum de-escalation ESS` = src$ess_deescalation_min,
  `Risk set rebuilt` = src$risk_set_rebuilt, `Estimate (95% CI)` = src$estimate_95ci,
  check.names = FALSE
)
outputs[["ST-07"]] <- out; sources[["ST-07"]] <- src

# ST-08
src <- read_locked_csv("09_输出结果/06_敏感性分析/controlled_reporting/ST-09_unmeasured_confounding_evalue.csv")
out <- data.frame(
  `Trial` = trial_label(src$trial), `Outcome` = src$outcome, `Model` = model_label(src$model),
  `Estimand` = src$estimand, `Risk ratio` = src$risk_ratio, `RR 95% CI lower` = src$rr_ci_lower,
  `RR 95% CI upper` = src$rr_ci_upper, `Confidence level` = src$confidence_level,
  `Quantile type` = src$quantile_type, `Successful bootstrap replicates` = src$bootstrap_success_n,
  `CI limit closest to the null` = src$ci_limit_closest_to_null,
  `CI includes the null` = src$ci_includes_null, `E-value for point estimate` = src$evalue_point,
  `E-value for CI limit` = src$evalue_ci_limit, `Interpretation` = src$interpretation,
  check.names = FALSE
)
outputs[["ST-08"]] <- out; sources[["ST-08"]] <- src

# ST-09
src <- read_locked_csv("09_输出结果/06_敏感性分析/ST-10_aki_phenotype_components.csv")
strategy <- c(`0` = "Continuation", `1` = "De-escalation")[as.character(src$treatment_strategy)]
if (anyNA(strategy)) stop("ST-09 contains an unmapped treatment strategy.")
out <- data.frame(
  `Trial` = trial_label(src$trial), `Treatment strategy` = unname(strategy), `N` = src$n,
  `Any AKI` = src$aki_any_n, `Creatinine pathway` = src$creatinine_path_n,
  `Urine-output pathway` = src$urine_output_path_n, `Renal replacement therapy pathway` = src$rrt_path_n,
  `Creatinine only` = src$creatinine_only_n, `Urine output only` = src$urine_only_n,
  `Renal replacement therapy only` = src$rrt_only_n, `Creatinine and urine output` = src$creatinine_urine_n,
  `Creatinine and renal replacement therapy` = src$creatinine_rrt_n,
  `Urine output and renal replacement therapy` = src$urine_rrt_n,
  `All three pathways` = src$all_three_n, check.names = FALSE
)
outputs[["ST-09"]] <- out; sources[["ST-09"]] <- src

expected_ids <- c(sprintf("MT-%02d", 1:3), sprintf("ST-%02d", 1:9))
stopifnot(setequal(names(outputs), expected_ids))

manifest_rows <- list()
qa_rows <- list()
for (id in expected_ids) {
  out <- sanitize_reader_text(outputs[[id]])
  outputs[[id]] <- out
  assert_reader_ready(out, id)
  target_dir <- if (startsWith(id, "MT")) main_dir else supp_dir
  stem <- c(
    "MT-01" = "MT-01_weighted_baseline_characteristics", "MT-02" = "MT-02_anti_mrsa_outcomes",
    "MT-03" = "MT-03_anti_psa_outcomes", "ST-01" = "ST-01_target_trial_mapping",
    "ST-02" = "ST-02_antibacterial_dictionary", "ST-03" = "ST-03_variable_definitions",
    "ST-04" = "ST-04_unweighted_baseline_characteristics", "ST-05" = "ST-05_multiple_imputation_specification",
    "ST-06" = "ST-06_weight_diagnostics", "ST-07" = "ST-07_prespecified_sensitivity_results",
    "ST-08" = "ST-08_unmeasured_confounding_evalue", "ST-09" = "ST-09_aki_phenotype_components"
  )[[id]]
  csv_path <- file.path(target_dir, paste0(stem, ".csv"))
  docx_path <- file.path(target_dir, paste0(stem, ".docx"))
  write_csv(out, csv_path)
  note <- if (id %in% c("MT-02", "MT-03", "ST-07", "ST-08")) {
    "Confidence intervals are two-sided 95% percentile intervals calculated with R quantile type 7 from the locked successful bootstrap replicates."
  } else if (id == "ST-03") {
    "Reader-facing labels follow publication_display_labels_v1.0.csv. Machine fields and Chinese labels are retained only in internal audit files."
  } else {
    "Abbreviations: AKI, acute kidney injury; ATE, average treatment effect; CBPS, covariate balancing propensity score; ESS, effective sample size; SMD, standardized mean difference."
  }
  if (id %in% c("MT-02", "MT-03")) {
    compact <- compact_outcome_table(out)
    make_docx(out, docx_path, title_map[[id]], note, subtables = list(compact), font_size = 7.5)
  } else if (id == "ST-06") {
    sub1 <- out[, c(1:2, 3:12), drop = FALSE]; attr(sub1, "subtitle") <- "Propensity-score overlap diagnostics"
    sub2 <- out[, c(1:2, 13:22), drop = FALSE]; attr(sub2, "subtitle") <- "Weight-distribution diagnostics"
    sub3 <- out[, c(1:2, 23:27), drop = FALSE]; attr(sub3, "subtitle") <- "Effective sample size and balance diagnostics"
    make_docx(out, docx_path, title_map[[id]], note, subtables = list(sub1, sub2, sub3), font_size = 7.5)
  } else {
    size <- if (id %in% c("ST-02", "ST-03", "ST-07")) 6.5 else 8
    make_docx(out, docx_path, title_map[[id]], note, font_size = size)
  }
  paths <- c(csv_path, docx_path)
  manifest_rows[[length(manifest_rows) + 1L]] <- data.frame(
    artifact_id = id, format = c("csv", "docx"), absolute_path = normalizePath(paths, winslash = "/", mustWork = TRUE),
    size_bytes = file.info(paths)$size,
    sha256 = vapply(paths, digest::digest, character(1), algo = "sha256", file = TRUE, serialize = FALSE),
    stringsAsFactors = FALSE
  )
  qa_rows[[length(qa_rows) + 1L]] <- data.frame(
    artifact_id = id, source_rows = nrow(sources[[id]]), publication_rows = nrow(out),
    publication_columns = ncol(out), reader_label_check = TRUE,
    csv_exists = file.exists(csv_path) && file.info(csv_path)$size > 0,
    docx_exists = file.exists(docx_path) && file.info(docx_path)$size > 0,
    stringsAsFactors = FALSE
  )
}

manifest <- do.call(rbind, manifest_rows)
qa <- do.call(rbind, qa_rows)
if (!all(qa$reader_label_check & qa$csv_exists & qa$docx_exists)) stop("Publication table generation QA failed.")
write_csv(manifest, file.path(audit_dir, "publication_table_manifest.csv"))
write_csv(qa, file.path(audit_dir, "publication_table_generation_QA.csv"))

cat("PUBLICATION_TABLE_IDS=", length(expected_ids), "\n", sep = "")
cat("PUBLICATION_TABLE_FILES=", nrow(manifest), "\n", sep = "")
cat("PUBLICATION_TABLE_GENERATION=PASS\n")
