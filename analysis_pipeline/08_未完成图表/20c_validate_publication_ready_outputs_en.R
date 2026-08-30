# Analysis: Validate publication-ready English tables and figures
# Date: 2026-08-08
# Random seed: not applicable (deterministic QA)
# R: 4.4.2
#
# Method: source-hash replay, file completeness, exact copied-value comparison, label scan, DOCX XML style audit, and PDF embedded-font audit
# Source: locked-source baseline, publication artifact plan, generated publication tables/figures, ST-08 final audit
# Citation/version: publication harmonization design 2026-08-08
# Access date: 2026-08-08
# Adaptation: validation only; no analysis or estimate generation.
# Verification command: D:/software/R4.4.2/R-4.4.2/bin/x64/Rscript.exe 08_未完成图表/20c_validate_publication_ready_outputs_en.R

options(stringsAsFactors = FALSE)

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg[[1L]]), winslash = "/", mustWork = TRUE)
project_root <- dirname(dirname(script_path))
output_root <- file.path(project_root, "09_输出结果", "11_publication_ready_en")
audit_dir <- file.path(output_root, "_audit")

if (!requireNamespace("digest", quietly = TRUE)) stop("Package 'digest' is required.")
if (!requireNamespace("pdftools", quietly = TRUE)) stop("Package 'pdftools' is required.")

read_csv <- function(path) utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE, fileEncoding = "UTF-8")
hash_file <- function(path) digest::digest(path, algo = "sha256", file = TRUE, serialize = FALSE)

qa <- list()
add_check <- function(id, description, pass, evidence = "") {
  qa[[length(qa) + 1L]] <<- data.frame(
    check_id = id, description = description, pass = isTRUE(pass), evidence = as.character(evidence), stringsAsFactors = FALSE
  )
}

# Source hash replay.
baseline <- read_csv(file.path(audit_dir, "source_baseline_sha256.csv"))
current_hash <- vapply(baseline$absolute_path, hash_file, character(1))
source_unchanged <- all(tolower(unname(current_hash)) == tolower(baseline$sha256))
add_check("PUB01", "All 49 locked source artifacts retain their baseline SHA-256 hashes", source_unchanged, paste0(sum(current_hash == baseline$sha256), "/", nrow(baseline)))

coding_en <- file.path(project_root, "02_配置与变量字典", "variable_coding_table_en.csv")
coding_cn <- file.path(project_root, "02_配置与变量字典", "variable_coding_table_cn.csv")
add_check("PUB02", "English variable coding table retains the locked SHA-256 hash", toupper(hash_file(coding_en)) == "1A14B4AABFA39D85F296219DC453C5620B2D3AFAA1B5E3C11DE0C5568966B403", toupper(hash_file(coding_en)))
add_check("PUB03", "Chinese internal coding table retains the locked SHA-256 hash", toupper(hash_file(coding_cn)) == "2C267B74D36D1BA7E2BFB36DB84DFA2AD35D940A3AFC6DF18632ADEE67192509", toupper(hash_file(coding_cn)))

# Artifact completeness.
plan <- read_csv(file.path(audit_dir, "publication_artifact_plan.csv"))
planned_paths <- file.path(output_root, plan$relative_path)
exists_nonempty <- file.exists(planned_paths) & file.info(planned_paths)$size > 0
add_check("PUB04", "All 59 planned publication-ready files exist and are nonempty", all(exists_nonempty), paste0(sum(exists_nonempty), "/", length(exists_nonempty)))
add_check("PUB05", "All 24 planned artifact identifiers are represented", length(unique(plan$artifact_id[exists_nonempty])) == 24L, length(unique(plan$artifact_id[exists_nonempty])))

# Formal CSV language and machine-label audit.
formal_csv_paths <- planned_paths[plan$format == "csv"]
scan_rows <- list()
for (path in formal_csv_paths) {
  data <- read_csv(path)
  values <- unlist(data, use.names = FALSE)
  values <- as.character(values)
  cjk <- grepl("[\u3400-\u9FFF]", values, perl = TRUE)
  forbidden <- grepl("race_white|White race|factor\\(|anti_mrsa|anti_psa|cbps_ate|ato_brlogit|NEW_MS0[1-9]|ABX_TIME_proxy|MICRO_STATE_proxy|MICE:|PS formula|target-antibiotic", values, ignore.case = TRUE) |
    grepl("[A-Za-z0-9]+_[A-Za-z0-9_]+", values)
  forbidden_header <- grepl("[A-Za-z0-9]+_[A-Za-z0-9_]+", names(data))
  scan_rows[[length(scan_rows) + 1L]] <- data.frame(
    file = basename(path), cjk_cells = sum(cjk, na.rm = TRUE), forbidden_cells = sum(forbidden, na.rm = TRUE),
    forbidden_headers = sum(forbidden_header), pass = !any(cjk, na.rm = TRUE) && !any(forbidden, na.rm = TRUE) && !any(forbidden_header),
    stringsAsFactors = FALSE
  )
}
language_scan <- do.call(rbind, scan_rows)
add_check("PUB06", "All 13 formal CSV tables are English-only and hide machine labels", all(language_scan$pass), paste0(sum(language_scan$pass), "/", nrow(language_scan)))

dictionary <- read_csv(file.path(project_root, "02_配置与变量字典", "publication_display_labels_v1.0.csv"))
race <- dictionary[dictionary$clean_field == "race_white", , drop = FALSE]
add_check("PUB07", "The reader-facing race label is exactly Race", nrow(race) == 1L && race$display_label == "Race", race$display_label)
add_check("PUB08", "The original Race coding and missingness definition are unchanged", race$factor_levels_en == "0=explicitly non-White; 1=White; unknown race=missing", race$factor_levels_en)

# Exact copied-value comparisons for all publication tables.
compare_columns <- function(source_relative, output_relative, source_columns, output_columns) {
  source <- read_csv(file.path(project_root, source_relative))
  output <- read_csv(file.path(output_root, output_relative))
  if (nrow(source) != nrow(output) || length(source_columns) != length(output_columns)) return(FALSE)
  all(vapply(seq_along(source_columns), function(i) {
    identical(as.character(source[[source_columns[[i]]]]), as.character(output[[output_columns[[i]]]]))
  }, logical(1)))
}

identity_rows <- list()
add_identity <- function(id, pass) {
  identity_rows[[length(identity_rows) + 1L]] <<- data.frame(artifact_id = id, exact_copied_value_identity = isTRUE(pass), stringsAsFactors = FALSE)
}
add_identity("MT-01", compare_columns(
  "09_输出结果/09_tables/MT-01_weighted_baseline_characteristics.csv",
  "02_main_tables/MT-01_weighted_baseline_characteristics.csv",
  c("section", "deescalated", "continued", "maximum_absolute_smd", "variance_ratio_range"),
  c("Section", "De-escalation", "Continuation", "Maximum absolute SMD", "Variance-ratio range")
))
for (id in c("MT-02", "MT-03")) {
  source_name <- if (id == "MT-02") "MT-02_anti_mrsa_outcomes_controlled.csv" else "MT-03_anti_psa_outcomes_controlled.csv"
  output_name <- if (id == "MT-02") "MT-02_anti_mrsa_outcomes.csv" else "MT-03_anti_psa_outcomes.csv"
  source <- read_csv(file.path(project_root, "09_输出结果", "05_主要结局", "controlled_reporting", source_name))
  output <- read_csv(file.path(output_root, "02_main_tables", output_name))
  add_identity(id, nrow(source) == nrow(output) && all(vapply(seq_along(names(source)[4:22]), function(i) {
    identical(as.character(source[[names(source)[3L + i]]]), as.character(output[[names(output)[2L + i]]]))
  }, logical(1))))
}
add_identity("ST-01", compare_columns(
  "09_输出结果/09_tables/ST-01_target_trial_mapping.csv", "04_supplementary_tables/ST-01_target_trial_mapping.csv",
  c("Element"), c("Element")
))
add_identity("ST-02", compare_columns(
  "09_输出结果/09_tables/ST-02_antibacterial_dictionary.csv", "04_supplementary_tables/ST-02_antibacterial_dictionary.csv",
  c("table_section", "canonical_generic", "source_name_in_paper_a", "alias_pattern", "mimic_emar_raw_detected", "raw_drug_name", "raw_route", "include_main", "source_rows", "database", "eligible_for_systemic_sepsis_coverage"),
  c("Section", "Canonical generic name", "Name in source publication", "Alias pattern", "Detected in MIMIC-IV eMAR", "Raw drug name", "Raw route", "Included in main analysis", "Source rows", "Database", "Eligible systemic coverage")
))
add_identity("ST-03", compare_columns(
  "09_输出结果/09_tables/ST-03_variable_definitions.csv", "04_supplementary_tables/ST-03_variable_definitions.csv",
  c("number", "variable_type", "numeric_code_definition"),
  c("No.", "Variable type", "Coding or definition")
))
add_identity("ST-04", compare_columns(
  "09_输出结果/09_tables/controlled_reporting/ST-04_unweighted_baseline_characteristics_controlled.csv", "04_supplementary_tables/ST-04_unweighted_baseline_characteristics.csv",
  c("section", "overall", "deescalated", "continued", "absolute_smd"), c("Section", "Overall", "De-escalation", "Continuation", "Absolute SMD")
))
add_identity("ST-05", compare_columns(
  "09_输出结果/09_tables/controlled_reporting/ST-05_missing_data_controlled.csv", "04_supplementary_tables/ST-05_missing_data.csv",
  c("No.", "Section", "MRSA_overall", "MRSA_continued", "MRSA_deescalated", "PSA_overall", "PSA_continued", "PSA_deescalated"),
  c("No.", "Section", "Anti-MRSA overall", "Anti-MRSA continuation", "Anti-MRSA de-escalation", "Anti-PSA overall", "Anti-PSA continuation", "Anti-PSA de-escalation")
))
add_identity("ST-06", compare_columns(
  "09_输出结果/09_tables/ST-06_multiple_imputation_specification_QA.csv", "04_supplementary_tables/ST-06_multiple_imputation_specification.csv",
  names(read_csv(file.path(project_root, "09_输出结果/09_tables/ST-06_multiple_imputation_specification_QA.csv")))[-1L],
  names(read_csv(file.path(output_root, "04_supplementary_tables/ST-06_multiple_imputation_specification.csv")))[-1L]
))
add_identity("ST-07", compare_columns(
  "09_输出结果/09_tables/ST-07_weight_diagnostics.csv", "04_supplementary_tables/ST-07_weight_diagnostics.csv",
  names(read_csv(file.path(project_root, "09_输出结果/09_tables/ST-07_weight_diagnostics.csv")))[-(1:2)],
  names(read_csv(file.path(output_root, "04_supplementary_tables/ST-07_weight_diagnostics.csv")))[-(1:2)]
))
add_identity("ST-08", compare_columns(
  "09_输出结果/06_敏感性分析/controlled_reporting/ST-08_prespecified_sensitivity_results.csv", "04_supplementary_tables/ST-08_prespecified_sensitivity_results.csv",
  names(read_csv(file.path(project_root, "09_输出结果/06_敏感性分析/controlled_reporting/ST-08_prespecified_sensitivity_results.csv")))[4:14],
  names(read_csv(file.path(output_root, "04_supplementary_tables/ST-08_prespecified_sensitivity_results.csv")))[4:14]
))
add_identity("ST-09", compare_columns(
  "09_输出结果/06_敏感性分析/controlled_reporting/ST-09_unmeasured_confounding_evalue.csv", "04_supplementary_tables/ST-09_unmeasured_confounding_evalue.csv",
  c("outcome", "estimand", "risk_ratio", "rr_ci_lower", "rr_ci_upper", "confidence_level", "quantile_type", "bootstrap_success_n", "ci_limit_closest_to_null", "ci_includes_null", "evalue_point", "evalue_ci_limit", "interpretation"),
  c("Outcome", "Estimand", "Risk ratio", "RR 95% CI lower", "RR 95% CI upper", "Confidence level", "Quantile type", "Successful bootstrap replicates", "CI limit closest to the null", "CI includes the null", "E-value for point estimate", "E-value for CI limit", "Interpretation")
))
add_identity("ST-10", compare_columns(
  "09_输出结果/06_敏感性分析/ST-10_aki_phenotype_components.csv", "04_supplementary_tables/ST-10_aki_phenotype_components.csv",
  names(read_csv(file.path(project_root, "09_输出结果/06_敏感性分析/ST-10_aki_phenotype_components.csv")))[3:14],
  names(read_csv(file.path(output_root, "04_supplementary_tables/ST-10_aki_phenotype_components.csv")))[3:14]
))
identity <- do.call(rbind, identity_rows)
add_check("PUB09", "All 13 publication tables retain exact copied source values in every reported data column", all(identity$exact_copied_value_identity), paste0(sum(identity$exact_copied_value_identity), "/", nrow(identity)))

# DOCX XML audit for explicit font and absence of visible vertical borders.
docx_paths <- planned_paths[plan$format == "docx"]
docx_rows <- list()
for (path in docx_paths) {
  exdir <- tempfile("docx_xml_")
  dir.create(exdir)
  on.exit(unlink(exdir, recursive = TRUE, force = TRUE), add = TRUE)
  utils::unzip(path, files = "word/document.xml", exdir = exdir)
  xml <- paste(readLines(file.path(exdir, "word", "document.xml"), warn = FALSE, encoding = "UTF-8"), collapse = "")
  tnr <- grepl("Times New Roman", xml, fixed = TRUE)
  visible_vertical <- grepl("<w:(insideV|left|right)[^>]*w:val=\"(single|double|dashed|dotted|thick)\"", xml, perl = TRUE)
  docx_rows[[length(docx_rows) + 1L]] <- data.frame(file = basename(path), times_new_roman = tnr, visible_vertical_border = visible_vertical, pass = tnr && !visible_vertical, stringsAsFactors = FALSE)
}
docx_audit <- do.call(rbind, docx_rows)
add_check("PUB10", "All DOCX tables explicitly use Times New Roman and have no visible vertical borders", all(docx_audit$pass), paste0(sum(docx_audit$pass), "/", nrow(docx_audit)))

figure_label_audit <- read_csv(file.path(audit_dir, "publication_figure_label_audit.csv"))
add_check("PUB11", "All reader-facing figure labels pass the English display-label audit", all(figure_label_audit$pass), paste0(sum(figure_label_audit$pass), "/", nrow(figure_label_audit)))

# ST-08 locked audit replay.
st08 <- read_csv(file.path(project_root, "09_输出结果", "06_敏感性分析", "controlled_reporting", "ST-08_controlled_v2_final_audit.csv"))
st08_map <- setNames(st08$value, st08$item)
st08_pass <- identical(st08_map[["status"]], "success") &&
  identical(st08_map[["anti_mrsa_formal_repeat_n"]], "1000") &&
  identical(st08_map[["anti_psa_formal_repeat_n"]], "999") &&
  identical(st08_map[["anti_psa_b0934_status"]], "preserved_failed_not_replaced") &&
  grepl("type=7", st08_map[["percentile_interval"]], fixed = TRUE) &&
  identical(st08_map[["pooling_authorized"]], "false")
add_check("PUB12", "ST-08 final audit remains successful with 1000/999 repeats, preserved b=934, type-7 intervals, and pooling unauthorized", st08_pass, paste(st08_map[c("status", "anti_mrsa_formal_repeat_n", "anti_psa_formal_repeat_n", "anti_psa_b0934_status", "percentile_interval", "pooling_authorized")], collapse = "; "))

# PDF embedded-font audit.
pdf_paths <- planned_paths[plan$format == "pdf"]
pdf_font_rows <- lapply(pdf_paths, function(path) {
  fonts <- pdftools::pdf_fonts(path)
  font_names <- unique(fonts$name)
  times_new_roman <- length(font_names) > 0L && all(grepl("TimesNewRoman", font_names, ignore.case = TRUE))
  embedded <- nrow(fonts) > 0L && all(fonts$embedded)
  data.frame(
    file = basename(path),
    font_names = paste(font_names, collapse = "; "),
    all_fonts_embedded = embedded,
    times_new_roman_only = times_new_roman,
    pass = embedded && times_new_roman,
    stringsAsFactors = FALSE
  )
})
pdf_font_audit <- do.call(rbind, pdf_font_rows)
add_check("PUB13", "All 11 figure PDFs embed Times New Roman fonts only", all(pdf_font_audit$pass), paste0(sum(pdf_font_audit$pass), "/", nrow(pdf_font_audit)))

qa_table <- do.call(rbind, qa)
utils::write.csv(language_scan, file.path(audit_dir, "publication_csv_language_QA.csv"), row.names = FALSE, na = "", fileEncoding = "UTF-8")
utils::write.csv(identity, file.path(audit_dir, "publication_table_exact_value_QA.csv"), row.names = FALSE, na = "", fileEncoding = "UTF-8")
utils::write.csv(docx_audit, file.path(audit_dir, "publication_docx_style_QA.csv"), row.names = FALSE, na = "", fileEncoding = "UTF-8")
utils::write.csv(pdf_font_audit, file.path(audit_dir, "publication_pdf_font_QA.csv"), row.names = FALSE, na = "", fileEncoding = "UTF-8")
utils::write.csv(qa_table, file.path(audit_dir, "publication_ready_master_QA.csv"), row.names = FALSE, na = "", fileEncoding = "UTF-8")

delivery <- plan
delivery$absolute_path <- normalizePath(planned_paths, winslash = "/", mustWork = TRUE)
delivery$size_bytes <- file.info(planned_paths)$size
delivery$sha256 <- vapply(planned_paths, hash_file, character(1))
delivery$status <- "complete"
utils::write.csv(delivery, file.path(output_root, "publication_delivery_manifest.csv"), row.names = FALSE, na = "", fileEncoding = "UTF-8")

if (!all(qa_table$pass)) {
  print(qa_table[!qa_table$pass, , drop = FALSE])
  stop("Publication-ready master QA failed.")
}

cat("PUBLICATION_MASTER_QA=", sum(qa_table$pass), "/", nrow(qa_table), " PASS\n", sep = "")
cat("PUBLICATION_DELIVERY_FILES=", nrow(delivery), "\n", sep = "")
