# ============================================================================== 
# Analysis: Rebuild Supplementary Table 8 with ATO E-values for all primary outcomes
# Date: 2026-08-27
# Random seed: not applicable (deterministic calculation from reported risk ratios)
# R: 4.4.2
# Key packages: EValue 4.1.4; flextable 0.9.11; officer 0.7.3
# Method: EValue::evalues.RR with an independently checked closed-form calculation.
# Purpose: Report E-values for the two prespecified primary outcomes (7-day AKI
#          and 30-day all-cause mortality) in each trial, using the same
#          overlap-weighting (ATO) risk ratios and percentile-bootstrap 95% CIs
#          displayed in Tables 2 and 3.
# Source type: GitHub + project reported aggregate results
# Source name or URL: https://github.com/mayamathur/evalue_package
# Citation or commit/version: EValue commit 333d022041951903744e097e052430700eabbbca;
#                             Tables 2 and 3, ATO outputs regenerated 2026-08-27.
# Access date: 2026-08-27
# Adaptation notes: Replaces the prior ST-09 input route, which extracted CBPS-ATE
#                   rows while the published supplementary table was labelled ATO.
#                   This script reads only aggregate publication tables, performs no
#                   bootstrap, does not pool trials, and does not read patient-level data.
# Verification command:
#   D:/software/R4.4.2/R-4.4.2/bin/x64/Rscript.exe --vanilla
#   12b_rebuild_ato_evalues_all_primary_outcomes_2026-08-27.R
# ============================================================================== 

options(stringsAsFactors = FALSE)

if (!requireNamespace("EValue", quietly = TRUE) ||
    !requireNamespace("flextable", quietly = TRUE) ||
    !requireNamespace("officer", quietly = TRUE)) {
  stop("Packages 'EValue', 'flextable', and 'officer' are required.")
}

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("Run this script with Rscript.")
script_path <- normalizePath(sub("^--file=", "", script_arg), winslash = "/", mustWork = TRUE)
project_root <- normalizePath(file.path(dirname(script_path), "..", "..", ".."), winslash = "/", mustWork = TRUE)
data_dir <- file.path(project_root, "data")

input_registry <- data.frame(
  trial = c("anti_mrsa", "anti_psa"),
  trial_label = c("Anti-MRSA trial", "Anti-pseudomonal trial"),
  input_file = c("Table 2_Anti-MRSA Outcomes.csv", "Table 3_Anti-PSA Outcomes.csv"),
  bootstrap_success_n = c(1000L, 999L),
  stringsAsFactors = FALSE
)
required_outcomes <- c("7-day AKI", "30-day all-cause mortality")
required_columns <- c(
  "Outcome", "Model", "Risk ratio", "RR 95% CI lower", "RR 95% CI upper",
  "Successful bootstrap replicates", "Confidence level", "Quantile type"
)

atomic_write_csv <- function(data, path) {
  temporary_path <- paste0(path, ".tmp")
  utils::write.csv(data, temporary_path, row.names = FALSE, na = "", fileEncoding = "UTF-8")
  if (file.exists(path)) unlink(path)
  if (!file.rename(temporary_path, path)) stop("Unable to atomically write: ", path)
}

atomic_write_docx <- function(document, path) {
  temporary_path <- paste0(path, ".tmp.docx")
  if (file.exists(temporary_path)) unlink(temporary_path)
  print(document, target = temporary_path)
  if (file.exists(path)) unlink(path)
  if (!file.rename(temporary_path, path)) stop("Unable to atomically write: ", path)
}

evalue_closed_form <- function(rr) {
  if (!is.finite(rr) || rr <= 0) stop("Risk ratio must be finite and positive.")
  strength <- max(rr, 1 / rr)
  strength + sqrt(strength * (strength - 1))
}

extract_trial <- function(registry_row) {
  input_path <- file.path(data_dir, registry_row$input_file)
  if (!file.exists(input_path)) stop("Missing ATO publication table: ", input_path)
  input <- utils::read.csv(input_path, check.names = FALSE, stringsAsFactors = FALSE)
  missing_columns <- setdiff(required_columns, names(input))
  if (length(missing_columns) > 0L) stop("Missing required columns in ", input_path, ": ", paste(missing_columns, collapse = ", "))
  selected <- input[
    input$Outcome %in% required_outcomes & input$Model == "Overlap weighting (ATO)",
    required_columns,
    drop = FALSE
  ]
  selected <- selected[match(required_outcomes, selected$Outcome), , drop = FALSE]
  if (nrow(selected) != length(required_outcomes) || anyNA(selected$Outcome)) {
    stop("Each ATO publication table must contain exactly one row for each primary outcome: ", input_path)
  }
  if (any(selected$`Successful bootstrap replicates` != registry_row$bootstrap_success_n) ||
      any(selected$`Confidence level` != 0.95) ||
      any(selected$`Quantile type` != 7)) {
    stop("Bootstrap or CI metadata do not match the registered ATO publication table: ", input_path)
  }
  data.frame(
    trial = registry_row$trial,
    trial_label = registry_row$trial_label,
    outcome = selected$Outcome,
    model = selected$Model,
    estimand = "ATO",
    risk_ratio = as.numeric(selected$`Risk ratio`),
    rr_ci_lower = as.numeric(selected$`RR 95% CI lower`),
    rr_ci_upper = as.numeric(selected$`RR 95% CI upper`),
    bootstrap_success_n = as.integer(selected$`Successful bootstrap replicates`),
    confidence_level = as.numeric(selected$`Confidence level`),
    quantile_type = as.integer(selected$`Quantile type`),
    input_file = registry_row$input_file,
    input_md5 = toupper(unname(tools::md5sum(input_path))),
    stringsAsFactors = FALSE
  )
}

evalue_table <- do.call(rbind, lapply(seq_len(nrow(input_registry)), function(i) extract_trial(input_registry[i, , drop = FALSE])))
if (nrow(evalue_table) != 4L || any(!is.finite(as.matrix(evalue_table[, c("risk_ratio", "rr_ci_lower", "rr_ci_upper")]))) ||
    any(evalue_table$rr_ci_lower > evalue_table$risk_ratio) || any(evalue_table$risk_ratio > evalue_table$rr_ci_upper)) {
  stop("The four ATO primary risk-ratio inputs failed validation.")
}

evalue_table$ci_includes_null <- evalue_table$rr_ci_lower <= 1 & evalue_table$rr_ci_upper >= 1
evalue_table$ci_limit_closest_to_null <- ifelse(
  evalue_table$risk_ratio < 1,
  evalue_table$rr_ci_upper,
  evalue_table$rr_ci_lower
)
package_results <- lapply(seq_len(nrow(evalue_table)), function(i) {
  EValue::evalues.RR(
    est = evalue_table$risk_ratio[i],
    lo = evalue_table$rr_ci_lower[i],
    hi = evalue_table$rr_ci_upper[i],
    true = 1
  )
})
evalue_table$evalue_point <- vapply(package_results, function(x) unname(x["E-values", "point"]), numeric(1))
evalue_table$evalue_ci_limit <- vapply(seq_len(nrow(evalue_table)), function(i) {
  if (evalue_table$ci_includes_null[i]) return(1)
  if (evalue_table$risk_ratio[i] < 1) {
    unname(package_results[[i]]["E-values", "upper"])
  } else {
    unname(package_results[[i]]["E-values", "lower"])
  }
}, numeric(1))
closed_point <- vapply(evalue_table$risk_ratio, evalue_closed_form, numeric(1))
closed_ci <- ifelse(evalue_table$ci_includes_null, 1, vapply(evalue_table$ci_limit_closest_to_null, evalue_closed_form, numeric(1)))
if (any(abs(evalue_table$evalue_point - closed_point) > 1e-10) ||
    any(abs(evalue_table$evalue_ci_limit - closed_ci) > 1e-10)) {
  stop("EValue package output did not match the independently checked closed-form calculation.")
}
evalue_table$interpretation <- paste(
  "The E-value is the minimum association strength an unmeasured confounder would need with both treatment strategy and outcome, conditional on measured covariates, to move the estimate to the null; it does not establish absence of residual confounding."
)

out_csv <- file.path(data_dir, "Supplementary Table 8_E-values for Unmeasured Confounding.csv")
out_docx <- file.path(data_dir, "Supplementary Table 8_E-values for Unmeasured Confounding.docx")
audit_csv <- file.path(data_dir, "Supplementary Table 8_E-values for Unmeasured Confounding_audit_2026-08-27.csv")
atomic_write_csv(evalue_table, out_csv)

display_table <- data.frame(
  Trial = evalue_table$trial_label,
  Outcome = evalue_table$outcome,
  `RR (95% CI)` = sprintf("%.2f (%.2f to %.2f)", evalue_table$risk_ratio, evalue_table$rr_ci_lower, evalue_table$rr_ci_upper),
  `E-value, point estimate` = sprintf("%.2f", evalue_table$evalue_point),
  `E-value, CI limit closest to null` = sprintf("%.2f", evalue_table$evalue_ci_limit),
  `Successful bootstrap replicates` = evalue_table$bootstrap_success_n,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
ft <- flextable::flextable(display_table)
ft <- flextable::theme_booktabs(ft)
ft <- flextable::font(ft, fontname = "Arial", part = "all")
ft <- flextable::fontsize(ft, size = 9, part = "all")
ft <- flextable::align(ft, j = 1:2, align = "left", part = "all")
ft <- flextable::align(ft, j = 3:6, align = "center", part = "all")
ft <- flextable::valign(ft, valign = "center", part = "all")
ft <- flextable::padding(ft, padding.top = 3, padding.bottom = 3, padding.left = 4, padding.right = 4, part = "all")
ft <- flextable::width(ft, j = 1, width = 1.25)
ft <- flextable::width(ft, j = 2, width = 1.55)
ft <- flextable::width(ft, j = 3, width = 1.55)
ft <- flextable::width(ft, j = 4, width = 1.1)
ft <- flextable::width(ft, j = 5, width = 1.55)
ft <- flextable::width(ft, j = 6, width = 1.0)
ft <- flextable::set_table_properties(ft, layout = "fixed", width = 1)

doc <- officer::read_docx()
doc <- officer::body_add_par(doc, "Supplementary Table 8. E-values for Unmeasured Confounding of the Primary ATO Outcome Risk Ratios", style = "Normal")
doc <- flextable::body_add_flextable(doc, value = ft)
note <- paste(
  "Abbreviations: AKI, acute kidney injury; ATO, average treatment effect in the overlap population; CI, confidence interval; RR, risk ratio.",
  "All inputs are the reported overlap-weighting (ATO) primary outcome RRs and 95% percentile-bootstrap CIs from Tables 2 and 3.",
  "Because each 95% CI includes 1.00, every CI-limit E-value is 1.00.",
  "E-values quantify the minimum strength of association that an unmeasured confounder would need with both treatment strategy and outcome, conditional on measured covariates, to move an estimate to the null; they do not demonstrate absence of residual confounding."
)
doc <- officer::body_add_fpar(doc, officer::fpar(officer::ftext(paste0("Note. ", note), officer::fp_text(font.family = "Arial", font.size = 8))))
doc <- officer::body_set_default_section(
  doc,
  officer::prop_section(
    page_size = officer::page_size(orient = "landscape"),
    page_margins = officer::page_mar(top = 0.6, bottom = 0.6, left = 0.6, right = 0.6)
  )
)
atomic_write_docx(doc, out_docx)

audit_table <- evalue_table[, c(
  "trial", "outcome", "model", "estimand", "input_file", "input_md5",
  "risk_ratio", "rr_ci_lower", "rr_ci_upper", "bootstrap_success_n",
  "confidence_level", "quantile_type", "ci_includes_null",
  "evalue_point", "evalue_ci_limit"
)]
audit_table$method_source <- "mayamathur/evalue_package commit 333d022041951903744e097e052430700eabbbca"
audit_table$input_scope <- "aggregate publication Table 2/Table 3 ATO primary outcomes only"
audit_table$pooling_authorized <- FALSE
atomic_write_csv(audit_table, audit_csv)

cat("EVALUE_INPUT_ESTIMAND=ATO\n")
cat("EVALUE_PRIMARY_OUTCOMES=4\n")
cat("EVALUE_BOOTSTRAP_RERUN=NO\n")
print(display_table, row.names = FALSE)
