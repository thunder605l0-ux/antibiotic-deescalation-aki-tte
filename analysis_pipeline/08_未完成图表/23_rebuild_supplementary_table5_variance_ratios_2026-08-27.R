# ============================================================================== 
# Analysis: Rebuild Supplementary Table 5 with ATO continuous-covariate variance ratios
# Date: 2026-08-27
# Random seed: not applicable (deterministic formatting of locked aggregate outputs)
# R: 4.4.2
# Key packages: flextable 0.9.11; officer 0.7.3
# Method: Reformat locked aggregate diagnostics; no propensity-score, imputation,
#         weighting, bootstrap, or outcome model is refit.
# Purpose: Make the continuous-covariate variance-ratio diagnostic cited in the
#          Limitations directly verifiable and explicitly identify it as primary ATO.
# Source type: project locked aggregate outputs + self-written formatting code
# Source name or URL: MT-01_weighted_baseline_characteristics.csv; Supplementary Table 5
# Citation or commit/version: D-4.7N-ANALYSIS-LOCK-V2.0; D-4.6-ANALYSIS-OUTPUT-PLAN-V1.0
# Access date: 2026-08-27
# Adaptation notes: Reads only aggregate outputs. The variance-ratio range is across
#                   20 completed data sets after primary overlap weighting (ATO).
# Verification command:
#   D:/software/R4.4.2/R-4.4.2/bin/x64/Rscript.exe --vanilla
#   23_rebuild_supplementary_table5_variance_ratios_2026-08-27.R
# ============================================================================== 

options(stringsAsFactors = FALSE)
if (!requireNamespace("flextable", quietly = TRUE) || !requireNamespace("officer", quietly = TRUE)) {
  stop("Packages 'flextable' and 'officer' are required.")
}

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("Run this script with Rscript.")
script_path <- normalizePath(sub("^--file=", "", script_arg), winslash = "/", mustWork = TRUE)
project_root <- normalizePath(Sys.getenv("TTE_PROJECT_ROOT", unset = "."), winslash = "/", mustWork = TRUE)
runtime_root <- normalizePath(Sys.getenv("TTE_RUNTIME_ROOT", unset = "."), winslash = "/", mustWork = TRUE)
data_dir <- file.path(project_root, "data")
mt01_path <- file.path(runtime_root, "09_输出结果", "09_tables", "MT-01_weighted_baseline_characteristics.csv")
existing_st05_path <- file.path(data_dir, "Supplementary Table 5_Weighting Diagnostics.csv")
table1_path <- file.path(data_dir, "Table 1_Weighted Baseline Characteristics.csv")

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

if (!file.exists(mt01_path) || !file.exists(existing_st05_path) || !file.exists(table1_path)) stop("Required aggregate diagnostic input is missing.")
mt01 <- utils::read.csv(mt01_path, check.names = FALSE, stringsAsFactors = FALSE)
panel_a_source <- utils::read.csv(existing_st05_path, check.names = FALSE, stringsAsFactors = FALSE)
table1 <- utils::read.csv(table1_path, check.names = FALSE, stringsAsFactors = FALSE)
required_mt01 <- c("trial", "row_type", "characteristic", "variance_ratio_range")
if (length(setdiff(required_mt01, names(mt01))) > 0L) stop("MT-01 is missing required variance-ratio fields.")

# 仅从已发布的聚合表重取总体最大 SMD，保证与 Table 1 使用同一口径。
required_table1 <- c("Anti-MRSA: SMD", "Anti-pseudomonal: SMD")
if (length(setdiff(required_table1, names(table1))) > 0L) stop("Table 1 is missing SMD columns.")
max_smd <- vapply(required_table1, function(column) {
  values <- suppressWarnings(as.numeric(table1[[column]]))
  max(abs(values[is.finite(values)]))
}, numeric(1))
if (any(!is.finite(max_smd))) stop("Unable to recover maximum SMD from Table 1.")

panel_a_out <- panel_a_source[panel_a_source$Panel == "A. Overall weighting diagnostics", c("Panel", "Trial", "Characteristic", "Value", "Detail"), drop = FALSE]
if (nrow(panel_a_out) != 4L) stop("Supplementary Table 5 panel A is incomplete.")
for (trial in c("Anti-MRSA trial", "Anti-PSA trial")) {
  row <- panel_a_out$Trial == trial & grepl("Overlap weighting", panel_a_out$Characteristic, fixed = TRUE)
  if (sum(row) != 1L) stop("Primary ATO row is not uniquely identified for ", trial)
  smd_name <- if (trial == "Anti-MRSA trial") "Anti-MRSA: SMD" else "Anti-pseudomonal: SMD"
  panel_a_out$Value[row] <- sprintf("%.3f", max_smd[[smd_name]])
  panel_a_out$Characteristic[row] <- "Overlap weighting (ATO; primary analysis; maximum |SMD| across 20 completed data sets)"
}

panel_b <- mt01[mt01$row_type == "continuous", c("trial", "characteristic", "variance_ratio_range"), drop = FALSE]
trial_map <- c(anti_mrsa = "Anti-MRSA trial", anti_psa = "Anti-pseudomonal trial")
panel_b$Trial <- unname(trial_map[panel_b$trial])
if (anyNA(panel_b$Trial) || any(!nzchar(panel_b$variance_ratio_range))) stop("Continuous ATO variance-ratio input is incomplete.")
split_range <- do.call(rbind, strsplit(panel_b$variance_ratio_range, "-", fixed = TRUE))
panel_b$range_lower <- as.numeric(split_range[, 1])
panel_b$range_upper <- as.numeric(split_range[, 2])
if (any(!is.finite(panel_b$range_lower)) || any(!is.finite(panel_b$range_upper))) stop("Variance-ratio range parsing failed.")
panel_b$`Within 0.5–2.0 in all imputations` <- ifelse(panel_b$range_lower >= 0.5 & panel_b$range_upper <= 2, "Yes", "No")
panel_b <- panel_b[, c("Trial", "characteristic", "variance_ratio_range", "Within 0.5–2.0 in all imputations")]
names(panel_b)[2:3] <- c("Continuous covariate", "Variance-ratio range across 20 imputations")

out_csv <- file.path(data_dir, "Supplementary Table 5_Weighting Diagnostics.csv")
out_docx <- file.path(data_dir, "Supplementary Table 5_Weighting Diagnostics.docx")
audit_csv <- file.path(data_dir, "Supplementary Table 5_Weighting Diagnostics_audit_2026-08-27.csv")
panel_b_out <- panel_b
panel_b_out$Panel <- "B. Continuous-covariate variance ratios after primary ATO"
combined_csv <- rbind(
  panel_a_out,
  data.frame(Panel = panel_b_out$Panel, Trial = panel_b_out$Trial, Characteristic = panel_b_out$`Continuous covariate`, Value = panel_b_out$`Variance-ratio range across 20 imputations`, Detail = panel_b_out$`Within 0.5–2.0 in all imputations`, stringsAsFactors = FALSE)
)
atomic_write_csv(combined_csv, out_csv)

ft_a <- flextable::flextable(panel_a_out)
ft_a <- flextable::theme_booktabs(ft_a)
ft_b <- flextable::flextable(panel_b)
ft_b <- flextable::theme_booktabs(ft_b)
for (name in c("ft_a", "ft_b")) {
  current <- get(name)
  current <- flextable::font(current, fontname = "Arial", part = "all")
  current <- flextable::fontsize(current, size = 8, part = "all")
  current <- flextable::valign(current, valign = "center", part = "all")
  current <- flextable::padding(current, padding.top = 2, padding.bottom = 2, padding.left = 3, padding.right = 3, part = "all")
  assign(name, current)
}
ft_a <- flextable::autofit(ft_a)
ft_b <- flextable::autofit(ft_b)
doc <- officer::read_docx()
doc <- officer::body_set_default_section(doc, officer::prop_section(page_size = officer::page_size(orient = "landscape"), page_margins = officer::page_mar(top = 0.55, bottom = 0.55, left = 0.55, right = 0.55)))
doc <- officer::body_add_par(doc, "Supplementary Table 5. Weighting Diagnostics for the Primary ATO Analysis", style = "Normal")
doc <- officer::body_add_par(doc, "A. Overall weighting diagnostics", style = "Normal")
doc <- flextable::body_add_flextable(doc, ft_a)
doc <- officer::body_add_par(doc, "B. Continuous-covariate variance ratios after primary overlap weighting", style = "Normal")
doc <- flextable::body_add_flextable(doc, ft_b)
note <- paste(
  "Abbreviations: ATO, average treatment effect in the overlap population; ESS, effective sample size; PS, propensity score; SMD, standardised mean difference.",
  "Panel A summarises the primary ATO and complementary CBPS-ATE diagnostics. Panel B reports the range of the weighted variance ratio (de-escalation divided by continuation) across 20 imputed data sets for every continuous covariate in the primary ATO analysis.",
  "The prespecified variance-ratio target was 0.5 to 2.0; a range outside that interval indicates residual distribution-shape imbalance despite a small SMD."
)
doc <- officer::body_add_fpar(doc, officer::fpar(officer::ftext(paste0("Note. ", note), officer::fp_text(font.family = "Arial", font.size = 8))))
atomic_write_docx(doc, out_docx)

audit <- panel_b
audit$source_file <- basename(mt01_path)
audit$source_md5 <- toupper(unname(tools::md5sum(mt01_path)))
audit$estimand <- "ATO"
audit$imputations <- 20L
atomic_write_csv(audit, audit_csv)
cat("ST05_ESTIMAND=ATO\nST05_CONTINUOUS_COVARIATES=", nrow(panel_b), "\n", sep = "")
