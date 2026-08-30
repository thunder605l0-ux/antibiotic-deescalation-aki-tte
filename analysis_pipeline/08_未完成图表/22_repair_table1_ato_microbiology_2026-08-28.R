# ==============================================================================
# 脚本用途：在不改变主分析权重或结局估计的前提下，为主表 1 补充 T0 微生物学状态。
# 方法来源：项目锁定的 bias-reduced logistic propensity score + overlap weights；
#   二分类/多分类水平的加权标准化差异按 cobalt 常用的 pooled-denominator 口径计算。
# 引用 / commit：RD-4.6-ANALYSIS-OUTPUT-PLAN-V1.0；锁定 weights_*_imp*.rds 对象。
# 访问日期：2026-08-28。
# 适配说明：仅对既有 ATO 权重进行描述性汇总，不重估倾向评分、不运行 bootstrap、
#   不改变任何结局估计；组内权重比例乘以原始组样本量，仅用于与现行 Table 1 的
#   “加权计数（百分比）”展示格式保持一致。
# 验证命令：Rscript --vanilla 22_repair_table1_ato_microbiology_2026-08-28.R
# ==============================================================================

options(stringsAsFactors = FALSE, warn = 1)

project_root <- normalizePath(Sys.getenv("TTE_PROJECT_ROOT", unset = "."), winslash = "/", mustWork = TRUE)
runtime_root <- normalizePath(Sys.getenv("TTE_RUNTIME_ROOT", unset = "."), winslash = "/", mustWork = TRUE)
table_path <- file.path(project_root, "data", "Table 1_Weighted Baseline Characteristics.csv")
audit_path <- file.path(project_root, "data", "Table 1_Microbiology ATO audit_2026-08-28.csv")
weight_dir <- file.path(runtime_root, "09_输出结果", "03_权重")
raw_dir <- file.path(runtime_root, "01_原始分析数据")

trial_spec <- data.frame(
  trial = c("anti_mrsa", "anti_psa"),
  raw_file = c("analysis_dataset_anti_mrsa.csv", "analysis_dataset_anti_psa.csv"),
  deescalated_n = c(83, 57),
  continued_n = c(185, 315),
  stringsAsFactors = FALSE
)

level_codes <- 0:3
level_labels <- c(
  "  All final culture results negative",
  "  Non-target organism identified",
  "  Unresolved at T0",
  "  No culture obtained"
)

weighted_level_smd <- function(level_value, exposure, weights, microbiology) {
  indicator <- as.numeric(microbiology == level_value)
  p1 <- stats::weighted.mean(indicator[exposure == 1], weights[exposure == 1])
  p0 <- stats::weighted.mean(indicator[exposure == 0], weights[exposure == 0])
  denominator <- sqrt((p1 * (1 - p1) + p0 * (1 - p0)) / 2)
  if (!is.finite(denominator) || denominator == 0) return(0)
  abs((p1 - p0) / denominator)
}

audit_rows <- list()

for (trial_index in seq_len(nrow(trial_spec))) {
  spec <- trial_spec[trial_index, ]
  raw <- utils::read.csv(
    file.path(raw_dir, spec$raw_file),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  required_raw <- c("icu_stay_id", "microbiology_state_group4")
  if (!all(required_raw %in% names(raw))) {
    stop("Required microbiology fields are absent from ", spec$raw_file)
  }
  if (anyDuplicated(raw$icu_stay_id)) stop("Duplicate ICU-stay identifiers in ", spec$raw_file)

  imputation_rows <- list()
  for (imputation in seq_len(20L)) {
    weight_path <- file.path(
      weight_dir,
      sprintf("weights_%s_imp%02d.rds", spec$trial, imputation)
    )
    if (!file.exists(weight_path)) stop("Missing locked weight object: ", weight_path)
    weight_object <- readRDS(weight_path)
    exposure <- as.numeric(as.character(weight_object$completed_data$deescalation))
    weights <- as.numeric(weight_object$ato_brlogit$weight)
    icu_stay_id <- as.character(weight_object$key_map$icu_stay_id)
    raw_match <- match(icu_stay_id, as.character(raw$icu_stay_id))
    if (anyNA(raw_match)) stop("Could not align raw microbiology fields for ", spec$trial)
    microbiology <- as.integer(raw$microbiology_state_group4[raw_match])
    if (anyNA(microbiology) || !all(microbiology %in% level_codes)) {
      stop("Unexpected microbiology-state code in ", spec$trial)
    }

    for (level_index in seq_along(level_codes)) {
      level_value <- level_codes[level_index]
      indicator <- as.numeric(microbiology == level_value)
      p_deesc <- stats::weighted.mean(indicator[exposure == 1], weights[exposure == 1])
      p_cont <- stats::weighted.mean(indicator[exposure == 0], weights[exposure == 0])
      imputation_rows[[length(imputation_rows) + 1L]] <- data.frame(
        trial = spec$trial,
        imputation = imputation,
        level_code = level_value,
        level_label = level_labels[level_index],
        deescalation_proportion = p_deesc,
        continuation_proportion = p_cont,
        absolute_smd = weighted_level_smd(level_value, exposure, weights, microbiology),
        stringsAsFactors = FALSE
      )
    }
  }

  imputation_data <- do.call(rbind, imputation_rows)
  for (level_index in seq_along(level_codes)) {
    selected <- imputation_data$level_code == level_codes[level_index]
    audit_rows[[length(audit_rows) + 1L]] <- data.frame(
      trial = spec$trial,
      level_code = level_codes[level_index],
      level_label = level_labels[level_index],
      mean_deescalation_proportion = mean(imputation_data$deescalation_proportion[selected]),
      mean_continuation_proportion = mean(imputation_data$continuation_proportion[selected]),
      maximum_absolute_smd_across_imputations = max(imputation_data$absolute_smd[selected]),
      deescalated_n = spec$deescalated_n,
      continued_n = spec$continued_n,
      stringsAsFactors = FALSE
    )
  }
}

audit <- do.call(rbind, audit_rows)
utils::write.csv(audit, audit_path, row.names = FALSE, fileEncoding = "UTF-8")

table_data <- utils::read.csv(table_path, check.names = FALSE, stringsAsFactors = FALSE)
if (sum(table_data$Characteristic == "Microbiology status at T0") > 1L) {
  stop("Table 1 contains duplicate microbiology headers.")
}
if (any(table_data$Characteristic == "Microbiology status at T0")) {
  remove_rows <- table_data$Characteristic == "Microbiology status at T0" |
    table_data$Characteristic %in% level_labels
  table_data <- table_data[!remove_rows, , drop = FALSE]
}

trial_summary <- function(trial) {
  subset_data <- audit[audit$trial == trial, , drop = FALSE]
  subset_data[match(level_codes, subset_data$level_code), , drop = FALSE]
}

mrsa <- trial_summary("anti_mrsa")
psa <- trial_summary("anti_psa")

format_weighted <- function(proportion, n) {
  sprintf("%.1f (%.1f)", proportion * n, proportion * 100)
}

micro_rows <- data.frame(
  Characteristic = c("Microbiology status at T0", level_labels),
  `Anti-MRSA: De-escalation` = c("", vapply(seq_along(level_codes), function(i) {
    format_weighted(mrsa$mean_deescalation_proportion[i], mrsa$deescalated_n[i])
  }, character(1))),
  `Anti-MRSA: Continuation` = c("", vapply(seq_along(level_codes), function(i) {
    format_weighted(mrsa$mean_continuation_proportion[i], mrsa$continued_n[i])
  }, character(1))),
  `Anti-MRSA: SMD` = c(sprintf("%.3f", max(mrsa$maximum_absolute_smd_across_imputations)), rep("", 4)),
  `Antipseudomonal: De-escalation` = c("", vapply(seq_along(level_codes), function(i) {
    format_weighted(psa$mean_deescalation_proportion[i], psa$deescalated_n[i])
  }, character(1))),
  `Antipseudomonal: Continuation` = c("", vapply(seq_along(level_codes), function(i) {
    format_weighted(psa$mean_continuation_proportion[i], psa$continued_n[i])
  }, character(1))),
  `Antipseudomonal: SMD` = c(sprintf("%.3f", max(psa$maximum_absolute_smd_across_imputations)), rep("", 4)),
  check.names = FALSE,
  stringsAsFactors = FALSE
)

insert_after <- match("Source-control procedure", table_data$Characteristic)
if (is.na(insert_after)) stop("Could not locate the Table 1 insertion point.")
table_data <- rbind(
  table_data[seq_len(insert_after), , drop = FALSE],
  micro_rows,
  table_data[seq.int(insert_after + 1L, nrow(table_data)), , drop = FALSE]
)

utils::write.csv(
  table_data,
  table_path,
  row.names = FALSE,
  na = "",
  fileEncoding = "UTF-8"
)

cat("TABLE1_ROWS=", nrow(table_data), "\n", sep = "")
cat("MICROBIOLOGY_AUDIT_ROWS=", nrow(audit), "\n", sep = "")
cat("ANTI_MRSA_MAX_MICROBIOLOGY_SMD=", sprintf("%.3f", max(mrsa$maximum_absolute_smd_across_imputations)), "\n", sep = "")
cat("ANTIPSEUDOMONAL_MAX_MICROBIOLOGY_SMD=", sprintf("%.3f", max(psa$maximum_absolute_smd_across_imputations)), "\n", sep = "")
