# ==============================================================================
# Analysis: MT-02/MT-03 受控主要结局表
# Date: 2026-08-01
# Random seed: 不使用随机数
# R: 4.4.2
# Analysis Lock: D-4.7N-ANALYSIS-LOCK-V2.0
# Output Plan: D-4.6-ANALYSIS-OUTPUT-PLAN-V1.0
#
# Method purpose:
#   使用原样本20份插补的锁定权重对象计算点估计，并连接06b产生的受控
#   Bootstrap双侧95%百分位区间。anti-MRSA使用1000个成功重复；anti-PSA
#   使用999个成功重复，b=934保持failed且不生成成功checkpoint。
#
# Method sources:
#   1) R survival::survfit documentation (Aalen-Johansen estimator).
#      https://stat.ethz.ch/R-manual/R-devel/library/survival/html/survfit.formula.html
#   2) R stats::quantile documentation (type=7 percentile endpoints).
#      https://stat.ethz.ch/R-manual/R-devel/library/stats/html/quantile.html
#   3) Hyndman RJ, Fan Y. Sample Quantiles in Statistical Packages.
#      The American Statistician. 1996;50(4):361-365.
#      https://doi.org/10.1080/00031305.1996.10473566
#   4) Pawel S, Bartoš F, Siepe BS, Lohmann A. Handling Missingness,
#      Failures, and Non-Convergence in Simulation Studies.
#      The American Statistician. 2026;80(1):31-48.
#      https://doi.org/10.1080/00031305.2025.2540002
#   5) flextable project, commit 5f357b490ef202e78480e60994922aa864b06a93.
#      https://github.com/davidgohel/flextable
#
# Citation/version and access date:
#   survival 3.7-0; R 4.4.2 stats::quantile(type=7); sources accessed
#   2026-08-01. Project methods are additionally governed by the two locks above.
#
# Adaptation notes:
#   点估计仍按原07脚本规则：每份插补先估计策略风险，再对20份风险取均值
#   并派生RD/RR。置信区间只读取06b受控池化结果，不重抽样、不重新插补。
#   风险差方向固定为de-escalation minus continuation。数值追溯CSV保留完整
#   精度；DOCX仅作显示格式化。旧07脚本和旧输出均不覆盖。
#
# Verification commands:
#   Rscript --vanilla -e \
#     "invisible(parse(file='06_主要结局分析/07b_primary_outcome_tables_controlled.R'))"
#   Rscript --vanilla 06_主要结局分析/07b_primary_outcome_tables_controlled.R
#
# RStudio execution:
#   source("06_主要结局分析/07b_primary_outcome_tables_controlled.R",
#          encoding="UTF-8")
# ==============================================================================

options(stringsAsFactors = FALSE)

# 同时兼容Rscript和RStudio source()；不支持逐行执行整个文件。
get_current_script_path <- function() {
  full_arguments <- commandArgs(trailingOnly = FALSE)
  file_arguments <- full_arguments[startsWith(full_arguments, "--file=")]
  if (length(file_arguments) == 1L) {
    return(normalizePath(
      sub("--file=", "", file_arguments, fixed = TRUE),
      winslash = "/",
      mustWork = TRUE
    ))
  }

  frame_files <- unlist(lapply(
    sys.frames(),
    function(current_frame) {
      if (is.null(current_frame$ofile)) {
        return(character(0))
      }
      as.character(current_frame$ofile)
    }
  ), use.names = FALSE)
  if (length(frame_files) > 0L) {
    return(normalizePath(
      frame_files[[length(frame_files)]],
      winslash = "/",
      mustWork = TRUE
    ))
  }
  stop("无法定位当前脚本；请使用Rscript或source()运行本文件。")
}

script_path <- get_current_script_path()
package_root_from_script <- dirname(dirname(script_path))
engine_path <- file.path(
  package_root_from_script,
  "05_Bootstrap",
  "04_bootstrap_single_repeat.R"
)
if (!file.exists(engine_path)) {
  stop("缺少锁定结局估计引擎：", engine_path)
}

# 加载统一配置、原子写出函数、权重和结局估计函数。
source(engine_path, local = globalenv(), encoding = "UTF-8")
setwd(package_root)
require_packages(c("survival", "flextable", "officer"))

expected_lock <- "D-4.7N-ANALYSIS-LOCK-V2.0"
expected_plan <- "D-4.6-ANALYSIS-OUTPUT-PLAN-V1.0"
expected_confidence_level <- 0.95
expected_quantile_type <- 7L
expected_models <- c("cbps_ate", "ato")
effect_fields <- c(
  "risk_aki_0",
  "risk_aki_1",
  "rd_aki",
  "rr_aki",
  "risk_death30_0",
  "risk_death30_1",
  "rd_death30",
  "rr_death30"
)
expected_success_n <- c(anti_mrsa = 1000L, anti_psa = 999L)

controlled_pooling_dir <- file.path(
  output_root,
  "04_Bootstrap",
  "controlled_pooling"
)
ci_path <- file.path(
  controlled_pooling_dir,
  "bootstrap_percentile_ci_controlled.csv"
)
pooling_audit_path <- file.path(
  controlled_pooling_dir,
  "bootstrap_pooling_audit.csv"
)
failure_path <- file.path(
  output_root,
  "04_Bootstrap",
  "failures",
  "anti_psa_b0934_failure.rds"
)
failed_checkpoint_path <- file.path(
  output_root,
  "04_Bootstrap",
  "checkpoints",
  "anti_psa_b0934_checkpoint.rds"
)

required_inputs <- c(ci_path, pooling_audit_path, failure_path)
if (!all(file.exists(required_inputs))) {
  stop(
    "受控主要结局表输入不完整：",
    paste(required_inputs[!file.exists(required_inputs)], collapse = "; ")
  )
}
if (file.exists(failed_checkpoint_path)) {
  stop("anti-PSA b=934不应存在成功checkpoint。")
}

# 保存失败证据摘要，确保本脚本执行期间不改变既有证据。
failure_md5_before <- unname(as.character(tools::md5sum(failure_path)))

ci_table <- utils::read.csv(
  ci_path,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
pooling_audit <- utils::read.csv(
  pooling_audit_path,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_ci_columns <- c(
  "trial",
  "model",
  "estimand_field",
  "n_success",
  "confidence_level",
  "quantile_type",
  "ci_lower",
  "ci_upper"
)
if (!all(required_ci_columns %in% names(ci_table))) {
  stop("受控Bootstrap区间表缺少必需字段。")
}
if (nrow(ci_table) != 32L) {
  stop("受控Bootstrap标量区间表应为32行。")
}

ci_keys <- do.call(
  paste,
  c(ci_table[c("trial", "model", "estimand_field")], sep = "|")
)
expected_ci_grid <- expand.grid(
  trial = trial_registry$trial,
  model = expected_models,
  estimand_field = effect_fields,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)
expected_ci_keys <- do.call(paste, c(expected_ci_grid, sep = "|"))
if (!(
  !anyDuplicated(ci_keys) &&
    setequal(ci_keys, expected_ci_keys) &&
    all(ci_table$n_success == expected_success_n[ci_table$trial]) &&
    all(ci_table$confidence_level == expected_confidence_level) &&
    all(as.integer(ci_table$quantile_type) == expected_quantile_type) &&
    all(is.finite(ci_table$ci_lower)) &&
    all(is.finite(ci_table$ci_upper)) &&
    all(ci_table$ci_lower <= ci_table$ci_upper)
)) {
  stop("受控Bootstrap区间的网格、成功数或端点校验失败。")
}

# 校验风险、风险差和风险比区间的自然取值范围。
risk_rows <- startsWith(ci_table$estimand_field, "risk_")
rd_rows <- startsWith(ci_table$estimand_field, "rd_")
rr_rows <- startsWith(ci_table$estimand_field, "rr_")
if (!(
  all(ci_table$ci_lower[risk_rows] >= 0 & ci_table$ci_upper[risk_rows] <= 1) &&
    all(ci_table$ci_lower[rd_rows] >= -1 & ci_table$ci_upper[rd_rows] <= 1) &&
    all(ci_table$ci_lower[rr_rows] >= 0)
)) {
  stop("受控Bootstrap区间超出估计量自然范围。")
}

# 从06b审计表取得唯一指标值。
get_audit_value <- function(trial, metric) {
  current <- pooling_audit[
    pooling_audit$trial == trial & pooling_audit$metric == metric,
    "value",
    drop = TRUE
  ]
  if (length(current) != 1L) {
    stop("受控池化审计指标缺失或重复：", trial, " / ", metric)
  }
  as.character(current)
}

if (!(
  get_audit_value("overall", "analysis_lock") == expected_lock &&
    get_audit_value("overall", "pooling_authorized") == "TRUE" &&
    as.integer(get_audit_value("anti_mrsa", "success_n")) == 1000L &&
    as.integer(get_audit_value("anti_psa", "success_n")) == 999L &&
    as.integer(get_audit_value("anti_psa", "preserved_failed_n")) == 1L &&
    get_audit_value("anti_psa", "preserved_failed_indices") == "934"
)) {
  stop("06b受控池化审计与批准方案不一致。")
}

imputation_rows <- list()
raw_count_rows <- list()
raw_count_reference <- list()

# 使用原样本20份锁定权重对象计算点估计，并逐份核验原始事件数不漂移。
for (trial in trial_registry$trial) {
  registry_row <- trial_registry[trial_registry$trial == trial, , drop = FALSE]
  for (imputation_index in seq_len(mice_m)) {
    object_path <- file.path(
      output_root,
      "03_权重",
      paste0(
        "weights_",
        trial,
        "_imp",
        sprintf("%02d", imputation_index),
        ".rds"
      )
    )
    if (!file.exists(object_path)) {
      stop("缺少锁定权重对象：", object_path)
    }
    saved <- readRDS(object_path)
    if (!(
      identical(saved$trial, trial) &&
        as.integer(saved$imputation) == imputation_index &&
        nrow(saved$completed_data) == registry_row$expected_n
    )) {
      stop("权重对象身份或样本量异常：", basename(object_path))
    }

    treatment <- as.integer(as.character(saved$completed_data$deescalation))
    if (!(
      !anyNA(treatment) &&
        sum(treatment == 1L) == registry_row$expected_deescalation_n &&
        sum(treatment == 0L) == registry_row$expected_continuation_n
    )) {
      stop("原样本策略组人数异常：", basename(object_path))
    }
    for (field in c("aki_7d", "death_30d")) {
      # 锁定对象中二元结局可能保存为因子；先经字符层安全转换为0/1整数。
      values <- suppressWarnings(as.integer(as.character(
        saved$completed_data[[field]]
      )))
      if (anyNA(values) || !all(values %in% c(0, 1))) {
        stop("原样本结局不是完整0/1字段：", basename(object_path), " / ", field)
      }
    }

    aki_values <- as.integer(as.character(saved$completed_data$aki_7d))
    death30_values <- as.integer(as.character(saved$completed_data$death_30d))

    current_counts <- data.frame(
      trial = trial,
      continuation_n = sum(treatment == 0L),
      deescalation_n = sum(treatment == 1L),
      continuation_aki_events = sum(aki_values[treatment == 0L]),
      deescalation_aki_events = sum(aki_values[treatment == 1L]),
      continuation_death30_events = sum(
        death30_values[treatment == 0L]
      ),
      deescalation_death30_events = sum(
        death30_values[treatment == 1L]
      ),
      stringsAsFactors = FALSE
    )
    if (imputation_index == 1L) {
      raw_count_reference[[trial]] <- current_counts
      raw_count_rows[[length(raw_count_rows) + 1L]] <- current_counts
    } else if (!identical(current_counts, raw_count_reference[[trial]])) {
      stop(trial, "的原始事件数在不同插补对象间发生漂移。")
    }

    for (model in expected_models) {
      weight <- if (model == "cbps_ate") {
        saved$cbps_ate$weight
      } else {
        saved$ato_brlogit$weight
      }
      estimates <- estimate_bootstrap_outcomes(saved$completed_data, weight)
      imputation_rows[[length(imputation_rows) + 1L]] <- data.frame(
        trial = trial,
        imputation = imputation_index,
        model = model,
        risk_aki_0 = estimates$risk_aki_0,
        risk_aki_1 = estimates$risk_aki_1,
        risk_death30_0 = estimates$risk_death30_0,
        risk_death30_1 = estimates$risk_death30_1,
        stringsAsFactors = FALSE
      )
    }
  }
}

imputation_table <- do.call(rbind, imputation_rows)
raw_count_table <- do.call(rbind, raw_count_rows)
point_rows <- list()

for (trial in trial_registry$trial) {
  for (model in expected_models) {
    current <- imputation_table[
      imputation_table$trial == trial & imputation_table$model == model,
      ,
      drop = FALSE
    ]
    if (nrow(current) != mice_m) {
      stop(trial, " / ", model, "不是20份插补点估计。")
    }
    pooled <- pool_bootstrap_imputations(current)
    point_rows[[length(point_rows) + 1L]] <- data.frame(
      trial = trial,
      model = model,
      risk_aki_0 = unname(pooled["risk_aki_0"]),
      risk_aki_1 = unname(pooled["risk_aki_1"]),
      rd_aki = unname(pooled["rd_aki"]),
      rr_aki = unname(pooled["rr_aki"]),
      risk_death30_0 = unname(pooled["risk_death30_0"]),
      risk_death30_1 = unname(pooled["risk_death30_1"]),
      rd_death30 = unname(pooled["rd_death30"]),
      rr_death30 = unname(pooled["rr_death30"]),
      stringsAsFactors = FALSE
    )
  }
}
point_table <- do.call(rbind, point_rows)

if (!(
  nrow(imputation_table) == 80L &&
    nrow(point_table) == 4L &&
    all(is.finite(as.matrix(point_table[effect_fields]))) &&
    all(as.matrix(point_table[c(
      "risk_aki_0",
      "risk_aki_1",
      "risk_death30_0",
      "risk_death30_1"
    )]) >= 0) &&
    all(as.matrix(point_table[c(
      "risk_aki_0",
      "risk_aki_1",
      "risk_death30_0",
      "risk_death30_1"
    )]) <= 1) &&
    all(point_table$rr_aki > 0) &&
    all(point_table$rr_death30 > 0)
)) {
  stop("原样本点估计结构或有限值校验失败。")
}

outcome_specification <- data.frame(
  outcome = c("7-day AKI", "30-day all-cause mortality"),
  outcome_code = c("aki", "death30"),
  continuation_event_field = c(
    "continuation_aki_events",
    "continuation_death30_events"
  ),
  deescalation_event_field = c(
    "deescalation_aki_events",
    "deescalation_death30_events"
  ),
  risk_0_field = c("risk_aki_0", "risk_death30_0"),
  risk_1_field = c("risk_aki_1", "risk_death30_1"),
  rd_field = c("rd_aki", "rd_death30"),
  rr_field = c("rr_aki", "rr_death30"),
  stringsAsFactors = FALSE
)

get_ci <- function(trial, model, field) {
  current <- ci_table[
    ci_table$trial == trial &
      ci_table$model == model &
      ci_table$estimand_field == field,
    ,
    drop = FALSE
  ]
  if (!(
    nrow(current) == 1L &&
      current$n_success == expected_success_n[[trial]] &&
      current$confidence_level == expected_confidence_level &&
      as.integer(current$quantile_type) == expected_quantile_type
  )) {
    stop(trial, " / ", model, " / ", field, "的受控区间异常。")
  }
  c(current$ci_lower, current$ci_upper)
}

numeric_result_tables <- list()
display_tables <- list()

format_percent_ci <- function(estimate, lower, upper) {
  sprintf("%.1f%% (%.1f%% to %.1f%%)", 100 * estimate, 100 * lower, 100 * upper)
}
format_percentage_point_ci <- function(estimate, lower, upper) {
  sprintf("%.1f (%.1f to %.1f)", 100 * estimate, 100 * lower, 100 * upper)
}
format_ratio_ci <- function(estimate, lower, upper) {
  sprintf("%.2f (%.2f to %.2f)", estimate, lower, upper)
}

for (trial in trial_registry$trial) {
  point_row <- point_table[
    point_table$trial == trial & point_table$model == "cbps_ate",
    ,
    drop = FALSE
  ]
  raw_row <- raw_count_table[raw_count_table$trial == trial, , drop = FALSE]
  if (nrow(point_row) != 1L || nrow(raw_row) != 1L) {
    stop(trial, "的主模型点估计或原始计数不唯一。")
  }

  table_rows <- list()
  for (outcome_index in seq_len(nrow(outcome_specification))) {
    specification <- outcome_specification[outcome_index, , drop = FALSE]
    risk_0_ci <- get_ci(trial, "cbps_ate", specification$risk_0_field)
    risk_1_ci <- get_ci(trial, "cbps_ate", specification$risk_1_field)
    rd_ci <- get_ci(trial, "cbps_ate", specification$rd_field)
    rr_ci <- get_ci(trial, "cbps_ate", specification$rr_field)

    table_rows[[length(table_rows) + 1L]] <- data.frame(
      trial = trial,
      outcome = specification$outcome,
      model = "cbps_ate",
      continuation_events = raw_row[[specification$continuation_event_field]],
      continuation_n = raw_row$continuation_n,
      continuation_risk = point_row[[specification$risk_0_field]],
      continuation_ci_lower = risk_0_ci[[1]],
      continuation_ci_upper = risk_0_ci[[2]],
      deescalation_events = raw_row[[specification$deescalation_event_field]],
      deescalation_n = raw_row$deescalation_n,
      deescalation_risk = point_row[[specification$risk_1_field]],
      deescalation_ci_lower = risk_1_ci[[1]],
      deescalation_ci_upper = risk_1_ci[[2]],
      risk_difference = point_row[[specification$rd_field]],
      rd_ci_lower = rd_ci[[1]],
      rd_ci_upper = rd_ci[[2]],
      risk_ratio = point_row[[specification$rr_field]],
      rr_ci_lower = rr_ci[[1]],
      rr_ci_upper = rr_ci[[2]],
      bootstrap_success_n = expected_success_n[[trial]],
      confidence_level = expected_confidence_level,
      quantile_type = expected_quantile_type,
      stringsAsFactors = FALSE
    )
  }

  numeric_table <- do.call(rbind, table_rows)
  numeric_result_tables[[trial]] <- numeric_table
  display_tables[[trial]] <- data.frame(
    Outcome = numeric_table$outcome,
    `Continuation events/N` = paste0(
      numeric_table$continuation_events,
      "/",
      numeric_table$continuation_n
    ),
    `Continuation risk, % (95% CI)` = format_percent_ci(
      numeric_table$continuation_risk,
      numeric_table$continuation_ci_lower,
      numeric_table$continuation_ci_upper
    ),
    `De-escalation events/N` = paste0(
      numeric_table$deescalation_events,
      "/",
      numeric_table$deescalation_n
    ),
    `De-escalation risk, % (95% CI)` = format_percent_ci(
      numeric_table$deescalation_risk,
      numeric_table$deescalation_ci_lower,
      numeric_table$deescalation_ci_upper
    ),
    `Risk difference, percentage points (95% CI)` =
      format_percentage_point_ci(
        numeric_table$risk_difference,
        numeric_table$rd_ci_lower,
        numeric_table$rd_ci_upper
      ),
    `Risk ratio (95% CI)` = format_ratio_ci(
      numeric_table$risk_ratio,
      numeric_table$rr_ci_lower,
      numeric_table$rr_ci_upper
    ),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

all_numeric_results <- do.call(rbind, numeric_result_tables)
if (!(
  nrow(all_numeric_results) == 4L &&
    all(is.finite(as.matrix(all_numeric_results[c(
      "continuation_risk",
      "continuation_ci_lower",
      "continuation_ci_upper",
      "deescalation_risk",
      "deescalation_ci_lower",
      "deescalation_ci_upper",
      "risk_difference",
      "rd_ci_lower",
      "rd_ci_upper",
      "risk_ratio",
      "rr_ci_lower",
      "rr_ci_upper"
    )])))
)) {
  stop("MT-02/MT-03数值表结构或有限值校验失败。")
}

output_dir <- file.path(
  output_root,
  "05_主要结局",
  "controlled_reporting"
)
if (!dir.exists(output_dir) && !dir.create(output_dir, recursive = TRUE)) {
  stop("无法创建受控主要结局输出目录。")
}

output_names <- c(
  "primary_point_estimates_by_imputation_controlled.csv",
  "primary_point_estimates_pooled_controlled.csv",
  "primary_raw_event_counts_controlled.csv",
  "MT-02_anti_mrsa_outcomes_controlled.csv",
  "MT-02_anti_mrsa_outcomes_controlled.docx",
  "MT-03_anti_psa_outcomes_controlled.csv",
  "MT-03_anti_psa_outcomes_controlled.docx",
  "primary_outcome_reporting_audit.csv",
  "_analysis_outputs.md"
)
output_paths <- file.path(output_dir, output_names)
if (any(file.exists(output_paths))) {
  stop("一个或多个07b输出已存在；为保护审计链，拒绝覆盖。")
}

staging_dir <- tempfile(pattern = ".07b_stage_", tmpdir = output_dir)
if (!dir.create(staging_dir)) {
  stop("无法创建07b staging目录。")
}
staged_paths <- file.path(staging_dir, output_names)

utils::write.csv(
  imputation_table,
  staged_paths[[1]],
  row.names = FALSE,
  na = "",
  fileEncoding = "UTF-8"
)
utils::write.csv(
  point_table,
  staged_paths[[2]],
  row.names = FALSE,
  na = "",
  fileEncoding = "UTF-8"
)
utils::write.csv(
  raw_count_table,
  staged_paths[[3]],
  row.names = FALSE,
  na = "",
  fileEncoding = "UTF-8"
)
utils::write.csv(
  numeric_result_tables$anti_mrsa,
  staged_paths[[4]],
  row.names = FALSE,
  na = "",
  fileEncoding = "UTF-8"
)
utils::write.csv(
  numeric_result_tables$anti_psa,
  staged_paths[[6]],
  row.names = FALSE,
  na = "",
  fileEncoding = "UTF-8"
)

make_outcome_flextable <- function(trial, table_id) {
  current <- flextable::flextable(display_tables[[trial]])
  current <- flextable::set_caption(
    current,
    caption = paste0(
      table_id,
      ". Outcomes in the ",
      if (trial == "anti_mrsa") "anti-MRSA" else "anti-PSA",
      " target trial"
    )
  )
  current <- flextable::theme_booktabs(current)
  current <- flextable::fontsize(current, size = 9, part = "all")
  current <- flextable::align(current, j = 1, align = "left", part = "all")
  current <- flextable::align(
    current,
    j = 2:ncol(display_tables[[trial]]),
    align = "center",
    part = "all"
  )
  current <- flextable::autofit(current)
  current <- flextable::add_footer_lines(
    current,
    values = "Risk difference is de-escalation minus continuation."
  )
  bootstrap_note <- if (trial == "anti_mrsa") {
    "Two-sided 95% percentile CIs use 1000 successful Bootstrap replicates."
  } else {
    paste0(
      "Two-sided 95% percentile CIs use 999 successful Bootstrap replicates; ",
      "replicate b=934 was preserved as failed and was not replaced."
    )
  }
  current <- flextable::add_footer_lines(current, values = bootstrap_note)
  flextable::fontsize(current, size = 8, part = "footer")
}

flextable::save_as_docx(
  `MT-02` = make_outcome_flextable("anti_mrsa", "MT-02"),
  path = staged_paths[[5]]
)
flextable::save_as_docx(
  `MT-03` = make_outcome_flextable("anti_psa", "MT-03"),
  path = staged_paths[[7]]
)

reporting_audit <- data.frame(
  trial = c("anti_mrsa", "anti_psa", "overall", "overall", "overall"),
  metric = c(
    "bootstrap_success_n",
    "bootstrap_success_n",
    "analysis_lock",
    "output_plan",
    "anti_psa_preserved_failed_index"
  ),
  value = c("1000", "999", expected_lock, expected_plan, "934"),
  stringsAsFactors = FALSE
)
utils::write.csv(
  reporting_audit,
  staged_paths[[8]],
  row.names = FALSE,
  na = "",
  fileEncoding = "UTF-8"
)

manifest_lines <- c(
  "# Analysis Outputs",
  "",
  "Generated: 2026-08-01",
  "Study type: controlled primary-outcome reporting",
  "Output plan: D-4.6-ANALYSIS-OUTPUT-PLAN-V1.0",
  "Analysis Lock: D-4.7N-ANALYSIS-LOCK-V2.0",
  "",
  "## Main tables",
  "",
  "- `MT-02_anti_mrsa_outcomes_controlled.csv` / `.docx` — anti-MRSA outcomes.",
  "- `MT-03_anti_psa_outcomes_controlled.csv` / `.docx` — anti-PSA outcomes.",
  "",
  "## Traceability data",
  "",
  "- `primary_point_estimates_by_imputation_controlled.csv` — 20-imputation risks.",
  "- `primary_point_estimates_pooled_controlled.csv` — pooled point estimates.",
  "- `primary_raw_event_counts_controlled.csv` — unweighted events and denominators.",
  "- `primary_outcome_reporting_audit.csv` — lock and Bootstrap inclusion audit.",
  "",
  "anti-MRSA uses 1000 successful Bootstrap replicates. anti-PSA uses 999",
  "successful replicates, with b=934 preserved as failed and not replaced."
)
writeLines(enc2utf8(manifest_lines), staged_paths[[9]], useBytes = TRUE)

if (!all(file.exists(staged_paths))) {
  unlink(staging_dir, recursive = TRUE, force = TRUE)
  stop("07b staging输出不完整。")
}

# 写出前确认b=934失败证据未改变且成功checkpoint仍不存在。
if (!(
  identical(
    unname(as.character(tools::md5sum(failure_path))),
    failure_md5_before
  ) &&
    !file.exists(failed_checkpoint_path)
)) {
  unlink(staging_dir, recursive = TRUE, force = TRUE)
  stop("写出前anti-PSA b=934证据保护复核失败。")
}

renamed_paths <- character(0)
for (output_index in seq_along(output_paths)) {
  if (!file.rename(staged_paths[[output_index]], output_paths[[output_index]])) {
    if (length(renamed_paths) > 0L) {
      unlink(renamed_paths, force = TRUE)
    }
    unlink(staging_dir, recursive = TRUE, force = TRUE)
    stop("07b原子重命名失败：", output_paths[[output_index]])
  }
  renamed_paths <- c(renamed_paths, output_paths[[output_index]])
}
unlink(staging_dir, recursive = TRUE, force = TRUE)

if (!(
  identical(
    unname(as.character(tools::md5sum(failure_path))),
    failure_md5_before
  ) &&
    !file.exists(failed_checkpoint_path)
)) {
  stop("写出后anti-PSA b=934证据保护复核失败。")
}

cat("CONTROLLED_PRIMARY_TABLE_STATUS=success\n")
cat("ANTI_MRSA_BOOTSTRAP_SUCCESS_N=1000\n")
cat("ANTI_PSA_BOOTSTRAP_SUCCESS_N=999\n")
cat("ANTI_PSA_PRESERVED_FAILED_INDEX=934\n")
cat("MT02_ROW_N=", nrow(numeric_result_tables$anti_mrsa), "\n", sep = "")
cat("MT03_ROW_N=", nrow(numeric_result_tables$anti_psa), "\n", sep = "")
cat("OUTPUT_DIR=", normalizePath(output_dir, winslash = "/"), "\n", sep = "")
print(all_numeric_results, row.names = FALSE)
