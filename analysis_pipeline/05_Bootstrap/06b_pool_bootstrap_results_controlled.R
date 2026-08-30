# ==============================================================================
# Analysis: anti-MRSA/anti-PSA 受控 Bootstrap 汇总与百分位 95% 置信区间
# Date: 2026-08-01
# Random seed: 不使用随机数
# R: 4.4.2
# Analysis Lock: D-4.7N-ANALYSIS-LOCK-V2.0
#
# Method purpose:
#   汇总 anti-MRSA 的 1000 个成功 Bootstrap checkpoint，以及 anti-PSA 的
#   999 个成功 checkpoint（b=934 保留为失败），计算 type=7 双侧 95%
#   百分位置信区间。脚本不重抽样、不重新插补、不重新拟合模型。
#
# Method sources:
#   1) R Core Team. stats::quantile documentation, type=7.
#      https://stat.ethz.ch/R-manual/R-devel/library/stats/html/quantile.html
#   2) Hyndman RJ, Fan Y. Sample Quantiles in Statistical Packages.
#      The American Statistician. 1996;50(4):361-365.
#      https://doi.org/10.1080/00031305.1996.10473566
#   3) Pawel S, Bartoš F, Siepe BS, Lohmann A. Handling Missingness, Failures,
#      and Non-Convergence in Simulation Studies. The American Statistician.
#      2026;80(1):31-48.
#      https://doi.org/10.1080/00031305.2025.2540002
#
# Citation/version and access date:
#   R 4.4.2 stats::quantile(type=7); sources accessed 2026-08-01.
#
# Adaptation notes:
#   anti-MRSA 的正式纳入集合为 b=1–1000；anti-PSA 的正式纳入集合为
#   b=1–933 与 b=935–1000。b=934 不赋值、不插补、不替代，也不使用
#   na.rm=TRUE 隐藏失败。曲线区间是逐时间点的点态区间，不是同时置信带。
#   所有输出写入新的 controlled_pooling 目录，旧 06 脚本及既有证据只读。
#
# Verification commands:
#   Rscript --vanilla -e \
#     "invisible(parse(file='05_Bootstrap/06b_pool_bootstrap_results_controlled.R'))"
#   Rscript --vanilla 05_Bootstrap/06b_pool_bootstrap_results_controlled.R
#
# RStudio execution:
#   source("05_Bootstrap/06b_pool_bootstrap_results_controlled.R", encoding="UTF-8")
# ==============================================================================

options(stringsAsFactors = FALSE)

# 同时兼容 Rscript 与 RStudio source()，从当前脚本定位项目根目录。
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
  stop("无法定位当前脚本；请使用 Rscript 或 source() 运行本文件。")
}

script_path <- get_current_script_path()
package_root <- dirname(dirname(script_path))
bootstrap_root <- file.path(package_root, "09_输出结果", "04_Bootstrap")
checkpoint_dir <- file.path(bootstrap_root, "checkpoints")
status_dir <- file.path(bootstrap_root, "status")
failure_dir <- file.path(bootstrap_root, "failures")
output_dir <- file.path(bootstrap_root, "controlled_pooling")

expected_lock <- "D-4.7N-ANALYSIS-LOCK-V2.0"
confidence_level <- 0.95
confidence_probs <- c(0.025, 0.975)
quantile_type <- 7L
expected_models <- c("cbps_ate", "ato")
expected_strategies <- c(0L, 1L)
expected_times <- 0:7
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
required_components <- c(
  "pooled_estimates",
  "imputation_estimates",
  "weight_diagnostics",
  "imputation_curves",
  "pooled_curves"
)
optional_audit_numeric_fields <- "ato_brglm_retry_maxit"

trial_plan <- list(
  anti_mrsa = list(
    attempted_indices = seq_len(1000L),
    included_indices = seq_len(1000L),
    failed_indices = integer(0)
  ),
  anti_psa = list(
    attempted_indices = seq_len(1000L),
    included_indices = setdiff(seq_len(1000L), 934L),
    failed_indices = 934L
  )
)

if (!all(dir.exists(c(bootstrap_root, checkpoint_dir, status_dir, failure_dir)))) {
  stop("Bootstrap 输入目录不完整。")
}

# 只检查结果对象中的数值列；未触发 ATO 重试时 maxit=NA 是允许的审计值。
numeric_values_finite <- function(object) {
  if (is.data.frame(object)) {
    numeric_columns <- (
      vapply(object, is.numeric, logical(1)) &
        !names(object) %in% optional_audit_numeric_fields
    )
    values <- unlist(object[numeric_columns], use.names = FALSE)
  } else if (is.list(object)) {
    values <- unlist(lapply(
      object,
      function(element) {
        if (is.data.frame(element)) {
          numeric_columns <- (
            vapply(element, is.numeric, logical(1)) &
              !names(element) %in% optional_audit_numeric_fields
          )
          return(unlist(element[numeric_columns], use.names = FALSE))
        }
        if (is.numeric(element)) {
          return(element)
        }
        numeric(0)
      }
    ), use.names = FALSE)
  } else if (is.numeric(object)) {
    values <- object
  } else {
    values <- numeric(0)
  }
  length(values) == 0L || all(is.finite(values))
}

file_md5_or_missing <- function(path) {
  if (!file.exists(path)) {
    return("MISSING")
  }
  unname(as.character(tools::md5sum(path)))
}

collapse_key_counts <- function(data, key_fields, count_field = NULL) {
  if (!is.data.frame(data) || nrow(data) == 0L) {
    return("none")
  }
  if (!all(key_fields %in% names(data))) {
    stop("审计摘要缺少分组字段：", paste(key_fields, collapse = ", "))
  }
  keys <- do.call(paste, c(data[key_fields], sep = "|"))
  counts <- if (is.null(count_field)) {
    table(keys)
  } else {
    if (!count_field %in% names(data)) {
      stop("审计摘要缺少计数字段：", count_field)
    }
    tapply(data[[count_field]], keys, sum)
  }
  paste(paste(names(counts), as.integer(counts), sep = "="), collapse = "; ")
}

# 合并可能来自不同审计版本、列集合不同的数据框。
row_bind_fill <- function(data_list) {
  if (length(data_list) == 0L) {
    return(data.frame())
  }
  all_columns <- unique(unlist(lapply(data_list, names), use.names = FALSE))
  normalized <- lapply(
    data_list,
    function(current) {
      missing_columns <- setdiff(all_columns, names(current))
      for (column in missing_columns) {
        current[[column]] <- NA
      }
      current[all_columns]
    }
  )
  do.call(rbind, normalized)
}

old_pooling_script <- file.path(
  package_root,
  "05_Bootstrap",
  "06_pool_bootstrap_results.R"
)
if (!file.exists(old_pooling_script)) {
  stop("旧 06_pool_bootstrap_results.R 缺失。")
}
old_pooling_script_md5_before <- file_md5_or_missing(old_pooling_script)

# 旧正式输出属于保护对象；无论当前是否存在，本脚本都不改变其状态。
legacy_output_paths <- file.path(
  bootstrap_root,
  c(
    "bootstrap_pooled_replicates.csv",
    "bootstrap_percentile_ci.csv",
    "bootstrap_curve_replicates.csv",
    "bootstrap_curve_percentile_ci.csv"
  )
)
legacy_output_md5_before <- vapply(
  legacy_output_paths,
  file_md5_or_missing,
  character(1)
)

psa_failed_task_id <- "anti_psa_b0934"
psa_failed_status_path <- file.path(
  status_dir,
  paste0(psa_failed_task_id, "_status.csv")
)
psa_failed_evidence_path <- file.path(
  failure_dir,
  paste0(psa_failed_task_id, "_failure.rds")
)
psa_failed_checkpoint_path <- file.path(
  checkpoint_dir,
  paste0(psa_failed_task_id, "_checkpoint.rds")
)
if (!file.exists(psa_failed_status_path) || !file.exists(psa_failed_evidence_path)) {
  stop("anti-PSA b=934 的 status 或 failure.rds 缺失。")
}
if (file.exists(psa_failed_checkpoint_path)) {
  stop("anti-PSA b=934 不应存在成功 checkpoint。")
}
psa_failed_status_md5_before <- file_md5_or_missing(psa_failed_status_path)
psa_failed_evidence_md5_before <- file_md5_or_missing(psa_failed_evidence_path)

# 失败文件集合必须与已批准方案完全一致。
mrsa_failure_files <- list.files(
  failure_dir,
  pattern = "^anti_mrsa_b[0-9]{4}_failure[.]rds$",
  full.names = FALSE
)
psa_failure_files <- list.files(
  failure_dir,
  pattern = "^anti_psa_b[0-9]{4}_failure[.]rds$",
  full.names = FALSE
)
if (length(mrsa_failure_files) != 0L) {
  stop("anti-MRSA 存在未批准的 failure.rds。")
}
if (!identical(psa_failure_files, "anti_psa_b0934_failure.rds")) {
  stop("anti-PSA failure.rds 集合不是仅含 b=934。")
}

# 逐试验验证全部 status，再读取正式纳入 checkpoint。
collect_trial <- function(trial, plan) {
  status_rows <- vector("list", length(plan$attempted_indices))
  replicate_rows <- vector("list", length(plan$included_indices))
  curve_rows <- vector("list", length(plan$included_indices))
  event_rows <- list()
  retry_rows <- list()
  warning_rows <- list()
  logged_event_total <- 0L
  logreg_retry_total <- 0L
  ato_retry_checkpoint_n <- 0L
  ato_retry_row_n <- 0L

  for (bootstrap_index in plan$attempted_indices) {
    task_id <- paste0(trial, "_b", sprintf("%04d", bootstrap_index))
    status_path <- file.path(status_dir, paste0(task_id, "_status.csv"))
    if (!file.exists(status_path)) {
      stop("缺少 status CSV：", task_id)
    }
    current_status <- utils::read.csv(
      status_path,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    expected_status <- if (bootstrap_index %in% plan$failed_indices) {
      "failed"
    } else {
      "success"
    }
    if (!(
      nrow(current_status) == 1L &&
        identical(current_status$trial, trial) &&
        identical(as.integer(current_status$bootstrap_index), bootstrap_index) &&
        identical(current_status$status, expected_status)
    )) {
      stop("status 身份或状态不符合纳入计划：", task_id)
    }
    status_rows[[bootstrap_index]] <- data.frame(
      trial = trial,
      bootstrap_index = bootstrap_index,
      status = expected_status,
      stringsAsFactors = FALSE
    )
  }

  for (position in seq_along(plan$included_indices)) {
    bootstrap_index <- plan$included_indices[[position]]
    task_id <- paste0(trial, "_b", sprintf("%04d", bootstrap_index))
    checkpoint_path <- file.path(
      checkpoint_dir,
      paste0(task_id, "_checkpoint.rds")
    )
    if (!file.exists(checkpoint_path)) {
      stop("缺少正式纳入 checkpoint：", task_id)
    }
    checkpoint <- readRDS(checkpoint_path)
    if (!(
      identical(checkpoint$status, "success") &&
        identical(checkpoint$trial, trial) &&
        identical(as.integer(checkpoint$bootstrap_index), bootstrap_index) &&
        identical(checkpoint$analysis_lock, expected_lock)
    )) {
      stop("checkpoint 身份或 Analysis Lock 异常：", task_id)
    }
    if (!(
      all(required_components %in% names(checkpoint)) &&
        all(vapply(checkpoint[required_components], length, integer(1)) > 0L)
    )) {
      stop("checkpoint 核心结构不完整：", task_id)
    }
    if (!all(vapply(
      checkpoint[required_components],
      numeric_values_finite,
      logical(1)
    ))) {
      stop("checkpoint 存在非有限关键结果：", task_id)
    }
    if (!(
      "imputation" %in% names(checkpoint) &&
        identical(as.integer(checkpoint$imputation$m), 20L)
    )) {
      stop("checkpoint 的插补次数不是 m=20：", task_id)
    }

    current_estimates <- checkpoint$pooled_estimates
    if (!(
      is.data.frame(current_estimates) &&
        nrow(current_estimates) == 2L &&
        all(c("trial", "bootstrap_index", "model", effect_fields) %in%
          names(current_estimates)) &&
        setequal(current_estimates$model, expected_models) &&
        !anyDuplicated(current_estimates$model) &&
        all(current_estimates$trial == trial) &&
        all(as.integer(current_estimates$bootstrap_index) == bootstrap_index) &&
        all(is.finite(as.matrix(current_estimates[effect_fields])))
    )) {
      stop("pooled_estimates 结构或数值异常：", task_id)
    }
    replicate_rows[[position]] <- current_estimates[
      ,
      c("trial", "bootstrap_index", "model", effect_fields),
      drop = FALSE
    ]

    current_curves <- checkpoint$pooled_curves
    required_curve_columns <- c(
      "trial",
      "bootstrap_index",
      "model",
      "strategy",
      "time_days",
      "aki_risk"
    )
    if (!(
      is.data.frame(current_curves) &&
        nrow(current_curves) == 32L &&
        all(required_curve_columns %in% names(current_curves))
    )) {
      stop("pooled_curves 结构异常：", task_id)
    }
    curve_grid <- expand.grid(
      model = expected_models,
      strategy = expected_strategies,
      time_days = expected_times,
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    )
    current_curve_keys <- do.call(
      paste,
      c(current_curves[c("model", "strategy", "time_days")], sep = "|")
    )
    expected_curve_keys <- do.call(paste, c(curve_grid, sep = "|"))
    if (!(
      all(current_curves$trial == trial) &&
        all(as.integer(current_curves$bootstrap_index) == bootstrap_index) &&
        !anyDuplicated(current_curve_keys) &&
        setequal(current_curve_keys, expected_curve_keys) &&
        all(is.finite(current_curves$aki_risk)) &&
        all(current_curves$aki_risk >= 0 & current_curves$aki_risk <= 1)
    )) {
      stop("pooled_curves 网格、身份或数值异常：", task_id)
    }
    curve_rows[[position]] <- current_curves[required_curve_columns]

    current_logged_event_n <- if (
      length(checkpoint$mice_logged_event_n) == 1L
    ) {
      as.integer(checkpoint$mice_logged_event_n)
    } else {
      0L
    }
    logged_event_total <- logged_event_total + current_logged_event_n
    if (
      is.data.frame(checkpoint$mice_logged_event_summary) &&
        nrow(checkpoint$mice_logged_event_summary) > 0L
    ) {
      current_events <- checkpoint$mice_logged_event_summary
      current_events$trial <- trial
      current_events$bootstrap_index <- bootstrap_index
      event_rows[[length(event_rows) + 1L]] <- current_events
    }
    if (length(checkpoint$mice_warning_messages) > 0L) {
      warning_rows[[length(warning_rows) + 1L]] <- data.frame(
        trial = trial,
        bootstrap_index = bootstrap_index,
        warning_message = as.character(checkpoint$mice_warning_messages),
        stringsAsFactors = FALSE
      )
    }

    current_retry_n <- if (length(checkpoint$mice_logreg_retry_n) == 1L) {
      as.integer(checkpoint$mice_logreg_retry_n)
    } else {
      0L
    }
    current_retry_summary_n <- if (
      is.data.frame(checkpoint$mice_logreg_retry_summary)
    ) {
      nrow(checkpoint$mice_logreg_retry_summary)
    } else {
      0L
    }
    if (current_retry_n != current_retry_summary_n) {
      stop("MICE logreg 重试计数与摘要不一致：", task_id)
    }
    logreg_retry_total <- logreg_retry_total + current_retry_n
    if (current_retry_summary_n > 0L) {
      current_retries <- checkpoint$mice_logreg_retry_summary
      if (!"target_field" %in% names(current_retries)) {
        current_retries$target_field <- "nephrotoxic_drug_exposure"
      }
      if (!"locked_method" %in% names(current_retries)) {
        current_retries$locked_method <- "logreg"
      }
      current_retries$trial <- trial
      current_retries$bootstrap_index <- bootstrap_index
      retry_rows[[length(retry_rows) + 1L]] <- current_retries
    }

    ato_retry_flags <- if (
      "ato_brglm_retry_used" %in% names(checkpoint$weight_diagnostics)
    ) {
      checkpoint$weight_diagnostics$ato_brglm_retry_used %in% TRUE
    } else {
      logical(0)
    }
    ato_retry_checkpoint_n <- (
      ato_retry_checkpoint_n + as.integer(any(ato_retry_flags))
    )
    ato_retry_row_n <- ato_retry_row_n + sum(ato_retry_flags)
  }

  status_table <- do.call(rbind, status_rows)
  replicate_table <- do.call(rbind, replicate_rows)
  curve_table <- do.call(rbind, curve_rows)
  event_table <- row_bind_fill(event_rows)
  retry_table <- row_bind_fill(retry_rows)
  warning_table <- row_bind_fill(warning_rows)

  if (nrow(replicate_table) != length(plan$included_indices) * 2L) {
    stop(trial, " 标量重复行数错误。")
  }
  if (nrow(curve_table) != length(plan$included_indices) * 32L) {
    stop(trial, " 曲线重复行数错误。")
  }

  list(
    status = status_table,
    replicates = replicate_table,
    curves = curve_table,
    events = event_table,
    retries = retry_table,
    warnings = warning_table,
    logged_event_total = logged_event_total,
    logreg_retry_total = logreg_retry_total,
    ato_retry_checkpoint_n = ato_retry_checkpoint_n,
    ato_retry_row_n = ato_retry_row_n
  )
}

trial_results <- lapply(
  names(trial_plan),
  function(trial) collect_trial(trial, trial_plan[[trial]])
)
names(trial_results) <- names(trial_plan)

replicate_table <- do.call(rbind, lapply(trial_results, `[[`, "replicates"))
curve_replicate_table <- do.call(rbind, lapply(trial_results, `[[`, "curves"))

expected_replicate_rows <- c(anti_mrsa = 2000L, anti_psa = 1998L)
expected_curve_rows <- c(anti_mrsa = 32000L, anti_psa = 31968L)
for (trial in names(trial_plan)) {
  if (nrow(trial_results[[trial]]$replicates) != expected_replicate_rows[[trial]]) {
    stop(trial, " 的标量重复总行数不符合 spec。")
  }
  if (nrow(trial_results[[trial]]$curves) != expected_curve_rows[[trial]]) {
    stop(trial, " 的曲线重复总行数不符合 spec。")
  }
}

# 按正式纳入集合计算标量 type=7 双侧 95% 百分位区间。
interval_rows <- list()
for (trial in names(trial_plan)) {
  expected_indices <- trial_plan[[trial]]$included_indices
  for (model in expected_models) {
    current <- replicate_table[
      replicate_table$trial == trial & replicate_table$model == model,
      ,
      drop = FALSE
    ]
    if (!(
      nrow(current) == length(expected_indices) &&
        !anyDuplicated(current$bootstrap_index) &&
        setequal(as.integer(current$bootstrap_index), expected_indices)
    )) {
      stop(trial, " ", model, " 的标量 Bootstrap 集合不完整。")
    }
    for (field in effect_fields) {
      values <- current[[field]]
      if (
        length(values) != length(expected_indices) ||
          anyNA(values) ||
          !all(is.finite(values))
      ) {
        stop(trial, " ", model, " ", field, " 不是完整有限向量。")
      }
      current_ci <- as.numeric(stats::quantile(
        values,
        probs = confidence_probs,
        type = quantile_type,
        na.rm = FALSE,
        names = FALSE
      ))
      interval_rows[[length(interval_rows) + 1L]] <- data.frame(
        trial = trial,
        model = model,
        estimand_field = field,
        n_success = length(values),
        confidence_level = confidence_level,
        quantile_type = quantile_type,
        ci_lower = current_ci[[1]],
        ci_upper = current_ci[[2]],
        stringsAsFactors = FALSE
      )
    }
  }
}
interval_table <- do.call(rbind, interval_rows)

# 计算 AKI 曲线的逐时间点点态百分位区间。
curve_interval_rows <- list()
for (trial in names(trial_plan)) {
  expected_n <- length(trial_plan[[trial]]$included_indices)
  for (model in expected_models) {
    for (strategy in expected_strategies) {
      for (time_days in expected_times) {
        values <- curve_replicate_table$aki_risk[
          curve_replicate_table$trial == trial &
            curve_replicate_table$model == model &
            as.integer(curve_replicate_table$strategy) == strategy &
            curve_replicate_table$time_days == time_days
        ]
        if (
          length(values) != expected_n ||
            anyNA(values) ||
            !all(is.finite(values))
        ) {
          stop(
            trial,
            " ",
            model,
            " strategy=",
            strategy,
            " time=",
            time_days,
            " 的曲线 Bootstrap 集合不完整。"
          )
        }
        current_ci <- as.numeric(stats::quantile(
          values,
          probs = confidence_probs,
          type = quantile_type,
          na.rm = FALSE,
          names = FALSE
        ))
        curve_interval_rows[[length(curve_interval_rows) + 1L]] <- data.frame(
          trial = trial,
          model = model,
          strategy = strategy,
          time_days = time_days,
          curve_quantity = "aki_risk_pointwise",
          simultaneous_band = FALSE,
          n_success = length(values),
          confidence_level = confidence_level,
          quantile_type = quantile_type,
          ci_lower = current_ci[[1]],
          ci_upper = current_ci[[2]],
          stringsAsFactors = FALSE
        )
      }
    }
  }
}
curve_interval_table <- do.call(rbind, curve_interval_rows)

if (nrow(interval_table) != 32L) {
  stop("标量区间表不是 32 行。")
}
if (nrow(curve_interval_table) != 64L) {
  stop("曲线区间表不是 64 行。")
}
if (!(
  all(is.finite(interval_table$ci_lower)) &&
    all(is.finite(interval_table$ci_upper)) &&
    all(interval_table$ci_lower <= interval_table$ci_upper) &&
    all(is.finite(curve_interval_table$ci_lower)) &&
    all(is.finite(curve_interval_table$ci_upper)) &&
    all(curve_interval_table$ci_lower <= curve_interval_table$ci_upper)
)) {
  stop("正式区间表包含非有限值或端点次序错误。")
}

# 汇总各试验的审计事件与重试。
audit_rows <- list()
for (trial in names(trial_plan)) {
  current <- trial_results[[trial]]
  event_summary <- collapse_key_counts(
    current$events,
    c("dep", "meth", "out"),
    "event_n"
  )
  warning_summary <- collapse_key_counts(
    current$warnings,
    "warning_message"
  )
  retry_summary <- collapse_key_counts(
    current$retries,
    "target_field"
  )
  trial_metrics <- c(
    attempted_n = length(trial_plan[[trial]]$attempted_indices),
    success_n = length(trial_plan[[trial]]$included_indices),
    preserved_failed_n = length(trial_plan[[trial]]$failed_indices),
    preserved_failed_indices = if (length(trial_plan[[trial]]$failed_indices)) {
      paste(trial_plan[[trial]]$failed_indices, collapse = ";")
    } else {
      "none"
    },
    included_indices = if (trial == "anti_mrsa") {
      "1-1000"
    } else {
      "1-933;935-1000"
    },
    mice_logged_event_n = current$logged_event_total,
    mice_logged_event_summary = event_summary,
    mice_warning_message_n = nrow(current$warnings),
    mice_warning_summary = warning_summary,
    mice_logreg_retry_n = current$logreg_retry_total,
    mice_logreg_retry_summary = retry_summary,
    ato_retry_checkpoint_n = current$ato_retry_checkpoint_n,
    ato_retry_imputation_model_row_n = current$ato_retry_row_n
  )
  audit_rows[[length(audit_rows) + 1L]] <- data.frame(
    trial = trial,
    metric = names(trial_metrics),
    value = unname(trial_metrics),
    stringsAsFactors = FALSE
  )
}

global_metrics <- c(
  analysis_lock = expected_lock,
  r_version = R.version.string,
  confidence_level = confidence_level,
  quantile_type = quantile_type,
  scalar_replicate_row_n = nrow(replicate_table),
  scalar_interval_row_n = nrow(interval_table),
  curve_replicate_row_n = nrow(curve_replicate_table),
  curve_interval_row_n = nrow(curve_interval_table),
  curve_interval_scope = "pointwise_only_not_simultaneous_band",
  anti_psa_rank_sensitivity = "09_输出结果/04_Bootstrap/rank_sensitivity",
  b0934_status_md5_unchanged = TRUE,
  b0934_failure_md5_unchanged = TRUE,
  b0934_success_checkpoint_absent = TRUE,
  old_06_script_unchanged = TRUE,
  legacy_pooling_outputs_unchanged = TRUE,
  pooling_authorized = TRUE,
  execution_status = "completed",
  method_source_r_quantile =
    "https://stat.ethz.ch/R-manual/R-devel/library/stats/html/quantile.html",
  method_source_hyndman_fan =
    "https://doi.org/10.1080/00031305.1996.10473566",
  failure_handling_source_pawel_et_al =
    "https://doi.org/10.1080/00031305.2025.2540002",
  source_access_date = "2026-08-01"
)
audit_rows[[length(audit_rows) + 1L]] <- data.frame(
  trial = "overall",
  metric = names(global_metrics),
  value = unname(global_metrics),
  stringsAsFactors = FALSE
)
audit_table <- do.call(rbind, audit_rows)

summary_lines <- c(
  "anti-MRSA/anti-PSA 受控 Bootstrap 池化摘要",
  "===========================================",
  paste0("Analysis Lock: ", expected_lock),
  "anti-MRSA：纳入 b=1–1000，共 1000 个成功重复。",
  "anti-PSA：纳入 b=1–933 与 b=935–1000，共 999 个成功重复。",
  "anti-PSA b=934：保留 failed、failure.rds 存在且无成功 checkpoint。",
  "置信区间：双侧 95% Bootstrap 百分位区间，stats::quantile(type=7, na.rm=FALSE)。",
  "曲线区间：逐时间点点态区间，不是同时置信带。",
  paste0("标量重复行数：", nrow(replicate_table), "。"),
  paste0("标量区间行数：", nrow(interval_table), "。"),
  paste0("曲线重复行数：", nrow(curve_replicate_table), "。"),
  paste0("曲线区间行数：", nrow(curve_interval_table), "。"),
  paste0(
    "anti-MRSA MICE logged events / logreg retries / ATO retry checkpoints：",
    trial_results$anti_mrsa$logged_event_total,
    " / ",
    trial_results$anti_mrsa$logreg_retry_total,
    " / ",
    trial_results$anti_mrsa$ato_retry_checkpoint_n,
    "。"
  ),
  paste0(
    "anti-PSA MICE logged events / logreg retries / ATO retry checkpoints：",
    trial_results$anti_psa$logged_event_total,
    " / ",
    trial_results$anti_psa$logreg_retry_total,
    " / ",
    trial_results$anti_psa$ato_retry_checkpoint_n,
    "。"
  ),
  "证据保护：b=934 status/failure.rds、旧 06 脚本及旧正式输出均未改变。",
  "",
  "方法依据：",
  "- R Core Team, stats::quantile(type=7): https://stat.ethz.ch/R-manual/R-devel/library/stats/html/quantile.html",
  "- Hyndman RJ, Fan Y (1996): https://doi.org/10.1080/00031305.1996.10473566",
  "- Pawel S, Bartoš F, Siepe BS, Lohmann A (2026): https://doi.org/10.1080/00031305.2025.2540002"
)

manifest_lines <- c(
  "# Analysis Outputs",
  "",
  "Generated: 2026-08-01",
  "Study type: controlled Bootstrap percentile confidence intervals",
  "",
  "## Tables",
  "",
  "- `bootstrap_percentile_ci_controlled.csv` — 标量 type=7 双侧 95% 百分位区间。",
  "- `bootstrap_curve_percentile_ci_controlled.csv` — AKI 曲线逐时间点点态区间。",
  "- `bootstrap_pooling_audit.csv` — 纳入集合、失败证据与重试审计。",
  "",
  "## Data",
  "",
  "- `bootstrap_replicates_controlled.csv` — 正式纳入的标量 Bootstrap 重复。",
  "- `bootstrap_curve_replicates_controlled.csv` — 正式纳入的曲线 Bootstrap 重复。",
  "",
  "## Documentation",
  "",
  "- `bootstrap_pooling_summary.txt` — 中文方法与审计摘要。",
  "",
  "anti-MRSA uses 1000 successful replicates; anti-PSA uses 999 successful",
  "replicates with b=934 preserved as failed. Curve intervals are pointwise."
)

# 所有保护对象在写出前必须保持不变。
if (!(
  identical(file_md5_or_missing(psa_failed_status_path), psa_failed_status_md5_before) &&
    identical(file_md5_or_missing(psa_failed_evidence_path), psa_failed_evidence_md5_before) &&
    !file.exists(psa_failed_checkpoint_path) &&
    identical(file_md5_or_missing(old_pooling_script), old_pooling_script_md5_before) &&
    identical(
      vapply(legacy_output_paths, file_md5_or_missing, character(1)),
      legacy_output_md5_before
    )
)) {
  stop("写出前保护对象发生变化。")
}

if (!dir.exists(output_dir) && !dir.create(output_dir, recursive = TRUE)) {
  stop("无法创建 controlled_pooling 输出目录。")
}
output_names <- c(
  "bootstrap_replicates_controlled.csv",
  "bootstrap_percentile_ci_controlled.csv",
  "bootstrap_curve_replicates_controlled.csv",
  "bootstrap_curve_percentile_ci_controlled.csv",
  "bootstrap_pooling_audit.csv",
  "bootstrap_pooling_summary.txt",
  "_analysis_outputs.md"
)
output_paths <- file.path(output_dir, output_names)
if (any(file.exists(output_paths))) {
  stop("一个或多个受控池化输出已存在；为保护审计链，拒绝覆盖。")
}

# 先在同一文件系统的 staging 目录生成全部文件，再逐个原子重命名。
staging_dir <- tempfile(pattern = ".06b_stage_", tmpdir = output_dir)
if (!dir.create(staging_dir)) {
  stop("无法创建 staging 目录。")
}
staged_paths <- file.path(staging_dir, output_names)

utils::write.csv(
  replicate_table,
  staged_paths[[1]],
  row.names = FALSE,
  na = "",
  fileEncoding = "UTF-8"
)
utils::write.csv(
  interval_table,
  staged_paths[[2]],
  row.names = FALSE,
  na = "",
  fileEncoding = "UTF-8"
)
utils::write.csv(
  curve_replicate_table,
  staged_paths[[3]],
  row.names = FALSE,
  na = "",
  fileEncoding = "UTF-8"
)
utils::write.csv(
  curve_interval_table,
  staged_paths[[4]],
  row.names = FALSE,
  na = "",
  fileEncoding = "UTF-8"
)
utils::write.csv(
  audit_table,
  staged_paths[[5]],
  row.names = FALSE,
  na = "",
  fileEncoding = "UTF-8"
)
writeLines(enc2utf8(summary_lines), staged_paths[[6]], useBytes = TRUE)
writeLines(enc2utf8(manifest_lines), staged_paths[[7]], useBytes = TRUE)

if (!all(file.exists(staged_paths))) {
  unlink(staging_dir, recursive = TRUE, force = TRUE)
  stop("staging 输出不完整。")
}
renamed_paths <- character(0)
for (output_index in seq_along(output_paths)) {
  if (!file.rename(staged_paths[[output_index]], output_paths[[output_index]])) {
    if (length(renamed_paths) > 0L) {
      unlink(renamed_paths, force = TRUE)
    }
    unlink(staging_dir, recursive = TRUE, force = TRUE)
    stop("原子重命名失败：", output_paths[[output_index]])
  }
  renamed_paths <- c(renamed_paths, output_paths[[output_index]])
}
unlink(staging_dir, recursive = TRUE, force = TRUE)

# 写出后再次验证旧证据、旧脚本与旧正式输出未被改变。
if (!(
  identical(file_md5_or_missing(psa_failed_status_path), psa_failed_status_md5_before) &&
    identical(file_md5_or_missing(psa_failed_evidence_path), psa_failed_evidence_md5_before) &&
    !file.exists(psa_failed_checkpoint_path) &&
    identical(file_md5_or_missing(old_pooling_script), old_pooling_script_md5_before) &&
    identical(
      vapply(legacy_output_paths, file_md5_or_missing, character(1)),
      legacy_output_md5_before
    )
)) {
  stop("写出后保护对象复核失败。")
}

cat("CONTROLLED_POOLING_STATUS=success\n")
cat("ANTI_MRSA_SUCCESS_N=1000\n")
cat("ANTI_PSA_SUCCESS_N=999\n")
cat("ANTI_PSA_PRESERVED_FAILED_INDEX=934\n")
cat("CONFIDENCE_LEVEL=0.95\n")
cat("QUANTILE_TYPE=7\n")
cat("SCALAR_INTERVAL_ROW_N=", nrow(interval_table), "\n", sep = "")
cat("CURVE_INTERVAL_ROW_N=", nrow(curve_interval_table), "\n", sep = "")
cat("OUTPUT_DIR=", normalizePath(output_dir, winslash = "/"), "\n", sep = "")
print(interval_table, row.names = FALSE)
