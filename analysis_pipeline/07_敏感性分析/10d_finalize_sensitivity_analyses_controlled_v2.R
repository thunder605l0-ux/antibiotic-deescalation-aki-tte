# =============================================================================
# ST-08 / SF-09 受控V2最终汇总器
#
# 代码来源：
# 1. Efron B, Tibshirani RJ. An Introduction to the Bootstrap. 1993.
# 2. R stats::quantile(type = 7)，R 4.4.2；访问日期：2026-08-02。
# 3. Wickham H. ggplot2: Elegant Graphics for Data Analysis. Springer; 2016；
#    ggplot2 4.0.2，访问日期：2026-08-02。
# 4. 项目锁定脚本 07_敏感性分析/10_run_sensitivity_analyses.R，
#    MD5=9c20d440972610c41e193cc7bf31eab1。
# 5. 项目受控汇总脚本10b，
#    MD5=91245d6f22b53e99b968c1f33b161e42。
# 6. Imai K, Ratkovic M. Covariate Balancing Propensity Score. Journal of the
#    Royal Statistical Society: Series B. 2014;76(1):243-263.
#    doi:10.1111/rssb.12027。
# 7. WeightIt 1.7.0.9004，Git commit
#    595d4194545ec5f55b884ef6f90d3b6389747d27；访问日期：2026-08-03。
# 8. 受控扩展规范：docs/superpowers/specs/
#    2026-08-03-st08-cbps-single-retry-controlled-v2-design.md。
# 9. b=79 logreg100诊断与救援规范：docs/superpowers/specs/
#    2026-08-04-st08-b0079-logreg100-diagnostic-and-rescue-design.md。
# 10. b=373零风险比边界规范：docs/superpowers/specs/
#    2026-08-04-st08-b0373-zero-risk-ratio-boundary-design.md。
# 11. b=584 race_white有序对比列线性依赖规范：docs/superpowers/specs/
#    2026-08-05-st08-b0584-race-white-ordered-contrast-lindep-design.md。
# 12. mice 3.19.0 remove.lindep()/updateLog()；amices/mice commit
#    61083e667fa41cc4bc655492591b100bf8a7e443；访问日期：2026-08-05。
#
# 适配说明：不改变点估计、Bootstrap重复、种子或区间算法；仅在最终汇总前
# 强制验证逐重复V2审计，并用ggplot2绘制SF-09。仅接受规范中逐文件、逐哈希
# 登记的两份历史CBPS失败证据、b=79 logreg100失败证据及b=373零RR失败证据；
# 对应正式重复必须由受控规则产生并通过审计。anti-MRSA使用1000个正式重复；
# anti-PSA使用999个成功重复，b=934保持失败且不替换。
#
# 验证命令：
# Rscript 07_敏感性分析/10d_finalize_sensitivity_analyses_controlled_v2.R
# =============================================================================

# b=635事件路由修正规范：docs/superpowers/specs/
# 2026-08-05-st08-b0635-event-routing-correction-design.md。
get_script_path_st08_finalizer <- function() {
  command_line <- commandArgs(trailingOnly = FALSE)
  file_argument <- grep("^--file=", command_line, value = TRUE)
  if (length(file_argument) == 1L) {
    return(normalizePath(
      sub("^--file=", "", file_argument),
      winslash = "/",
      mustWork = TRUE
    ))
  }
  stop("无法定位10d脚本；请使用Rscript运行。")
}

finalizer_path <- get_script_path_st08_finalizer()
package_root <- normalizePath(
  file.path(dirname(finalizer_path), ".."),
  winslash = "/",
  mustWork = TRUE
)
worker_path <- file.path(
  package_root,
  "07_敏感性分析",
  "10c_run_sensitivity_analyses_controlled_v2.R"
)

old_library_only <- Sys.getenv(
  "ST08_CONTROLLED_V2_LIBRARY_ONLY",
  unset = NA_character_
)
Sys.setenv(ST08_CONTROLLED_V2_LIBRARY_ONLY = "1")
source(worker_path, chdir = TRUE, encoding = "UTF-8")
if (is.na(old_library_only)) {
  Sys.unsetenv("ST08_CONTROLLED_V2_LIBRARY_ONLY")
} else {
  Sys.setenv(ST08_CONTROLLED_V2_LIBRARY_ONLY = old_library_only)
}

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("缺少ggplot2，不能生成SF-09。")
}

# 代码来源：ggplot2官方geom_errorbar()参考文档；ggplot2 4.0.2；
# 访问日期：2026-08-06。
# 适配说明：仅将已弃用的水平误差线包装函数迁移至orientation接口；
# 不改变数据、点估计、Bootstrap重复或区间算法。
# 验证命令：Rscript 07_敏感性分析/10d_finalize_sensitivity_analyses_controlled_v2.R
reviewed_finalizer_warning_failure_path_v2 <- file.path(
  sensitivity_dir,
  "diagnostic_history",
  paste0(
    "ST-08_controlled_v2_finalization_failure_geom_errorbarh_",
    "20260806_204833.rds"
  )
)
expected_finalizer_warning_failure_md5_v2 <-
  "7696469b8611b5962b1942e4ebb5d187"
is_reviewed_finalizer_warning_failure_v2 <- function(path) {
  if (!file.exists(path)) {
    return(FALSE)
  }
  if (!identical(
    tolower(unname(tools::md5sum(path))),
    expected_finalizer_warning_failure_md5_v2
  )) {
    return(FALSE)
  }
  evidence <- tryCatch(
    readRDS(path),
    error = function(condition) NULL
  )
  if (!is.list(evidence)) {
    return(FALSE)
  }
  identical(evidence$status, "failed") &&
    identical(evidence$stage, "controlled_v2_finalization") &&
    identical(
      evidence$message,
      paste0(
        "V2最终汇总出现未分类警告：`geom_errorbarh()` was deprecated ",
        "in ggplot2 4.0.0.\n",
        "ℹ Please use the `orientation` argument of `geom_errorbar()` ",
        "instead."
      )
    ) &&
    identical(
      evidence$controlled_v2_worker_md5,
      "41b534755f84ebc4d9668af3922b54c6"
    ) &&
    identical(
      evidence$controlled_v2_finalizer_md5,
      "76399b12b0941ed0fb05150d7e0d92d9"
    ) &&
    identical(evidence$timestamp, "2026-08-06 20:48:33 CST")
}
if (!is_reviewed_finalizer_warning_failure_v2(
  reviewed_finalizer_warning_failure_path_v2
)) {
  stop("ST-08终审ggplot2兼容性失败证据缺失或不匹配。")
}

finalizer_failure_path_v2 <- file.path(
  failure_dir_v2,
  "ST-08_controlled_v2_finalization_failure.rds"
)
if (file.exists(finalizer_failure_path_v2)) {
  stop("存在尚未审核的V2最终汇总失败证据：", finalizer_failure_path_v2)
}

finalize_st08_v2 <- function() {
  if (!identical(
    tolower(unname(tools::md5sum(reviewed_b0031_failure_path_v2))),
    expected_md5_v2[["anti_mrsa_b0031_failure"]]
  )) {
    stop("anti-MRSA b=31原失败证据MD5不匹配。")
  }
  if (!identical(
    tolower(unname(tools::md5sum(preserved_failure_path))),
    expected_md5_v2[["anti_psa_b0934_failure"]]
  )) {
    stop("anti-PSA b=934保留失败证据MD5不匹配。")
  }
  if (!is_reviewed_b0079_failure_v2(reviewed_b0079_failure_path_v2)) {
    stop("anti-PSA b=79原始logreg100失败证据不匹配。")
  }
  if (!is_verified_b0079_diagnostic_v2(
    reviewed_b0079_diagnostic_path_v2
  )) {
    stop("anti-PSA b=79 logreg100诊断证据不匹配。")
  }
  if (!is_reviewed_b0373_failure_v2(reviewed_b0373_failure_path_v2)) {
    stop("anti-MRSA b=373原始零RR失败证据不匹配。")
  }
  if (!is_verified_b0373_diagnostic_v2(
    reviewed_b0373_diagnostic_path_v2
  )) {
    stop("anti-MRSA b=373零RR诊断证据不匹配。")
  }
  if (!is_reviewed_b0584_failure_v2(reviewed_b0584_failure_path_v2)) {
    stop("anti-MRSA b=584原始MICE事件失败证据不匹配。")
  }
  if (!is_verified_b0584_diagnostic_v2(
    reviewed_b0584_diagnostic_path_v2
  )) {
    stop("anti-MRSA b=584只读MICE事件诊断证据不匹配。")
  }
  if (!is_verified_b0584_preflight_v2(
    reviewed_b0584_preflight_path_v2
  )) {
    stop("anti-MRSA b=584候选精确白名单预检不匹配。")
  }
  if (!is_reviewed_b0635_failure_v2(reviewed_b0635_failure_path_v2)) {
    stop("anti-MRSA b=635原始MICE事件失败证据不匹配。")
  }
  if (!is_verified_b0635_diagnostic_v2(
    reviewed_b0635_diagnostic_path_v2
  )) {
    stop("anti-MRSA b=635事件路由只读诊断证据不匹配。")
  }
  if (
    file.exists(forbidden_success_path) ||
      file.exists(file.path(repeat_dir, "anti_psa_b0934_sensitivity.rds"))
  ) {
    stop("anti-PSA b=934出现禁止的成功对象。")
  }

  pending_v2_failures <- list.files(
    failure_dir_v2,
    pattern = "\\.rds$",
    full.names = TRUE
  )
  reviewed_failure_paths <- normalizePath(
    c(
      reviewed_cbps_failure_registry_v2$path,
      reviewed_b0079_failure_path_v2,
      reviewed_b0373_failure_path_v2,
      reviewed_b0584_failure_path_v2,
      reviewed_b0635_failure_path_v2
    ),
    winslash = "/",
    mustWork = TRUE
  )
  pending_failure_paths <- if (length(pending_v2_failures) == 0L) {
    character(0)
  } else {
    normalizePath(
      pending_v2_failures,
      winslash = "/",
      mustWork = TRUE
    )
  }
  unreviewed_v2_failures <- setdiff(
    pending_failure_paths,
    reviewed_failure_paths
  )
  if (length(unreviewed_v2_failures) > 0L) {
    stop(
      "存在未审核V2失败证据；不得汇总：",
      paste(unreviewed_v2_failures, collapse = " | ")
    )
  }

  ensure_existing_repeat_audits_v2()

  audit_index_rows <- list()
  race_retry_rows <- list()
  existing_retry_rows <- list()
  logreg100_escalation_rows <- list()
  cbps_retry_rows <- list()
  zero_rr_boundary_rows <- list()
  race_white_lindep_rows <- list()
  logged_event_rows <- list()
  event_route_rows <- list()
  timing_rows <- list()

  for (trial in names(formal_bootstrap_indices)) {
    current_indices <- formal_bootstrap_indices[[trial]]
    if (!identical(length(current_indices), expected_replicate_n[[trial]])) {
      stop(trial, "正式编号数与预期不一致。")
    }
    for (bootstrap_index in current_indices) {
      repeat_path <- file.path(
        repeat_dir,
        sprintf("%s_b%04d_sensitivity.rds", trial, bootstrap_index)
      )
      audit_path <- file.path(
        repeat_audit_dir_v2,
        sprintf("%s_b%04d_sensitivity_audit.rds", trial, bootstrap_index)
      )
      if (!file.exists(repeat_path) || !file.exists(audit_path)) {
        stop("缺少正式repeat或V2审计：", trial, " b=", bootstrap_index)
      }
      result <- readRDS(repeat_path)
      audit <- readRDS(audit_path)
      validate_repeat_result(result, trial, bootstrap_index)
      validate_repeat_audit_v2(audit, trial, bootstrap_index)
      current_result_md5 <- unname(tools::md5sum(repeat_path))
      if (!identical(tolower(audit$result_md5), tolower(current_result_md5))) {
        stop("repeat与审计MD5不一致：", trial, " b=", bootstrap_index)
      }

      audit_index_rows[[length(audit_index_rows) + 1L]] <- data.frame(
        trial = trial,
        bootstrap_index = bootstrap_index,
        execution_version = audit$execution_version,
        result_md5 = current_result_md5,
        result_specification_n = as.integer(audit$result_specification_n),
        all_logged_events_allowed = isTRUE(audit$all_logged_events_allowed),
        race_white_lindep_event_n = if (is.null(
          audit$race_white_lindep_event_n
        )) {
          0L
        } else {
          as.integer(audit$race_white_lindep_event_n)
        },
        all_race_white_lindep_events_verified = if (is.null(
          audit$all_race_white_lindep_events_verified
        )) {
          TRUE
        } else {
          isTRUE(audit$all_race_white_lindep_events_verified)
        },
        race_white_retry_n = nrow(audit$race_white_retry_summary),
        existing_logreg_retry_n = nrow(audit$existing_logreg_retry_summary),
        logreg100_escalation_n = if (
          is.null(audit$logreg100_escalation_summary)
        ) {
          0L
        } else {
          nrow(audit$logreg100_escalation_summary)
        },
        cbps_retry_n = if (is.null(audit$cbps_retry_summary)) {
          0L
        } else {
          nrow(audit$cbps_retry_summary)
        },
        zero_rr_boundary_n = if (is.null(
          audit$zero_rr_boundary_summary
        )) {
          0L
        } else {
          nrow(audit$zero_rr_boundary_summary)
        },
        all_zero_rr_boundaries_verified = if (is.null(
          audit$all_zero_rr_boundaries_verified
        )) {
          TRUE
        } else {
          isTRUE(audit$all_zero_rr_boundaries_verified)
        },
        mice_logged_event_n = nrow(audit$mice_logged_events),
        mice_event_route_existing_v2_n = if (is.null(
          audit$mice_logged_event_route
        )) {
          0L
        } else {
          sum(audit$mice_logged_event_route == "existing_v2")
        },
        mice_event_route_race_white_lindep_v3_n = if (is.null(
          audit$mice_logged_event_route
        )) {
          0L
        } else {
          sum(audit$mice_logged_event_route == "race_white_lindep_v3")
        },
        mice_event_route_disallowed_n = if (is.null(
          audit$mice_logged_event_route
        )) {
          0L
        } else {
          sum(audit$mice_logged_event_route == "disallowed")
        },
        elapsed_seconds = if (is.null(audit$elapsed_seconds)) {
          NA_real_
        } else {
          as.numeric(audit$elapsed_seconds)
        },
        stringsAsFactors = FALSE
      )

      if (nrow(audit$race_white_retry_summary) > 0L) {
        current <- audit$race_white_retry_summary
        current$trial <- trial
        current$bootstrap_index <- bootstrap_index
        race_retry_rows[[length(race_retry_rows) + 1L]] <- current
      }
      if (nrow(audit$existing_logreg_retry_summary) > 0L) {
        current <- audit$existing_logreg_retry_summary
        current$trial <- trial
        current$bootstrap_index <- bootstrap_index
        existing_retry_rows[[length(existing_retry_rows) + 1L]] <- current
      }
      if (
        !is.null(audit$logreg100_escalation_summary) &&
          nrow(audit$logreg100_escalation_summary) > 0L
      ) {
        current <- audit$logreg100_escalation_summary
        current$trial <- trial
        current$bootstrap_index <- bootstrap_index
        logreg100_escalation_rows[[
          length(logreg100_escalation_rows) + 1L
        ]] <- current
      }
      if (
        !is.null(audit$cbps_retry_summary) &&
          nrow(audit$cbps_retry_summary) > 0L
      ) {
        current <- audit$cbps_retry_summary
        current$trial <- trial
        current$bootstrap_index <- bootstrap_index
        cbps_retry_rows[[length(cbps_retry_rows) + 1L]] <- current
      }
      if (
        !is.null(audit$zero_rr_boundary_summary) &&
          nrow(audit$zero_rr_boundary_summary) > 0L
      ) {
        zero_rr_boundary_rows[[
          length(zero_rr_boundary_rows) + 1L
        ]] <- audit$zero_rr_boundary_summary
      }
      if (
        !is.null(audit$race_white_lindep_summary) &&
          nrow(audit$race_white_lindep_summary) > 0L
      ) {
        race_white_lindep_rows[[
          length(race_white_lindep_rows) + 1L
        ]] <- audit$race_white_lindep_summary
      }
      if (nrow(audit$mice_logged_event_summary) > 0L) {
        current <- audit$mice_logged_event_summary
        current$trial <- trial
        current$bootstrap_index <- bootstrap_index
        logged_event_rows[[length(logged_event_rows) + 1L]] <- current
      }
      if (
        !is.null(audit$mice_logged_event_route_summary) &&
          nrow(audit$mice_logged_event_route_summary) > 0L
      ) {
        current <- audit$mice_logged_event_route_summary
        current$trial <- trial
        current$bootstrap_index <- bootstrap_index
        event_route_rows[[length(event_route_rows) + 1L]] <- current
      }
      timing_rows[[length(timing_rows) + 1L]] <- data.frame(
        trial = trial,
        bootstrap_index = bootstrap_index,
        started_at = if (is.null(audit$started_at)) NA_character_ else audit$started_at,
        completed_at = if (is.null(audit$completed_at)) NA_character_ else audit$completed_at,
        elapsed_seconds = if (is.null(audit$elapsed_seconds)) NA_real_ else audit$elapsed_seconds,
        stringsAsFactors = FALSE
      )
    }
  }

  audit_index <- do.call(rbind, audit_index_rows)
  stopifnot(
    sum(audit_index$trial == "anti_mrsa") == 1000L,
    sum(audit_index$trial == "anti_psa") == 999L,
    nrow(audit_index) == 1999L,
    all(audit_index$result_specification_n == length(expected_specifications)),
    all(audit_index$all_logged_events_allowed),
    all(audit_index$mice_event_route_disallowed_n == 0L),
    all(audit_index$all_race_white_lindep_events_verified),
    all(audit_index$all_zero_rr_boundaries_verified)
  )

  empty_with_identity <- function(template) {
    template$trial <- character(0)
    template$bootstrap_index <- integer(0)
    template
  }
  race_retry <- if (length(race_retry_rows) == 0L) {
    empty_with_identity(empty_race_white_retry_audit_v2())
  } else {
    do.call(rbind, race_retry_rows)
  }
  existing_retry <- if (length(existing_retry_rows) == 0L) {
    empty_with_identity(empty_existing_logreg_retry_audit_v2())
  } else {
    do.call(rbind, existing_retry_rows)
  }
  logreg100_escalation <- if (
    length(logreg100_escalation_rows) == 0L
  ) {
    empty_with_identity(empty_logreg100_escalation_audit_v2())
  } else {
    do.call(rbind, logreg100_escalation_rows)
  }
  cbps_retry <- if (length(cbps_retry_rows) == 0L) {
    empty_cbps_retry_audit_v2()
  } else {
    do.call(rbind, cbps_retry_rows)
  }
  zero_rr_boundary <- if (length(zero_rr_boundary_rows) == 0L) {
    empty_zero_rr_boundary_summary_v2()
  } else {
    do.call(rbind, zero_rr_boundary_rows)
  }
  b0373_zero_rr <- zero_rr_boundary[
    zero_rr_boundary$trial == "anti_mrsa" &
      zero_rr_boundary$bootstrap_index == 373L,
    ,
    drop = FALSE
  ]
  stopifnot(
    nrow(b0373_zero_rr) == 1L,
    identical(b0373_zero_rr$specification, "complete_case"),
    identical(b0373_zero_rr$outcome, "death_30d"),
    isTRUE(b0373_zero_rr$boundary_verified)
  )
  race_white_lindep <- if (length(race_white_lindep_rows) == 0L) {
    empty_race_white_lindep_summary_v3()
  } else {
    do.call(rbind, race_white_lindep_rows)
  }
  b0584_race_lindep <- race_white_lindep[
    race_white_lindep$trial == "anti_mrsa" &
      race_white_lindep$bootstrap_index == 584L,
    ,
    drop = FALSE
  ]
  stopifnot(
    nrow(b0584_race_lindep) == 46L,
    all(b0584_race_lindep$strategy == race_white_lindep_strategy_v3),
    all(b0584_race_lindep$event_verified),
    sum(b0584_race_lindep$out ==
      "other_broad_spectrum_antibiotics.L") == 20L,
    sum(b0584_race_lindep$out ==
      "other_broad_spectrum_antibiotics.Q") == 26L
  )
  logged_events <- if (length(logged_event_rows) == 0L) {
    empty <- summarize_logged_events_v2(empty_mice_logged_events_v2())
    empty$trial <- character(0)
    empty$bootstrap_index <- integer(0)
    empty
  } else {
    do.call(rbind, logged_event_rows)
  }
  event_routes <- if (length(event_route_rows) == 0L) {
    empty <- empty_mice_event_route_summary_v4()
    empty$trial <- character(0)
    empty$bootstrap_index <- integer(0)
    empty
  } else {
    do.call(rbind, event_route_rows)
  }
  b0635_event_routes <- event_routes[
    event_routes$trial == "anti_mrsa" &
      event_routes$bootstrap_index == 635L,
    ,
    drop = FALSE
  ]
  stopifnot(
    is_exact_b0635_route_summary_v4(b0635_event_routes),
    sum(as.integer(b0635_event_routes$event_n)) == 5000L,
    sum(as.integer(b0635_event_routes$existing_allowed_n)) == 5000L,
    sum(as.integer(b0635_event_routes$race_lindep_allowed_n)) == 0L,
    sum(as.integer(b0635_event_routes$corrected_allowed_n)) == 5000L
  )
  timing <- do.call(rbind, timing_rows)

  atomic_write_csv(
    audit_index,
    file.path(sensitivity_dir, "ST-08_controlled_v2_repeat_audit_index.csv")
  )
  atomic_write_csv(
    race_retry,
    file.path(sensitivity_dir, "ST-08_controlled_v2_race_white_logreg_retries.csv")
  )
  atomic_write_csv(
    existing_retry,
    file.path(sensitivity_dir, "ST-08_controlled_v2_existing_logreg_retries.csv")
  )
  atomic_write_csv(
    logreg100_escalation,
    file.path(
      sensitivity_dir,
      "ST-08_controlled_v2_logreg100_escalations.csv"
    )
  )
  atomic_write_csv(
    cbps_retry,
    file.path(sensitivity_dir, "ST-08_controlled_v2_cbps_retries.csv")
  )
  atomic_write_csv(
    zero_rr_boundary,
    file.path(
      sensitivity_dir,
      "ST-08_controlled_v2_zero_rr_boundaries.csv"
    )
  )
  atomic_write_csv(
    race_white_lindep,
    file.path(
      sensitivity_dir,
      "ST-08_controlled_v2_race_white_lindep_events.csv"
    )
  )
  atomic_write_csv(
    logged_events,
    file.path(sensitivity_dir, "ST-08_controlled_v2_mice_logged_events.csv")
  )
  atomic_write_csv(
    event_routes,
    file.path(sensitivity_dir, "ST-08_controlled_v2_mice_event_routes.csv")
  )
  atomic_write_csv(
    timing,
    file.path(sensitivity_dir, "ST-08_controlled_v2_repeat_timing.csv")
  )

  # 复用10b已经审计的汇总后缀：type=7区间、原始样本点估计、ST-08 CSV/DOCX。
  controlled_v1_lines <- readLines(controlled_v1_path, warn = FALSE, encoding = "UTF-8")
  finalization_start <- which(
    trimws(controlled_v1_lines) == "# 合并通过单重复验证的正式结果。"
  )
  if (length(finalization_start) != 1L) {
    stop("10b最终汇总分界标记不唯一。")
  }
  eval(
    parse(text = controlled_v1_lines[finalization_start:length(controlled_v1_lines)]),
    envir = .GlobalEnv
  )

  interval_path <- file.path(sensitivity_dir, "ST-08_sensitivity_percentile_ci.csv")
  interval_table <- utils::read.csv(
    interval_path,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  stopifnot(
    all(c(
      "trial", "specification", "estimand_field", "bootstrap_B",
      "ci_lower", "bootstrap_median", "ci_upper", "point_estimate"
    ) %in% names(interval_table)),
    !anyNA(interval_table$point_estimate),
    all(is.finite(interval_table$point_estimate)),
    all(is.finite(interval_table$ci_lower)),
    all(is.finite(interval_table$ci_upper)),
    all(interval_table$ci_lower <= interval_table$ci_upper)
  )

  figure_data <- interval_table[
    interval_table$estimand_field == "rd_aki",
    ,
    drop = FALSE
  ]
  if (nrow(figure_data) != 2L * length(expected_specifications)) {
    stop("SF-09的AKI风险差行数不符合预期。")
  }
  figure_data$specification <- factor(
    figure_data$specification,
    levels = rev(expected_specifications)
  )
  figure_data$trial_label <- factor(
    figure_data$trial,
    levels = c("anti_mrsa", "anti_psa"),
    labels = c("Anti-MRSA cohort", "Anti-PSA cohort")
  )
  figure_data$point_percent <- 100 * figure_data$point_estimate
  figure_data$lower_percent <- 100 * figure_data$ci_lower
  figure_data$upper_percent <- 100 * figure_data$ci_upper

  sf09 <- ggplot2::ggplot(
    figure_data,
    ggplot2::aes(
      x = point_percent,
      y = specification
    )
  ) +
    ggplot2::geom_vline(
      xintercept = 0,
      linewidth = 0.35,
      linetype = "dashed",
      colour = "grey45"
    ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(xmin = lower_percent, xmax = upper_percent),
      width = 0.18,
      orientation = "y",
      linewidth = 0.45,
      colour = "#264653"
    ) +
    ggplot2::geom_point(size = 2.2, shape = 21, fill = "#2A9D8F", colour = "black") +
    ggplot2::facet_wrap(~trial_label, ncol = 1L, scales = "free_y") +
    ggplot2::labs(
      x = "Risk difference for 7-day AKI, percentage points (95% percentile CI)",
      y = NULL,
      caption = paste0(
        "Anti-MRSA: 1000 bootstrap replicates; anti-PSA: 999 successful ",
        "replicates (b=934 preserved as failed and not replaced)."
      )
    ) +
    ggplot2::theme_classic(base_size = 10) +
    ggplot2::theme(
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold"),
      axis.text.y = ggplot2::element_text(size = 8),
      plot.caption = ggplot2::element_text(hjust = 0, size = 8)
    )

  figure_base <- file.path(sensitivity_dir, "SF-09_sensitivity_aki_risk_difference")
  ggplot2::ggsave(
    paste0(figure_base, ".png"),
    sf09,
    width = 8.0,
    height = 8.5,
    units = "in",
    dpi = 300,
    bg = "white"
  )
  ggplot2::ggsave(
    paste0(figure_base, ".tiff"),
    sf09,
    width = 8.0,
    height = 8.5,
    units = "in",
    dpi = 600,
    compression = "lzw",
    bg = "white"
  )
  ggplot2::ggsave(
    paste0(figure_base, ".pdf"),
    sf09,
    width = 8.0,
    height = 8.5,
    units = "in",
    device = grDevices::cairo_pdf
  )

  final_audit <- data.frame(
    item = c(
      "status",
      "anti_mrsa_formal_repeat_n",
      "anti_psa_formal_repeat_n",
      "anti_psa_b0934_status",
      "anti_mrsa_b0031_original_failure_md5",
      "anti_psa_b0052_original_cbps_failure_md5",
      "anti_mrsa_b0140_original_cbps_failure_md5",
      "anti_psa_b0079_original_logreg100_failure_md5",
      "anti_psa_b0079_logreg100_diagnostic_md5",
      "anti_mrsa_b0373_original_zero_rr_failure_md5",
      "anti_mrsa_b0373_zero_rr_diagnostic_md5",
      "anti_mrsa_b0584_original_mice_event_failure_md5",
      "anti_mrsa_b0584_mice_event_diagnostic_md5",
      "anti_mrsa_b0584_candidate_preflight_md5",
      "anti_mrsa_b0635_original_mice_event_failure_md5",
      "anti_mrsa_b0635_event_routing_diagnostic_md5",
      "anti_mrsa_b0635_event_routing_diagnostic_script_md5",
      "anti_psa_b0934_failure_md5",
      "reviewed_finalizer_geom_errorbarh_failure_md5",
      "cbps_verified_retry_n",
      "logreg100_verified_escalation_n",
      "verified_zero_rr_boundary_n",
      "verified_race_white_lindep_event_n",
      "anti_mrsa_b0635_existing_v2_event_route_n",
      "anti_mrsa_b0635_race_white_lindep_v3_event_route_n",
      "locked_sensitivity_script_md5",
      "controlled_v1_script_md5",
      "controlled_v2_worker_md5",
      "controlled_v2_finalizer_md5",
      "percentile_interval",
      "pooling_authorized"
    ),
    value = c(
      "success",
      "1000",
      "999",
      "preserved_failed_not_replaced",
      unname(tools::md5sum(reviewed_b0031_failure_path_v2)),
      unname(tools::md5sum(
        reviewed_cbps_failure_registry_v2$path[
          reviewed_cbps_failure_registry_v2$trial == "anti_psa"
        ]
      )),
      unname(tools::md5sum(
        reviewed_cbps_failure_registry_v2$path[
          reviewed_cbps_failure_registry_v2$trial == "anti_mrsa"
        ]
      )),
      unname(tools::md5sum(reviewed_b0079_failure_path_v2)),
      unname(tools::md5sum(reviewed_b0079_diagnostic_path_v2)),
      unname(tools::md5sum(reviewed_b0373_failure_path_v2)),
      unname(tools::md5sum(reviewed_b0373_diagnostic_path_v2)),
      unname(tools::md5sum(reviewed_b0584_failure_path_v2)),
      unname(tools::md5sum(reviewed_b0584_diagnostic_path_v2)),
      unname(tools::md5sum(reviewed_b0584_preflight_path_v2)),
      unname(tools::md5sum(reviewed_b0635_failure_path_v2)),
      unname(tools::md5sum(reviewed_b0635_diagnostic_path_v2)),
      unname(tools::md5sum(reviewed_b0635_diagnostic_script_path_v2)),
      unname(tools::md5sum(preserved_failure_path)),
      unname(tools::md5sum(reviewed_finalizer_warning_failure_path_v2)),
      as.character(nrow(cbps_retry)),
      as.character(nrow(logreg100_escalation)),
      as.character(nrow(zero_rr_boundary)),
      as.character(nrow(race_white_lindep)),
      as.character(sum(as.integer(b0635_event_routes$existing_allowed_n))),
      as.character(sum(as.integer(b0635_event_routes$race_lindep_allowed_n))),
      unname(tools::md5sum(locked_script_path)),
      unname(tools::md5sum(controlled_v1_path)),
      unname(tools::md5sum(controlled_v2_worker_path)),
      unname(tools::md5sum(finalizer_path)),
      "two-sided 95% percentile CI; R stats::quantile(type=7)",
      "false"
    ),
    stringsAsFactors = FALSE
  )
  atomic_write_csv(
    final_audit,
    file.path(sensitivity_dir, "ST-08_controlled_v2_final_audit.csv")
  )

  required_outputs <- c(
    "ST-08_sensitivity_bootstrap_replicates.csv",
    "ST-08_sensitivity_point_estimates.csv",
    "ST-08_sensitivity_percentile_ci.csv",
    "ST-08_prespecified_sensitivity_results.csv",
    "ST-08_prespecified_sensitivity_results.docx",
    "ST-08_controlled_v2_repeat_audit_index.csv",
    "ST-08_controlled_v2_cbps_retries.csv",
    "ST-08_controlled_v2_logreg100_escalations.csv",
    "ST-08_controlled_v2_zero_rr_boundaries.csv",
    "ST-08_controlled_v2_race_white_lindep_events.csv",
    "ST-08_controlled_v2_mice_event_routes.csv",
    "ST-08_controlled_v2_final_audit.csv",
    "SF-09_sensitivity_aki_risk_difference.png",
    "SF-09_sensitivity_aki_risk_difference.tiff",
    "SF-09_sensitivity_aki_risk_difference.pdf"
  )
  required_paths <- file.path(sensitivity_dir, required_outputs)
  if (!all(file.exists(required_paths)) || any(file.info(required_paths)$size <= 0L)) {
    stop("ST-08/SF-09最终输出缺失或为空。")
  }

  message(
    "CONTROLLED_ST08_V2_FINAL_STATUS=success; anti_mrsa_B=1000; ",
    "anti_psa_B=999; anti_psa_b0934=preserved_failed; ",
    "reviewed_cbps_failures=2; verified_cbps_retries=", nrow(cbps_retry), "; ",
    "reviewed_logreg100_failures=1; verified_logreg100_escalations=",
    nrow(logreg100_escalation), "; reviewed_zero_rr_failures=1; ",
    "verified_zero_rr_boundaries=", nrow(zero_rr_boundary), "; ",
    "reviewed_race_white_lindep_failures=1; ",
    "verified_race_white_lindep_events=", nrow(race_white_lindep), "; ",
    "reviewed_event_routing_failures=1; ",
    "anti_mrsa_b0635_existing_v2_event_routes=",
    sum(as.integer(b0635_event_routes$existing_allowed_n)), "; ",
    "pooling_authorized=false"
  )
  invisible(TRUE)
}

finalization <- tryCatch(
  withCallingHandlers(
    finalize_st08_v2(),
    warning = function(condition) {
      stop("V2最终汇总出现未分类警告：", conditionMessage(condition), call. = FALSE)
    }
  ),
  error = function(condition) condition
)
if (inherits(finalization, "error")) {
  atomic_save_rds_controlled(
    list(
      status = "failed",
      stage = "controlled_v2_finalization",
      message = conditionMessage(finalization),
      controlled_v2_worker_md5 = unname(tools::md5sum(controlled_v2_worker_path)),
      controlled_v2_finalizer_md5 = unname(tools::md5sum(finalizer_path)),
      timestamp = format(Sys.time(), tz = "Asia/Shanghai", usetz = TRUE)
    ),
    finalizer_failure_path_v2
  )
  stop("受控V2最终汇总失败：", conditionMessage(finalization))
}
