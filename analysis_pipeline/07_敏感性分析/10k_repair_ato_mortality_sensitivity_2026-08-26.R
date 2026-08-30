# ==============================================================================
# ST-08 ATO 30-day mortality sensitivity-analysis repair
# 方法：复用锁定的 Bootstrap/MICE 检查点；对保留的死亡结局敏感性规格以
#       与主分析一致的 br.logit overlap weighting（ATO）重新估计。
# 代码来源：
#   1) 07_敏感性分析/10_run_sensitivity_analyses.R（锁定的规格和结局函数）
#   2) 05_Bootstrap/04_bootstrap_single_repeat.R（ATO br.logit 权重函数）
#   3) WeightIt 595d2812c2f4d56ff7c14aaf6fcde76c59aa8c8f；mice
#      61083e667fa41cc4bc655492591b100bf8a7e443
# 访问日期：2026-08-26
# 适配说明：不修改锁定脚本、原始数据、主分析或 CBPS-ATE 结果；仅在独立目录
#       以 ATO 重算可影响 30 天死亡估计的预设规格。AKI-only 定义和 ATE 截尾
#       规格不进入本次死亡结局敏感性分析。
# 验证命令：
#   D:/software/R4.4.2/R-4.4.2/bin/x64/Rscript.exe --vanilla
#   07_敏感性分析/10k_repair_ato_mortality_sensitivity_2026-08-26.R
# ==============================================================================

get_script_path_ato_repair <- function() {
  file_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_argument) != 1L) {
    stop("请通过 Rscript 运行本脚本，以便记录可追溯的脚本路径。")
  }
  normalizePath(sub("^--file=", "", file_argument), winslash = "/", mustWork = TRUE)
}

get_argument_ato_repair <- function(name, default = NULL) {
  prefix <- paste0("--", name, "=")
  hit <- commandArgs(trailingOnly = TRUE)
  hit <- hit[startsWith(hit, prefix)]
  if (length(hit) == 0L) return(default)
  if (length(hit) != 1L) stop("参数重复：", name)
  sub(prefix, "", hit, fixed = TRUE)
}

atomic_save_rds_ato_repair <- function(object, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary_path <- tempfile(pattern = paste0(basename(path), ".tmp_"), tmpdir = dirname(path))
  saveRDS(object, temporary_path, version = 3)
  if (!file.rename(temporary_path, path)) {
    stop("无法原子写入 RDS：", path)
  }
  invisible(path)
}

script_path <- get_script_path_ato_repair()
package_root <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = TRUE)
locked_script_path <- file.path(package_root, "07_敏感性分析", "10_run_sensitivity_analyses.R")
if (!file.exists(locked_script_path)) stop("缺少锁定敏感性分析脚本：", locked_script_path)

# 只载入锁定脚本中的配置和函数定义，不执行其原 CBPS-ATE 驱动程序。
locked_lines <- readLines(locked_script_path, warn = FALSE, encoding = "UTF-8")
driver_marker <- which(trimws(locked_lines) == "replicate_rows <- list()")
if (length(driver_marker) != 1L) stop("锁定脚本的执行分界标记不唯一。")
eval(parse(text = locked_lines[seq_len(driver_marker - 1L)]), envir = .GlobalEnv)

# 保存锁定函数，并仅在本脚本的局部重算路径中强制使用主分析 ATO 权重模型。
fit_bootstrap_weights_locked <- fit_bootstrap_weights
# 对既有函数中唯一且精确的 brglmFit 未收敛路径，将同一模型的单次 maxit
# 条件式重试上限由 200 提升至 2000；不改变公式、链接函数、估计量或随机流。
ato_brglm_retry_maxit <- 2000L
fit_bootstrap_weights <- function(data, formula, model) {
  if (!identical(model, "cbps_ate")) stop("修复脚本只接受锁定规格函数传入的 cbps_ate 占位符。")
  fit_bootstrap_weights_locked(data = data, formula = formula, model = "ato")
}

# NEW_MS03 中 race_white 的锁定 logreg 在极少数重复可出现唯一的未收敛警告。
# 仅对该精确情形，以相同增广数据和 quasi-binomial logit 将 maxit 从 25 提高至
# 2000 重试一次；其他警告或不收敛均停止。该适配复用受控 V2 的已审计规则。
mice.impute.racewhitelogreg <- function(y, ry, x, wy = NULL, ...) {
  if (is.null(wy)) wy <- !ry
  augmented <- mice:::augment(y, ry, x, wy)
  x <- cbind(1, as.matrix(augmented$x))
  y <- augmented$y
  ry <- augmented$ry
  wy <- augmented$wy
  weight <- augmented$w
  fit_warnings <- character(0)
  fit <- withCallingHandlers(
    stats::glm.fit(x = x[ry, , drop = FALSE], y = y[ry], family = stats::quasibinomial(link = "logit"), weights = weight[ry], control = stats::glm.control(maxit = 25)),
    warning = function(condition) {
      fit_warnings <<- c(fit_warnings, conditionMessage(condition))
      invokeRestart("muffleWarning")
    }
  )
  if (length(fit_warnings) == 1L && identical(fit_warnings, "glm.fit: algorithm did not converge") && !isTRUE(fit$converged)) {
    retry_warnings <- character(0)
    fit <- withCallingHandlers(
      stats::glm.fit(x = x[ry, , drop = FALSE], y = y[ry], family = stats::quasibinomial(link = "logit"), weights = weight[ry], control = stats::glm.control(maxit = 2000)),
      warning = function(condition) {
        retry_warnings <<- c(retry_warnings, conditionMessage(condition))
        invokeRestart("muffleWarning")
      }
    )
    if (length(retry_warnings) > 0L || !isTRUE(fit$converged) || any(!is.finite(stats::coef(fit)))) {
      stop("race_white logreg 条件式重试失败。")
    }
  } else if (length(fit_warnings) > 0L || !isTRUE(fit$converged) || any(!is.finite(stats::coef(fit)))) {
    stop("race_white logreg 出现未批准的异常。")
  }
  beta_star <- stats::coef(fit) + t(chol(mice:::sym(summary.glm(fit)$cov.unscaled))) %*% stats::rnorm(length(stats::coef(fit)))
  probability <- 1 / (1 + exp(-(x[wy, , drop = FALSE] %*% beta_star)))
  if (any(!is.finite(probability))) stop("race_white 产生非有限插补概率。")
  imputed <- stats::runif(nrow(probability)) <= probability
  imputed[imputed] <- 1
  if (is.factor(y)) imputed <- factor(imputed, c(0, 1), levels(y))
  imputed
}
rebuild_sensitivity_method_locked <- rebuild_sensitivity_method
rebuild_sensitivity_method <- function(trial, field_order) {
  method <- rebuild_sensitivity_method_locked(trial, field_order)
  if (!identical(unname(method["race_white"]), "logreg")) stop("NEW_MS03 的 race_white 锁定方法不是 logreg。")
  method["race_white"] <- "racewhitelogreg"
  method
}

# anti-PSA NEW_MS03 的 nephrotoxic_drug_exposure 使用项目锁定的 logreg100。
# 仅在 maxit=100 与 2000 均出现同一精确未收敛警告时，按受控 V2 规则执行
# 一次 maxit=20000 的最终重试；模型、数据、随机流和插补公式均不改变。
mice.impute.logreg100 <- function(y, ry, x, wy = NULL, ...) {
  if (is.null(wy)) wy <- !ry
  augmented <- mice:::augment(y, ry, x, wy)
  x <- cbind(1, as.matrix(augmented$x))
  y <- augmented$y
  ry <- augmented$ry
  wy <- augmented$wy
  weight <- augmented$w
  fit_once <- function(maxit) {
    warnings <- character(0)
    fit <- withCallingHandlers(
      stats::glm.fit(x = x[ry, , drop = FALSE], y = y[ry], family = stats::quasibinomial(link = "logit"), weights = weight[ry], control = stats::glm.control(maxit = maxit)),
      warning = function(condition) {
        warnings <<- c(warnings, conditionMessage(condition))
        invokeRestart("muffleWarning")
      }
    )
    list(fit = fit, warnings = warnings)
  }
  acceptable <- function(result) {
    length(result$warnings) == 0L && isTRUE(result$fit$converged) && all(is.finite(stats::coef(result$fit)))
  }
  exact_nonconvergence <- function(result) {
    length(result$warnings) == 1L && identical(result$warnings, "glm.fit: algorithm did not converge") && !isTRUE(result$fit$converged)
  }
  result <- fit_once(100L)
  if (exact_nonconvergence(result)) {
    result <- fit_once(2000L)
    if (exact_nonconvergence(result)) result <- fit_once(20000L)
  }
  if (!acceptable(result)) stop("logreg100 受控重试失败。")
  fit <- result$fit
  covariance <- summary.glm(fit)$cov.unscaled
  if (any(!is.finite(covariance)) || inherits(try(chol(mice:::sym(covariance)), silent = TRUE), "try-error")) {
    stop("logreg100 的后验协方差矩阵未通过有限性或 Cholesky 验证。")
  }
  beta_star <- stats::coef(fit) + t(chol(mice:::sym(covariance))) %*% stats::rnorm(length(stats::coef(fit)))
  probability <- stats::plogis(drop(x[wy, , drop = FALSE] %*% beta_star))
  if (any(!is.finite(probability))) stop("logreg100 产生非有限插补概率。")
  imputed <- stats::runif(length(probability)) <= probability
  imputed[imputed] <- 1
  if (is.factor(y)) imputed <- factor(imputed, c(0, 1), levels(y))
  imputed
}

# 受控 V2 对 NEW_MS03 的 logged events 采用逐事件、逐结构审计；本修复脚本
# 仅复用其中的允许条件。任何未在下列条件中明确定义的事件都会中止运行。
is_allowed_new_ms03_event_ato_repair <- function(events, work_data, method) {
  if (nrow(events) == 0L) return(logical(0))
  ordered_contrast <- (
    events$dep == "urine_output_rate_24h_ml_kg_h" &
      events$meth == "pmm" &
      grepl("^other_broad_spectrum_antibiotics\\.(L|Q)(, other_broad_spectrum_antibiotics\\.(L|Q))*$", events$out)
  )
  zero_count_death <- vapply(seq_len(nrow(events)), function(i) {
    dep <- events$dep[i]
    dep %in% names(work_data) && dep %in% names(method) &&
      events$out[i] == "aki_competing_statedeath_before_aki" &&
      events$meth[i] == unname(method[dep]) &&
      sum(!is.na(work_data[[dep]]) & as.character(work_data$aki_competing_state) == "death_before_aki") == 0L
  }, logical(1))
  race_lindep <- (
    events$dep == "race_white" &
      events$meth == "racewhitelogreg" &
      grepl("^other_broad_spectrum_antibiotics\\.(L|Q)(, other_broad_spectrum_antibiotics\\.(L|Q))*$", events$out) &
      identical(levels(work_data$race_white), c("0", "1")) &
      identical(levels(work_data$other_broad_spectrum_antibiotics), c("0", "1", "2")) &
      "race_white" %in% rownames(predictor_matrix_ato_repair) &
      "other_broad_spectrum_antibiotics" %in% colnames(predictor_matrix_ato_repair) &
      identical(as.numeric(predictor_matrix_ato_repair["race_white", "other_broad_spectrum_antibiotics"]), 1)
  )
  ordered_contrast | zero_count_death | race_lindep
}

run_new_ms03_mice_ato_repair <- function(work_data, method, predictor_matrix, visit_sequence, seed, trial, bootstrap_index) {
  # anti-PSA 的 source_control_status 已有项目锁定的条件式 logreg 重试适配；
  # 仅将该变量由锁定的 logreg 路由至该适配方法。
  if (identical(trial, "anti_psa")) {
    if (!identical(unname(method["source_control_status"]), "logreg")) {
      stop("NEW_MS03 anti-PSA 的 source_control_status 锁定方法不是 logreg。")
    }
    method["source_control_status"] <- "sourcecontrollogreg"
  }
  predictor_matrix_ato_repair <<- predictor_matrix
  captured_warnings <- character(0)
  set.seed(as.integer(seed))
  imputation <- withCallingHandlers(
    mice::mice(
      work_data, m = mice_m, maxit = mice_maxit, method = method,
      predictorMatrix = predictor_matrix, visitSequence = visit_sequence,
      donors = mice_donors, seed = seed, printFlag = FALSE
    ),
    warning = function(condition) {
      message <- conditionMessage(condition)
      if (grepl("^Number of logged events: [0-9]+$", message)) {
        captured_warnings <<- c(captured_warnings, message)
        invokeRestart("muffleWarning")
      } else {
        # 仅在未分类警告时记录当前 MICE 调用栈中的目标变量；不改变插补或估计。
        active_frames <- sys.frames()
        active_variables <- unique(unlist(lapply(active_frames, function(frame) {
          if (exists("j", envir = frame, inherits = FALSE)) {
            value <- get("j", envir = frame, inherits = FALSE)
            return(as.character(value))
          }
          character(0)
        }), use.names = FALSE))
        stop(
          "NEW_MS03 MICE 出现未分类警告：", message,
          "; MICE_active_j=", paste(active_variables, collapse = ","),
          call. = FALSE
        )
      }
    }
  )
  events <- imputation$loggedEvents
  if (is.null(events)) events <- data.frame(it = integer(0), im = integer(0), dep = character(0), meth = character(0), out = character(0))
  expected_warning <- if (nrow(events) == 0L) character(0) else paste0("Number of logged events: ", nrow(events))
  allowed <- is_allowed_new_ms03_event_ato_repair(events, work_data, method)
  if (!identical(captured_warnings, expected_warning) || (nrow(events) > 0L && !all(allowed))) {
    stop("NEW_MS03 MICE 含未批准 logged event；事件数=", nrow(events))
  }
  attr(imputation, "ato_repair_new_ms03_event_audit") <- list(
    trial = trial, bootstrap_index = as.integer(bootstrap_index), events = events, allowed = allowed
  )
  imputation$loggedEvents <- NULL
  imputation
}

# 使用锁定函数的全部结局和诊断逻辑，仅替换 NEW_MS03 的 MICE 调用及其全停规则。
estimate_start <- which(trimws(locked_lines) == "estimate_specification <- function(")
estimate_end_marker <- which(grepl("^recalculated_specifications <- c\\(", trimws(locked_lines)))
if (length(estimate_start) != 1L || length(estimate_end_marker) != 1L || estimate_start >= estimate_end_marker) stop("无法定位锁定 estimate_specification 函数。")
estimate_lines <- locked_lines[estimate_start:(estimate_end_marker - 1L)]
mice_call_start <- which(trimws(estimate_lines) == "analysis_imputation <- mice::mice(")
logged_stop_start <- which(trimws(estimate_lines) == "if (!is.null(analysis_imputation$loggedEvents)) {")
if (length(mice_call_start) != 1L || length(logged_stop_start) != 1L || logged_stop_start <= mice_call_start) stop("无法定位锁定 NEW_MS03 MICE 块。")
logged_stop_end <- logged_stop_start + 2L
replacement_mice_block <- c(
  "    analysis_imputation <- run_new_ms03_mice_ato_repair(",
  "      work_data = work_data, method = method, predictor_matrix = predictor_matrix,",
  "      visit_sequence = visit_sequence, seed = sensitivity_seed, trial = trial,",
  "      bootstrap_index = if (!is.null(checkpoint$bootstrap_index)) as.integer(checkpoint$bootstrap_index) else 0L",
  "    )"
)
estimate_lines <- c(
  estimate_lines[seq_len(mice_call_start - 1L)], replacement_mice_block,
  estimate_lines[(logged_stop_end + 1L):length(estimate_lines)]
)
eval(parse(text = estimate_lines), envir = .GlobalEnv)

formal_bootstrap_indices <- list(
  anti_mrsa = seq_len(1000L),
  anti_psa = c(seq_len(933L), 935L:1000L)
)
expected_replicate_n <- c(anti_mrsa = 1000L, anti_psa = 999L)
mortality_specifications <- c(
  "death_calendar_boundary",
  "NEW_MS03",
  "ABX_TIME_proxy",
  "MICRO_STATE_proxy",
  "complete_case"
)
requested_specification <- get_argument_ato_repair("specification", NULL)
if (!is.null(requested_specification)) {
  if (!requested_specification %in% mortality_specifications) {
    stop("未知或不适用于死亡结局的 specification：", requested_specification)
  }
  mortality_specifications <- requested_specification
}
effect_fields <- c(
  "risk_aki_0", "risk_aki_1", "rd_aki", "rr_aki",
  "risk_death30_0", "risk_death30_1", "rd_death30", "rr_death30"
)
diagnostic_fields <- c(
  "n_analysis", "truncated_fraction", "ess_continuation_min", "ess_deescalation_min"
)

run_label <- get_argument_ato_repair("run_label", "full")
if (!grepl("^[A-Za-z0-9_-]+$", run_label)) stop("run_label 只能包含英文字母、数字、下划线或连字符。")
repair_dir <- file.path(output_root, "06_敏感性分析", "ato_mortality_repair_20260826", run_label)
repeat_dir <- file.path(repair_dir, "repeat_checkpoints")
failure_dir <- file.path(repair_dir, "failures")
dir.create(repeat_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(failure_dir, recursive = TRUE, showWarnings = FALSE)

validate_ato_repair_repeat <- function(result, trial, bootstrap_index) {
  required_columns <- c("trial", "bootstrap_index", "specification", "estimand", "weight_model", effect_fields, diagnostic_fields)
  stopifnot(
    is.data.frame(result),
    nrow(result) == length(mortality_specifications),
    all(required_columns %in% names(result)),
    all(result$trial == trial),
    all(as.integer(result$bootstrap_index) == bootstrap_index),
    setequal(result$specification, mortality_specifications),
    all(result$estimand == "ATO"),
    all(result$weight_model == "overlap_weight_brlogit"),
    !anyDuplicated(result$specification)
  )
  numeric_fields <- c(effect_fields, diagnostic_fields)
  if (!all(vapply(result[numeric_fields], function(x) is.numeric(x) && all(is.finite(x)), logical(1)))) {
    stop("ATO 重复结果包含非有限数值。")
  }
  if (any(result$risk_death30_0 < 0 | result$risk_death30_0 > 1 | result$risk_death30_1 < 0 | result$risk_death30_1 > 1)) {
    stop("ATO 重复结果包含越界的 30 天死亡风险。")
  }
  if (!all.equal(result$rd_death30, result$risk_death30_1 - result$risk_death30_0, tolerance = 1e-12)) {
    stop("ATO 重复结果的 30 天死亡 RD 不满足风险差恒等式。")
  }
  invisible(TRUE)
}

calculate_ato_repair_repeat <- function(checkpoint, raw_data, trial, bootstrap_index) {
  rows <- lapply(mortality_specifications, function(specification) {
    estimate <- estimate_specification(checkpoint, raw_data, trial, specification)
    data.frame(
      trial = trial,
      bootstrap_index = as.integer(bootstrap_index),
      specification = specification,
      estimand = "ATO",
      weight_model = "overlap_weight_brlogit",
      as.list(estimate),
      stringsAsFactors = FALSE
    )
  })
  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  validate_ato_repair_repeat(result, trial, bootstrap_index)
  result
}

selected_trial <- get_argument_ato_repair("trial", NULL)
if (!is.null(selected_trial) && !selected_trial %in% names(formal_bootstrap_indices)) {
  stop("未知 trial：", selected_trial)
}
selected_trials <- if (is.null(selected_trial)) names(formal_bootstrap_indices) else selected_trial
b_start <- as.integer(get_argument_ato_repair("b_start", NA_character_))
b_end <- as.integer(get_argument_ato_repair("b_end", NA_character_))
if (xor(is.na(b_start), is.na(b_end))) stop("b_start 和 b_end 必须同时提供。")

for (trial in selected_trials) {
  raw_data <- utils::read.csv(file.path(data_dir, paste0("analysis_dataset_", trial, ".csv")), stringsAsFactors = FALSE, check.names = FALSE)
  raw_data <- augment_creatinine_only_outcome(raw_data, trial)
  current_indices <- formal_bootstrap_indices[[trial]]
  if (!is.na(b_start)) current_indices <- current_indices[current_indices >= b_start & current_indices <= b_end]
  if (length(current_indices) == 0L) stop("所选 Bootstrap 编号不属于正式成功重复集。")
  for (position in seq_along(current_indices)) {
    bootstrap_index <- current_indices[position]
    repeat_path <- file.path(repeat_dir, sprintf("%s_b%04d_ato_mortality.rds", trial, bootstrap_index))
    failure_path <- file.path(failure_dir, sprintf("%s_b%04d_ato_mortality_failure.rds", trial, bootstrap_index))
    if (file.exists(failure_path)) stop("存在待处理失败证据：", failure_path)
    if (file.exists(repeat_path)) {
      validate_ato_repair_repeat(readRDS(repeat_path), trial, bootstrap_index)
    } else {
      checkpoint_path <- file.path(checkpoint_dir, sprintf("%s_b%04d_checkpoint.rds", trial, bootstrap_index))
      if (!file.exists(checkpoint_path)) stop("缺少 Bootstrap checkpoint：", checkpoint_path)
      checkpoint <- readRDS(checkpoint_path)
      stopifnot(identical(checkpoint$status, "success"), identical(checkpoint$trial, trial), identical(as.integer(checkpoint$bootstrap_index), bootstrap_index))
      calculation <- tryCatch(
        withCallingHandlers(
          calculate_ato_repair_repeat(checkpoint, raw_data, trial, bootstrap_index),
          warning = function(warning_condition) stop("未分类警告：", conditionMessage(warning_condition), call. = FALSE)
        ),
        error = function(error_condition) error_condition
      )
      if (inherits(calculation, "error")) {
        atomic_save_rds_ato_repair(list(
          status = "failed", trial = trial, bootstrap_index = bootstrap_index,
          locked_script_md5 = unname(tools::md5sum(locked_script_path)),
          repair_script_md5 = unname(tools::md5sum(script_path)),
          message = conditionMessage(calculation),
          timestamp = format(Sys.time(), tz = "Asia/Shanghai", usetz = TRUE)
        ), failure_path)
        stop("ATO 死亡敏感性分析失败：", trial, " b=", bootstrap_index, "；", conditionMessage(calculation))
      }
      atomic_save_rds_ato_repair(calculation, repeat_path)
    }
    if (position %% 10L == 0L || position == length(current_indices)) {
      message("ATO_MORTALITY_REPAIR_PROGRESS trial=", trial, " completed=", position, "/", length(current_indices), " last_b=", bootstrap_index)
    }
  }
}

# 仅在两个试验的全部正式重复均完成后汇总，避免以不完整 Bootstrap 分布生成结果。
all_repeat_paths <- unlist(lapply(names(formal_bootstrap_indices), function(trial) {
  file.path(repeat_dir, sprintf("%s_b%04d_ato_mortality.rds", trial, formal_bootstrap_indices[[trial]]))
}), use.names = FALSE)
if (!all(file.exists(all_repeat_paths))) {
  message("ATO_MORTALITY_REPAIR_STATUS=partial; 完成后以相同命令继续运行。")
  quit(save = "no", status = 0L)
}

replicates <- do.call(rbind, lapply(all_repeat_paths, readRDS))
for (trial in names(formal_bootstrap_indices)) {
  current <- replicates[replicates$trial == trial, , drop = FALSE]
  stopifnot(
    length(unique(current$bootstrap_index)) == expected_replicate_n[[trial]],
    nrow(current) == expected_replicate_n[[trial]] * length(mortality_specifications)
  )
}
atomic_write_csv(replicates, file.path(repair_dir, "ST-08_ato_mortality_sensitivity_bootstrap_replicates.csv"))

# 点估计使用与 Bootstrap 相同的 ATO 规格函数和原始样本 MICE 对象。
point_rows <- list()
death_boundary_audit <- list()
for (trial in names(formal_bootstrap_indices)) {
  raw_data <- utils::read.csv(file.path(data_dir, paste0("analysis_dataset_", trial, ".csv")), stringsAsFactors = FALSE, check.names = FALSE)
  raw_data <- augment_creatinine_only_outcome(raw_data, trial)
  primary_mice_path <- file.path(output_root, "02_MICE", paste0("mids_", trial, "_m20_maxit50.rds"))
  primary_mice <- readRDS(primary_mice_path)
  checkpoint <- list(sampled_row = seq_len(nrow(raw_data)), imputation = primary_mice$imputation, mice_seed = master_seed)
  for (specification in mortality_specifications) {
    estimate <- estimate_specification(checkpoint, raw_data, trial, specification)
    point_rows[[length(point_rows) + 1L]] <- data.frame(
      trial = trial, specification = specification, estimand = "ATO", weight_model = "overlap_weight_brlogit",
      as.list(estimate), stringsAsFactors = FALSE
    )
  }
  death_main <- as.integer(as.character(raw_data$death_30d))
  death_boundary <- as.integer(as.character(raw_data$death_30d_boundary_sensitivity))
  death_boundary_audit[[length(death_boundary_audit) + 1L]] <- data.frame(
    trial = trial,
    n_analysis = nrow(raw_data),
    n_death_reclassified = sum(death_main != death_boundary),
    n_changed_from_0_to_1 = sum(death_main == 0L & death_boundary == 1L),
    n_changed_from_1_to_0 = sum(death_main == 1L & death_boundary == 0L),
    stringsAsFactors = FALSE
  )
}
points <- do.call(rbind, point_rows)
death_boundary_audit <- do.call(rbind, death_boundary_audit)
atomic_write_csv(points, file.path(repair_dir, "ST-08_ato_mortality_sensitivity_point_estimates.csv"))
atomic_write_csv(death_boundary_audit, file.path(repair_dir, "ST-08_ato_mortality_death_boundary_audit.csv"))

summary_rows <- list()
for (trial in names(formal_bootstrap_indices)) {
  for (specification in mortality_specifications) {
    current <- replicates[replicates$trial == trial & replicates$specification == specification, "rd_death30"]
    point <- points[points$trial == trial & points$specification == specification, "rd_death30"]
    stopifnot(length(current) == expected_replicate_n[[trial]], length(point) == 1L)
    ci <- stats::quantile(current, probs = c(0.025, 0.5, 0.975), type = 7, names = FALSE)
    summary_rows[[length(summary_rows) + 1L]] <- data.frame(
      trial = trial, specification = specification, estimand_field = "rd_death30",
      estimand = "ATO", weight_model = "overlap_weight_brlogit", bootstrap_B = length(current),
      ci_lower = ci[1], bootstrap_median = ci[2], ci_upper = ci[3], point_estimate = point,
      estimate_95ci = sprintf("%.1f (%.1f to %.1f)", point * 100, ci[1] * 100, ci[3] * 100),
      stringsAsFactors = FALSE
    )
  }
}
summary_table <- do.call(rbind, summary_rows)
atomic_write_csv(summary_table, file.path(repair_dir, "ST-08_ato_mortality_sensitivity_results.csv"))

provenance <- data.frame(
  item = c("repair_script", "repair_script_md5", "locked_sensitivity_script", "locked_sensitivity_script_md5", "estimand", "weight_model", "anti_mrsa_bootstrap_B", "anti_psa_bootstrap_B", "excluded_specifications"),
  value = c(normalizePath(script_path, winslash = "/"), unname(tools::md5sum(script_path)), normalizePath(locked_script_path, winslash = "/"), unname(tools::md5sum(locked_script_path)), "ATO", "WeightIt glm / br.logit / overlap weights", "1000", "999", "aki_cross_t0_uo_boundary; creatinine_only_phenotype; weight_truncation_1_99; weight_truncation_5_95"),
  stringsAsFactors = FALSE
)
atomic_write_csv(provenance, file.path(repair_dir, "ST-08_ato_mortality_repair_provenance.csv"))
message("ATO_MORTALITY_REPAIR_STATUS=success; output=", repair_dir)
