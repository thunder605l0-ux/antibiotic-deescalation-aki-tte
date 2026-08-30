# ==============================================================================
# Analysis: 单次Boot-MI引擎与br.logit ATO适配
# Date: 2026-07-29
# Random seed: 由试验代码和Bootstrap重复号确定性派生
# R: 4.4.2
#
# Method:
# Purpose: 加载锁定Boot-MI核心，并仅把ATO倾向评分实现统一为br.logit。
# Source type: local_knowledge_base / GitHub / project_locked_code
# Source name or URL:
#   D:/source/obsidian/method_wiki/causal_inference/g_methods_ipw_standardization.md
#   https://github.com/ngreifer/WeightIt
#   https://cran.r-project.org/package=brglm2
# Citation or commit/version:
#   WeightIt commit 595d2812c2f4d56ff7c14aaf6fcde76c59aa8c8f
#   brglm2 1.0.1
#   Li F, Morgan KL, Zaslavsky AM. JASA. 2018;113:390-400.
#   D-4.7I-ATO-BRLOGIT-SENSITIVITY-PROTOCOL-V2.0
# Access date: 2026-07-29
# Adaptation notes:
#   抽样、MICE、CBPS-ATE、Aalen–Johansen、死亡风险及MI合并均继承锁定核心；
#   ATO仅把普通logit替换为偏倚校正br.logit，并核对手工重叠权重公式。
# Verification command:
#   parse(file="05_Bootstrap/04_bootstrap_single_repeat.R")
# ==============================================================================

# 读取source()调用记录中的当前文件路径；直接运行时该值可能为空。
# 使用用户锁定的绝对路径加载统一配置，不依赖当前工作目录或调用方式。
source(
  file.path(Sys.getenv("TTE_RUNTIME_ROOT", unset = "."), "02_配置与变量字典", "analysis_config.R"),
  local = globalenv(),
  encoding = "UTF-8"
)
# 将工作目录固定到运行包根目录，保证交互式与命令行执行一致。
setwd(package_root)

# 核对单次Boot-MI需要的全部R包。
require_packages(c(
  "mice",
  "WeightIt",
  "brglm2",
  "cobalt",
  "survival"
))

# 构造锁定Boot-MI核心函数文件路径。
core_path <- normalizePath(
  file.path(package_root, "05_Bootstrap", "04_bootstrap_core_locked.R"),
  winslash = "/",
  mustWork = TRUE
)
# 加载经逐行注释的锁定核心函数。
source(core_path, local = globalenv(), encoding = "UTF-8")
source(
  file.path(config_dir, "mice_imputation_extensions.R"),
  local = globalenv(),
  encoding = "UTF-8"
)

# 用统一配置覆盖锁定核心中的同名常量，避免两个位置产生参数漂移。
bootstrap_master_seed <- master_seed
# 固定单次重复的插补数据集数量。
bootstrap_m <- mice_m
# 固定单次重复的MICE迭代次数。
bootstrap_maxit <- mice_maxit
# 固定预测均值匹配供体数。
bootstrap_donors <- mice_donors
# 固定AKI主要随访时点。
bootstrap_horizon_aki <- aki_horizon_days
# 固定数值交叉核验容差。
bootstrap_tolerance <- numeric_tolerance

# br.logit仍使用相同的偏倚校正估计方法；仅当默认100次迭代明确未收敛时，
# 允许在相同数据、公式和起始规则下把优化迭代上限提高到200次后重试一次。
ato_brglm_retry_maxit <- 200L

# 重新定义权重拟合函数，以锁定br.logit ATO并保留CBPS-ATE原实现。
fit_bootstrap_weights <- function(data, formula, model) {
  # 把暴露因子安全转换为0/1整数。
  treatment <- as.integer(as.character(data$deescalation))
  # 检查当前完成数据是否同时保留两个策略组。
  if (
    anyNA(treatment) ||
      !setequal(unique(treatment), c(0L, 1L))
  ) {
    # 缺少任一策略组时停止，避免不可识别的权重模型。
    stop("Bootstrap样本或完成数据失去一个策略组。")
  }

  # 初始化警告收集器，任何模型警告都不允许静默通过。
  captured_warnings <- character(0)
  ato_brglm_retry_used <- FALSE
  ato_brglm_initial_warning_n <- 0L
  # 主模型分支拟合线性just-identified CBPS-ATE。
  if (model == "cbps_ate") {
    # 在调用WeightIt时捕获完整警告正文。
    fit <- withCallingHandlers(
      # 使用锁定WeightIt实现拟合未稳定化CBPS-ATE权重。
      WeightIt::weightit(
        formula = formula,
        data = data,
        method = "cbps",
        estimand = "ATE",
        stabilize = FALSE,
        # 固定恰好识别CBPS；over不控制PS重叠或positivity检查。
        over = FALSE,
        include.obj = TRUE
      ),
      # 定义WeightIt警告处理函数。
      warning = function(condition) {
        # 把本次警告正文追加到字符向量。
        captured_warnings <<- c(
          captured_warnings,
          conditionMessage(condition)
        )
        # 暂时抑制控制台重复打印，随后由硬门槛统一停止。
        invokeRestart("muffleWarning")
      }
    )
    # 提取每位患者的倾向评分。
    propensity_score <- as.numeric(fit$ps)
    # 提取WeightIt返回的未稳定化ATE权重。
    raw_weight <- as.numeric(fit$weights)
    # 计算当前试验和插补数据中的边际降阶梯概率。
    marginal_treated <- mean(treatment == 1L)
    # 按策略组边际概率显式构造稳定化ATE权重。
    weight <- ifelse(
      treatment == 1L,
      raw_weight * marginal_treated,
      raw_weight * (1 - marginal_treated)
    )
    # 要求CBPS优化器明确返回收敛码0。
    converged <- identical(extract_cbps_convergence_boot(fit), 0L)
    # 主CBPS不使用ATO手工权重误差，因此记为缺失。
    manual_weight_error <- NA_real_
    # CBPS分支不依赖brglm2系数有限性，固定记为TRUE。
    coefficients_finite <- TRUE
  } else if (model == "ato") {
    # ATO敏感性分支使用偏倚校正二项logistic倾向评分。
    fit <- withCallingHandlers(
      # 调用WeightIt的br.logit接口估计ATO重叠权重。
      WeightIt::weightit(
        formula = formula,
        data = data,
        method = "glm",
        estimand = "ATO",
        link = "br.logit",
        stabilize = FALSE,
        include.obj = TRUE
      ),
      # 定义ATO拟合警告处理函数。
      warning = function(condition) {
        # 保存警告正文，避免会话结束时只留下警告数量。
        captured_warnings <<- c(
          captured_warnings,
          conditionMessage(condition)
        )
        # 暂时抑制重复打印，随后由统一门槛处理。
        invokeRestart("muffleWarning")
      }
    )
    exact_nonconvergence_warning <- (
      length(captured_warnings) > 0L &&
        all(startsWith(
          captured_warnings,
          "(from `glm()`): brglmFit: algorithm did not converge."
        )) &&
        !isTRUE(fit$obj$converged)
    )
    if (exact_nonconvergence_warning) {
      ato_brglm_initial_warning_n <- length(captured_warnings)
      captured_warnings <- character(0)
      ato_brglm_retry_used <- TRUE
      fit <- withCallingHandlers(
        WeightIt::weightit(
          formula = formula,
          data = data,
          method = "glm",
          estimand = "ATO",
          link = "br.logit",
          stabilize = FALSE,
          include.obj = TRUE,
          control = list(maxit = ato_brglm_retry_maxit)
        ),
        warning = function(condition) {
          captured_warnings <<- c(
            captured_warnings,
            conditionMessage(condition)
          )
          invokeRestart("muffleWarning")
        }
      )
    }
    # 提取偏倚校正模型给出的倾向评分。
    propensity_score <- as.numeric(fit$ps)
    # 提取WeightIt返回的ATO重叠权重。
    weight <- as.numeric(fit$weights)
    # 按ATO定义手工重建每位患者的理论重叠权重。
    manual_weight <- ifelse(
      treatment == 1L,
      1 - propensity_score,
      propensity_score
    )
    # 计算软件权重与手工公式之间的最大绝对误差。
    manual_weight_error <- max(abs(weight - manual_weight))
    # 检查偏倚校正logistic模型全部系数是否为有限值。
    coefficients_finite <- all(is.finite(stats::coef(fit$obj)))
    # 同时要求模型收敛标志和系数有限性通过。
    converged <- isTRUE(fit$obj$converged) && coefficients_finite
  } else {
    # 拒绝未在分析计划中登记的权重模型名称。
    stop("未知权重模型：", model)
  }

  # 汇总检查警告、收敛、PS范围、权重范围和ATO公式一致性。
  if (
    length(captured_warnings) > 0L ||
      !converged ||
      length(propensity_score) != nrow(data) ||
      any(!is.finite(propensity_score)) ||
      any(propensity_score <= 0 | propensity_score >= 1) ||
      length(weight) != nrow(data) ||
      any(!is.finite(weight)) ||
      any(weight <= 0) ||
      (
        model == "ato" &&
          manual_weight_error > sqrt(.Machine$double.eps)
      )
  ) {
    # 任一硬门槛失败时停止并保留所有模型警告正文。
    stop(
      model,
      "拟合或权重QA失败；warnings=",
      paste(captured_warnings, collapse = " | ")
    )
  }

  # 使用cobalt计算未加权与加权后的协变量平衡。
  balance <- cobalt::bal.tab(
    formula,
    data = data,
    weights = weight,
    estimand = if (model == "cbps_ate") "ATE" else "ATO",
    s.d.denom = "pooled",
    binary = "std",
    continuous = "std",
    un = TRUE
  )$Balance
  # 去除倾向评分距离本身，仅保留预设操作变量。
  balance <- balance[balance$Type != "Distance", , drop = FALSE]
  # 计算每个操作变量的加权绝对标准化差异。
  absolute_smd <- abs(balance$Diff.Adj)
  # 定位加权后绝对SMD最大的变量。
  maximum_smd_index <- which.max(absolute_smd)
  # 分策略计算有效样本量。
  ess <- vapply(
    c(0L, 1L),
    function(strategy) {
      # 使用Kish公式计算当前策略组ESS。
      effective_sample_size_boot(weight[treatment == strategy])
    },
    numeric(1)
  )

  # 返回权重、倾向评分和平衡诊断，不在此函数读取结局。
  list(
    weight = weight,
    propensity_score = propensity_score,
    max_absolute_smd = absolute_smd[maximum_smd_index],
    max_absolute_smd_term = rownames(balance)[maximum_smd_index],
    ess_continuation = unname(ess[1]),
    ess_deescalation = unname(ess[2]),
    max_weight = max(weight),
    warning_n = length(captured_warnings),
    ato_link = if (model == "ato") "br.logit" else NA_character_,
    coefficients_finite = coefficients_finite,
    manual_weight_max_absolute_error = manual_weight_error,
    ato_brglm_retry_used = ato_brglm_retry_used,
    ato_brglm_initial_warning_n = ato_brglm_initial_warning_n,
    ato_brglm_retry_maxit = if (ato_brglm_retry_used) {
      ato_brglm_retry_maxit
    } else {
      NA_integer_
    }
  )
}

# 保存锁定核心中仅计算第7天风险的结局函数，供扩展函数复用。
estimate_bootstrap_outcomes_horizon_only <- estimate_bootstrap_outcomes
# 保存锁定核心的单次Bootstrap主函数，供曲线扩展包装器复用。
run_single_bootstrap_repeat_horizon_only <- run_single_bootstrap_repeat

# 定义在固定日网格读取加权Aalen–Johansen AKI累积发生概率的函数。
weighted_aj_curve_boot <- function(time, event_factor, weight) {
  # 先使用锁定输入QA函数检查时间、事件和权重。
  validate_aki_boot(time, event_factor, weight)
  # 使用与第7天点估计相同参数拟合加权Aalen–Johansen。
  fit <- survival::survfit(
    survival::Surv(time, event_factor) ~ 1,
    data = data.frame(time, event_factor, weight),
    weights = weight,
    stype = 1,
    ctype = 1,
    se.fit = FALSE,
    time0 = TRUE,
    timefix = FALSE
  )
  # 固定图形时间网格为T0后0至7天，每1天一个点。
  time_grid <- seq(0, bootstrap_horizon_aki, by = 1)
  # 在固定时间网格读取状态概率，并允许末次状态向后保持。
  fit_summary <- summary(
    fit,
    times = time_grid,
    extend = TRUE
  )
  # 把状态概率强制转换为矩阵以稳定单时间点边界行为。
  probability_matrix <- as.matrix(fit_summary$pstate)
  # 取得AKI状态列的位置。
  aki_column <- match("aki", colnames(probability_matrix))
  # 找不到AKI状态列时停止。
  if (is.na(aki_column)) {
    # 报告曲线对象结构错误。
    stop("Aalen–Johansen曲线对象缺少AKI状态。")
  }
  # 提取每个时间网格的AKI累积发生概率。
  aki_risk <- probability_matrix[, aki_column]
  # 要求风险有限、位于0至1且不下降。
  if (
    any(!is.finite(aki_risk)) ||
      any(aki_risk < -bootstrap_tolerance) ||
      any(aki_risk > 1 + bootstrap_tolerance) ||
      any(diff(aki_risk) < -bootstrap_tolerance)
  ) {
    # 曲线违反概率或单调性时停止。
    stop("AKI累积发生曲线QA失败。")
  }
  # 返回固定时间网格及AKI风险。
  data.frame(
    time_days = time_grid,
    aki_risk = as.numeric(aki_risk),
    stringsAsFactors = FALSE
  )
}

# 初始化曲线收集环境；每个Bootstrap重复开始时由包装器重置。
.curve_collector <- new.env(parent = emptyenv())
# 初始化曲线调用计数。
.curve_collector$call_index <- 0L
# 初始化逐插补、逐模型曲线列表。
.curve_collector$rows <- list()

# 扩展结局估计函数：保留原第7天结果，同时收集完整0至7天曲线。
estimate_bootstrap_outcomes <- function(data, weight) {
  # 先调用锁定函数计算第7天AKI、竞争事件和30天死亡风险。
  horizon_result <- estimate_bootstrap_outcomes_horizon_only(data, weight)
  # 把暴露安全转换为0/1整数。
  treatment <- as.integer(as.character(data$deescalation))
  # 构造锁定四状态事件因子。
  event_factor <- make_aki_event_factor_boot(
    data$aki_7d,
    data$aki_competing_state
  )
  # 当前调用计数加1；调用顺序固定为每份插补先CBPS后ATO。
  .curve_collector$call_index <- .curve_collector$call_index + 1L
  # 根据调用序号推导插补编号。
  imputation_index <- ceiling(.curve_collector$call_index / 2)
  # 奇数调用对应CBPS-ATE，偶数调用对应ATO。
  model <- if (.curve_collector$call_index %% 2L == 1L) {
    # 返回主CBPS-ATE模型标识。
    "cbps_ate"
  } else {
    # 返回ATO敏感性模型标识。
    "ato"
  }
  # 分策略计算AKI累积发生曲线。
  for (strategy in c(0L, 1L)) {
    # 选择当前策略组患者。
    index <- treatment == strategy
    # 计算当前策略的固定日网格曲线。
    curve <- weighted_aj_curve_boot(
      data$aki_followup_days[index],
      event_factor[index],
      weight[index]
    )
    # 补充插补、模型和策略标识。
    curve$imputation <- imputation_index
    # 保存当前曲线的模型标识。
    curve$model <- model
    # 保存当前曲线的策略编码。
    curve$strategy <- strategy
    # 把当前曲线追加到环境列表。
    .curve_collector$rows[[length(.curve_collector$rows) + 1L]] <- curve
  }
  # 返回锁定的第7天和死亡风险结果。
  horizon_result
}

# 扩展单次Bootstrap函数：在原结果中附加跨20份插补合并的AKI曲线。
run_single_bootstrap_repeat <- function(
  trial,
  trial_code,
  bootstrap_index,
  raw_data,
  dictionary,
  method_specification,
  predictor_specification,
  ps_formula
) {
  # 每个Bootstrap重复开始前重置曲线调用计数。
  .curve_collector$call_index <- 0L
  # 每个Bootstrap重复开始前清空旧曲线。
  .curve_collector$rows <- list()
  # 调用锁定核心完成抽样、MICE、权重、结局和第7天合并。
  result <- run_single_bootstrap_repeat_horizon_only(
    trial = trial,
    trial_code = trial_code,
    bootstrap_index = bootstrap_index,
    raw_data = raw_data,
    dictionary = dictionary,
    method_specification = method_specification,
    predictor_specification = predictor_specification,
    ps_formula = ps_formula
  )
  # 合并本次重复的全部逐插补曲线。
  curve_rows <- do.call(rbind, .curve_collector$rows)
  # 要求20份插补×2模型×2策略×8时间点完整存在。
  expected_curve_n <- bootstrap_m * 2L * 2L * 8L
  # 曲线行数或缺失不符合合同时停止。
  if (nrow(curve_rows) != expected_curve_n || anyNA(curve_rows)) {
    # 报告实际和期望曲线行数。
    stop(
      "Bootstrap曲线收集不完整；实际=",
      nrow(curve_rows),
      "，期望=",
      expected_curve_n
    )
  }
  # 按模型、策略和时间对20份插补风险取均值。
  pooled_curves <- stats::aggregate(
    aki_risk ~ model + strategy + time_days,
    data = curve_rows,
    FUN = mean
  )
  # 补充当前试验标识。
  pooled_curves$trial <- trial
  # 补充当前Bootstrap编号。
  pooled_curves$bootstrap_index <- bootstrap_index
  # 调整列顺序以便后续批量汇总。
  pooled_curves <- pooled_curves[
    ,
    c(
      "trial",
      "bootstrap_index",
      "model",
      "strategy",
      "time_days",
      "aki_risk"
    )
  ]
  # 把逐插补曲线附加到检查点，支持独立复算。
  result$imputation_curves <- curve_rows
  # 把跨20份插补合并曲线附加到检查点。
  result$pooled_curves <- pooled_curves
  # 返回包含原效应与曲线扩展的单次Bootstrap结果。
  result
}
