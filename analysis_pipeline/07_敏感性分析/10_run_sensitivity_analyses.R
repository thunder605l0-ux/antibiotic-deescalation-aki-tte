# ==============================================================================
# ST-08：预设敏感性分析
# Method: 复用每次主Bootstrap中已经完成的MICE对象，重新估计边界和权重敏感性规格
# Purpose: 避免再次插补或改变抽样种子，并为SF-09/ST-08生成点估计与95%百分位区间
# Source type: local_knowledge_base + GitHub + self_written
# Source name or URL: D:/source/obsidian；https://github.com/ngreifer/WeightIt；https://github.com/amices/mice
# Citation or commit/version: WeightIt 595d2812c2f4d56ff7c14aaf6fcde76c59aa8c8f；mice 61083e667fa41cc4bc655492591b100bf8a7e443
# Access date: 2026-07-29
# Adaptation notes: 主分析与ATO复用检查点；NEW_MS03在同一重抽样键上按锁定31列规格重新插补；日期按锁定DATETIME格式显式解析，空白值为缺失。
# Verification command: Rscript 07_敏感性分析/10_run_sensitivity_analyses.R
# ==============================================================================

# 读取统一配置。
source(
  file.path(Sys.getenv("TTE_RUNTIME_ROOT", unset = "."), "02_配置与变量字典", "analysis_config.R"),
  encoding = "UTF-8"
)
setwd(package_root)

# 加载单次Bootstrap锁定函数和br.logit适配器。
source(
  file.path(package_root, "05_Bootstrap", "04_bootstrap_single_repeat.R"),
  encoding = "UTF-8"
)

# 检查所需依赖。
require_packages(c("mice", "flextable", "officer"))

# 定义Bootstrap检查点目录。
checkpoint_dir <- file.path(output_root, "04_Bootstrap", "checkpoints")

# 定义敏感性分析输出目录。
sensitivity_dir <- file.path(output_root, "06_敏感性分析")

# 创建输出目录。
dir.create(sensitivity_dir, recursive = TRUE, showWarnings = FALSE)

# 读取锁定MICE字典与两套插补规格，供NEW_MS03和完整病例分析使用。
mice_dictionary <- utils::read.csv(
  file.path(config_dir, "mice_working_variable_dictionary.csv"),
  stringsAsFactors = FALSE,
  na.strings = ""
)
mice_method_specification <- utils::read.csv(
  file.path(config_dir, "mice_method_specification.csv"),
  stringsAsFactors = FALSE,
  na.strings = ""
)
mice_predictor_specification <- utils::read.csv(
  file.path(config_dir, "mice_predictor_matrix_directed_uo.csv"),
  stringsAsFactors = FALSE
)

# NEW_MS03以RACE_W替代SEX，其余操作变量保持不变。
ps_formula_new_ms03 <- stats::reformulate(
  c(
    setdiff(
      setdiff(all.vars(ps_formula_main), "deescalation"),
      "male"
    ),
    "race_white"
  ),
  response = "deescalation"
)
# 治疗时机代理增强分析在NEW_MS01基础上增加ABX_TIME。
ps_formula_abx_time <- stats::update(
  ps_formula_main,
  . ~ . + sepsis_to_target_start_hours
)
# 微生物代理扰动分析在NEW_MS01基础上增加T0微生物信息状态。
ps_formula_micro_state <- stats::update(
  ps_formula_main,
  . ~ . + microbiology_state_group4
)

# 为NEW_MS03构建31列工作数据，不改变锁定主MICE函数。
build_new_ms03_working_data <- function(
  bootstrap_raw,
  reference_raw,
  trial
) {
  work <- build_bootstrap_working_data(
    bootstrap_raw,
    reference_raw,
    trial,
    mice_dictionary
  )
  work$race_white <- factor(
    bootstrap_raw$race_white,
    levels = c(0, 1)
  )
  work
}

# 按指定规格重建MICE方法向量。
rebuild_sensitivity_method <- function(trial, field_order) {
  rows <- mice_method_specification[
    mice_method_specification$trial == trial &
      mice_method_specification$specification == "sensitivity_new_ms03",
    ,
    drop = FALSE
  ]
  rows <- rows[match(field_order, rows$field), , drop = FALSE]
  if (anyNA(rows$field)) {
    stop("NEW_MS03 MICE方法规格字段无法闭合。")
  }
  method <- rows$active_method
  method[is.na(method)] <- ""
  names(method) <- rows$field
  method
}

# 按指定规格重建31×31定向尿量预测矩阵。
rebuild_sensitivity_predictor_matrix <- function(trial, field_order) {
  rows <- mice_predictor_specification[
    mice_predictor_specification$trial == trial &
      mice_predictor_specification$specification ==
        "sensitivity_new_ms03",
    ,
    drop = FALSE
  ]
  rows$imputation_target <- factor(
    rows$imputation_target,
    levels = field_order
  )
  rows$predictor <- factor(rows$predictor, levels = field_order)
  matrix_object <- unclass(stats::xtabs(
    include ~ imputation_target + predictor,
    data = rows,
    drop.unused.levels = FALSE
  ))
  storage.mode(matrix_object) <- "integer"
  matrix_object[field_order, field_order, drop = FALSE]
}

# 读取患者级三路径AKI中间表，用于肌酐单路径结局敏感性分析。
aki_component_data <- utils::read.csv(
  file.path(data_dir, "aki_three_path_patient_summary.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE,
  na.strings = c("", "NA")
)

# 将肌酐单路径事件、事件时间和竞争事件随访添加到当前试验原始数据。
augment_creatinine_only_outcome <- function(raw_data, trial) {
  components <- aki_component_data[
    aki_component_data$trial == trial,
    ,
    drop = FALSE
  ]
  match_index <- match(raw_data$icu_stay_id, components$stay_id)
  if (anyNA(match_index)) {
    stop(trial, "存在无法连接到AKI三路径中间表的患者键。")
  }
  components <- components[match_index, , drop = FALSE]
  index_time <- as.POSIXct(raw_data$index_time, tz = "UTC")
  event_time <- as.POSIXct(
    components$first_creatinine_aki_time,
    tz = "UTC"
  )
  followup_end <- as.POSIXct(components$renal_followup_end, tz = "UTC")
  event_days <- as.numeric(difftime(event_time, index_time, units = "days"))
  followup_days <- as.numeric(
    difftime(followup_end, index_time, units = "days")
  )
  event <- (
    components$creatinine_aki %in% c(TRUE, "TRUE", "t", "1") &
      is.finite(event_days) &
      event_days > 0 &
      event_days <= aki_horizon_days &
      event_time <= followup_end
  )
  raw_data$aki_7d_creatinine_only <- as.integer(event)
  raw_data$aki_7d_creatinine_only_onset_days <- event_days
  raw_data$aki_followup_days_creatinine_only <- ifelse(
    event,
    event_days,
    pmin(followup_days, aki_horizon_days)
  )
  raw_data$aki_competing_state_creatinine_only <- ifelse(
    event,
    "no_pre_aki_competing_event",
    ifelse(
      components$followup_end_type == "death",
      "death_before_aki",
      ifelse(
        components$followup_end_type == "live_discharge",
        "discharge_before_aki",
        "no_pre_aki_competing_event"
      )
    )
  )
  raw_data
}

# 定义对正权重按经验分位数截尾的函数。
truncate_weights <- function(weight, lower_probability, upper_probability) {
  # 计算下限与上限。
  limits <- stats::quantile(
    weight,
    probs = c(lower_probability, upper_probability),
    names = FALSE,
    type = 7
  )
  # 将低于下限或高于上限的权重压至边界。
  pmin(pmax(weight, limits[1]), limits[2])
}

# 定义将四项核心风险派生为RD和RR的函数。
derive_effects <- function(risk_vector) {
  # 要求两种策略风险均为有限概率。
  stopifnot(all(is.finite(risk_vector)))
  # 返回固定命名的风险和效应。
  c(
    risk_aki_0 = unname(risk_vector["risk_aki_0"]),
    risk_aki_1 = unname(risk_vector["risk_aki_1"]),
    rd_aki = unname(risk_vector["risk_aki_1"] - risk_vector["risk_aki_0"]),
    rr_aki = unname(risk_vector["risk_aki_1"] / risk_vector["risk_aki_0"]),
    risk_death30_0 = unname(risk_vector["risk_death30_0"]),
    risk_death30_1 = unname(risk_vector["risk_death30_1"]),
    rd_death30 = unname(
      risk_vector["risk_death30_1"] - risk_vector["risk_death30_0"]
    ),
    rr_death30 = unname(
      risk_vector["risk_death30_1"] / risk_vector["risk_death30_0"]
    )
  )
}

# 定义在20份插补上合并指定结局和权重规则的函数。
estimate_specification <- function(
  checkpoint,
  raw_data,
  trial,
  specification
) {
  # 按检查点保存的抽样行重建Bootstrap原始样本。
  bootstrap_raw <- raw_data[checkpoint$sampled_row, , drop = FALSE]
  # 默认复用主分析检查点中的MIDS对象。
  analysis_imputation <- checkpoint$imputation
  # NEW_MS03必须使用其锁定31列插补规格重新插补race_white。
  if (specification == "NEW_MS03") {
    work_data <- build_new_ms03_working_data(
      bootstrap_raw,
      raw_data,
      trial
    )
    method <- rebuild_sensitivity_method(trial, names(work_data))
    predictor_matrix <- rebuild_sensitivity_predictor_matrix(
      trial,
      names(work_data)
    )
    urine_fields <- c(
      "urine_output_rate_6h_ml_kg_h",
      "urine_output_rate_12h_ml_kg_h",
      "urine_output_rate_24h_ml_kg_h"
    )
    active_targets <- names(method)[nzchar(method)]
    visit_sequence <- c(
      setdiff(active_targets, urine_fields),
      urine_fields
    )
    sensitivity_seed <- if (!is.null(checkpoint$mice_seed)) {
      as.integer(checkpoint$mice_seed)
    } else {
      master_seed
    }
    analysis_imputation <- mice::mice(
      work_data,
      m = mice_m,
      maxit = mice_maxit,
      method = method,
      predictorMatrix = predictor_matrix,
      visitSequence = visit_sequence,
      donors = mice_donors,
      seed = sensitivity_seed,
      printFlag = FALSE
    )
    if (!is.null(analysis_imputation$loggedEvents)) {
      stop("NEW_MS03 MICE产生logged event。")
    }
  }
  # 初始化逐插补结果列表。
  imputation_rows <- vector("list", mice_m)
  # 初始化逐插补权重诊断列表。
  diagnostic_rows <- vector("list", mice_m)
  # 逐份完成数据计算敏感性估计。
  for (imputation_index in seq_len(mice_m)) {
    # 从检查点内MIDS对象提取当前完成数据。
    completed <- mice::complete(analysis_imputation, imputation_index)
    # 把边界敏感性结局从同序原始Bootstrap样本复制进完成数据。
    completed$aki_7d_cross_t0_uo_sensitivity <- factor(
      bootstrap_raw$aki_7d_cross_t0_uo_sensitivity,
      levels = c(0, 1)
    )
    # 把死亡边界敏感性结局复制进完成数据。
    completed$death_30d_boundary_sensitivity <- factor(
      bootstrap_raw$death_30d_boundary_sensitivity,
      levels = c(0, 1)
    )
    # 添加仅用于预设代理敏感性分析的T0前结构化字段。
    completed$sepsis_to_target_start_hours <- as.numeric(
      bootstrap_raw$sepsis_to_target_start_hours
    )
    completed$microbiology_state_group4 <- factor(
      bootstrap_raw$microbiology_state_group4,
      levels = c(0, 1, 2, 3)
    )
    # 根据当前规格选择倾向评分公式。
    current_formula <- if (specification == "NEW_MS03") {
      ps_formula_new_ms03
    } else if (specification == "ABX_TIME_proxy") {
      ps_formula_abx_time
    } else if (specification == "MICRO_STATE_proxy") {
      ps_formula_micro_state
    } else {
      ps_formula_main
    }
    # 完整病例分析按原始未插补主PS字段确定分母。
    if (specification == "complete_case") {
      raw_work <- build_bootstrap_working_data(
        bootstrap_raw,
        raw_data,
        trial,
        mice_dictionary
      )
      complete_index <- stats::complete.cases(
        raw_work[, all.vars(ps_formula_main), drop = FALSE]
      )
      completed <- completed[complete_index, , drop = FALSE]
      bootstrap_raw_current <- bootstrap_raw[
        complete_index,
        ,
        drop = FALSE
      ]
      if (
        nrow(completed) == 0L ||
          !setequal(
            unique(as.integer(as.character(completed$deescalation))),
            c(0L, 1L)
          )
      ) {
        stop("完整病例Bootstrap样本缺少任一策略组。")
      }
    } else {
      bootstrap_raw_current <- bootstrap_raw
    }
    # 拟合主CBPS-ATE权重。
    weight_fit <- fit_bootstrap_weights(
      completed,
      model = "cbps_ate",
      formula = current_formula
    )
    # 默认使用未截尾主权重，并保留截尾前副本用于审计。
    untruncated_weight <- weight_fit$weight
    current_weight <- untruncated_weight
    # 1/99规格进行相应分位数截尾。
    if (specification == "weight_truncation_1_99") {
      current_weight <- truncate_weights(current_weight, 0.01, 0.99)
    }
    # 5/95规格进行相应分位数截尾。
    if (specification == "weight_truncation_5_95") {
      current_weight <- truncate_weights(current_weight, 0.05, 0.95)
    }
    # 计算当前插补实际被截尾的权重比例。
    truncated_fraction <- mean(
      abs(current_weight - untruncated_weight) >
        sqrt(.Machine$double.eps)
    )
    # 分策略计算当前敏感性权重的Kish有效样本量。
    treatment <- as.integer(as.character(completed$deescalation))
    current_ess <- vapply(
      c(0L, 1L),
      function(strategy) {
        effective_sample_size_boot(
          current_weight[treatment == strategy]
        )
      },
      numeric(1)
    )
    # 跨T0尿量边界规格替换AKI二分类结局。
    if (specification == "aki_cross_t0_uo_boundary") {
      # 解析索引时点。
      index_time <- as.POSIXct(
        bootstrap_raw$index_time,
        format = "%Y-%m-%d %H:%M:%S",
        tz = "UTC"
      )
      # 解析跨T0尿量口径下的首次AKI时点。
      cross_event_time <- as.POSIXct(
        bootstrap_raw$aki_7d_cross_t0_uo_onset_time,
        format = "%Y-%m-%d %H:%M:%S",
        tz = "UTC"
      )
      # 计算跨T0口径事件距T0的天数。
      cross_event_days <- as.numeric(
        difftime(cross_event_time, index_time, units = "days")
      )
      # 标记跨T0口径发生AKI且事件时间有效的记录。
      cross_event_index <- (
        bootstrap_raw$aki_7d_cross_t0_uo_sensitivity == 1L &
          is.finite(cross_event_days) &
          cross_event_days > 0 &
          cross_event_days <= aki_horizon_days
      )
      # 对新增或更早的AKI事件更新首次事件/随访时间。
      completed$aki_followup_days[cross_event_index] <- pmin(
        completed$aki_followup_days[cross_event_index],
        cross_event_days[cross_event_index]
      )
      # AKI先发生时竞争状态必须重置为无AKI前竞争事件。
      completed$aki_competing_state[cross_event_index] <- factor(
        "no_pre_aki_competing_event",
        levels = levels(completed$aki_competing_state)
      )
      # 最后替换跨T0口径的AKI二分类结局。
      completed$aki_7d <- completed$aki_7d_cross_t0_uo_sensitivity
    }
    # 肌酐单路径规格替换AKI事件、时间及竞争状态。
    if (specification == "creatinine_only_phenotype") {
      completed$aki_7d <- factor(
        bootstrap_raw_current$aki_7d_creatinine_only,
        levels = c(0, 1)
      )
      completed$aki_followup_days <- as.numeric(
        bootstrap_raw_current$aki_followup_days_creatinine_only
      )
      completed$aki_competing_state <- factor(
        bootstrap_raw_current$aki_competing_state_creatinine_only,
        levels = c(
          "no_pre_aki_competing_event",
          "death_before_aki",
          "discharge_before_aki"
        )
      )
    }
    # 日期级死亡边界规格替换30天死亡二分类结局。
    if (specification == "death_calendar_boundary") {
      completed$death_30d <- completed$death_30d_boundary_sensitivity
    }
    # 用锁定加权Aalen–Johansen和加权风险函数估计结局。
    estimate <- estimate_bootstrap_outcomes_horizon_only(
      completed,
      current_weight
    )
    # 保存当前插补的四项核心风险。
    imputation_rows[[imputation_index]] <- data.frame(
      risk_aki_0 = estimate$risk_aki_0,
      risk_aki_1 = estimate$risk_aki_1,
      risk_death30_0 = estimate$risk_death30_0,
      risk_death30_1 = estimate$risk_death30_1,
      stringsAsFactors = FALSE
    )
    # 保存当前插补的样本量、截尾比例和两组ESS。
    diagnostic_rows[[imputation_index]] <- data.frame(
      n_analysis = nrow(completed),
      truncated_fraction = truncated_fraction,
      ess_continuation = unname(current_ess[1]),
      ess_deescalation = unname(current_ess[2]),
      stringsAsFactors = FALSE
    )
  }
  # 合并20份插补结果。
  imputation_table <- do.call(rbind, imputation_rows)
  # 对策略风险先取跨插补均值。
  pooled_risks <- vapply(imputation_table, mean, numeric(1))
  # 合并逐插补诊断；ESS同时保留最差值以反映最弱支持。
  diagnostic_table <- do.call(rbind, diagnostic_rows)
  # 从合并风险派生RD和RR，并追加可审计诊断。
  c(
    derive_effects(pooled_risks),
    n_analysis = unique(diagnostic_table$n_analysis),
    truncated_fraction = mean(diagnostic_table$truncated_fraction),
    ess_continuation_min = min(diagnostic_table$ess_continuation),
    ess_deescalation_min = min(diagnostic_table$ess_deescalation)
  )
}

# 定义需要由主Bootstrap检查点或同一重抽样键重新计算的预设规格。
recalculated_specifications <- c(
  "aki_cross_t0_uo_boundary",
  "death_calendar_boundary",
  "weight_truncation_1_99",
  "weight_truncation_5_95",
  "creatinine_only_phenotype",
  "NEW_MS03",
  "ABX_TIME_proxy",
  "MICRO_STATE_proxy",
  "complete_case"
)

# 初始化每次重复结果列表。
replicate_rows <- list()

# 逐试验读取1000个正式检查点。
for (trial in trial_registry$trial) {
  # 读取当前试验原始分析数据。
  raw_data <- utils::read.csv(
    file.path(data_dir, paste0("analysis_dataset_", trial, ".csv")),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  raw_data <- augment_creatinine_only_outcome(raw_data, trial)
  # 构造当前试验1000个检查点路径。
  checkpoint_paths <- file.path(
    checkpoint_dir,
    sprintf("%s_b%04d_checkpoint.rds", trial, seq_len(bootstrap_B))
  )
  # 要求全部1000个检查点存在。
  if (!all(file.exists(checkpoint_paths))) {
    stop(trial, "尚未完成1000个正式Bootstrap检查点。")
  }
  # 逐个检查点进行敏感性估计。
  for (bootstrap_index in seq_len(bootstrap_B)) {
    # 读取单个成功检查点。
    checkpoint <- readRDS(checkpoint_paths[bootstrap_index])
    # 要求检查点身份和状态正确。
    stopifnot(
      identical(checkpoint$status, "success"),
      identical(checkpoint$trial, trial),
      identical(as.integer(checkpoint$bootstrap_index), bootstrap_index)
    )
    # 从检查点直接提取主CBPS与ATO结果。
    for (model_name in c("cbps_ate", "ato")) {
      # 选择当前模型的合并估计。
      pooled <- checkpoint$pooled_estimates[
        checkpoint$pooled_estimates$model == model_name,
        ,
        drop = FALSE
      ]
      # 要求当前模型只有一行。
      stopifnot(nrow(pooled) == 1L)
      current_diagnostics <- checkpoint$weight_diagnostics[
        checkpoint$weight_diagnostics$model == model_name,
        ,
        drop = FALSE
      ]
      stopifnot(nrow(current_diagnostics) == mice_m)
      # 保存可直接复用的主/ATO估计。
      replicate_rows[[length(replicate_rows) + 1L]] <- data.frame(
        trial = trial,
        bootstrap_index = bootstrap_index,
        specification = ifelse(
          model_name == "cbps_ate",
          "main_cbps_ate",
          "overlap_weight_ato_brlogit"
        ),
        pooled[, c(
          "risk_aki_0",
          "risk_aki_1",
          "rd_aki",
          "rr_aki",
          "risk_death30_0",
          "risk_death30_1",
          "rd_death30",
          "rr_death30"
        )],
        n_analysis = length(checkpoint$sampled_row),
        truncated_fraction = 0,
        ess_continuation_min = min(
          current_diagnostics$ess_continuation
        ),
        ess_deescalation_min = min(
          current_diagnostics$ess_deescalation
        ),
        stringsAsFactors = FALSE
      )
    }
    # 逐项重新计算四个敏感性规格。
    for (specification in recalculated_specifications) {
      # 运行当前规格。
      estimate <- estimate_specification(
        checkpoint,
        raw_data,
        trial,
        specification
      )
      # 保存当前重复的敏感性估计。
      replicate_rows[[length(replicate_rows) + 1L]] <- data.frame(
        trial = trial,
        bootstrap_index = bootstrap_index,
        specification = specification,
        as.list(estimate),
        stringsAsFactors = FALSE
      )
    }
  }
}

# 合并全部重复。
sensitivity_replicates <- do.call(rbind, replicate_rows)

# 保存逐重复结果。
atomic_write_csv(
  sensitivity_replicates,
  file.path(sensitivity_dir, "ST-08_sensitivity_bootstrap_replicates.csv")
)

# 定义需要报告和形成区间的效应字段。
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

# ---- 原始样本敏感性点估计 ---------------------------------------------------
# Bootstrap分布的中位数不是原始样本点估计。以下严格按照主分析合同，
# 在原始样本20份插补数据上计算每项敏感性规格的点估计。
primary_point_path <- file.path(
  output_root,
  "05_主要结局",
  "primary_point_estimates_pooled.csv"
)
if (!file.exists(primary_point_path)) {
  stop("缺少原始样本主要点估计；请先运行07_primary_outcome_tables.R。")
}
primary_points <- utils::read.csv(
  primary_point_path,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
point_rows <- list()

# 主CBPS-ATE和ATO点估计直接复用已经按原始样本20份插补合并的正式结果。
for (trial in trial_registry$trial) {
  for (model_name in c("cbps_ate", "ato")) {
    current <- primary_points[
      primary_points$trial == trial &
        primary_points$model == model_name,
      ,
      drop = FALSE
    ]
    stopifnot(nrow(current) == 1L)
    model_diagnostics <- lapply(seq_len(mice_m), function(imputation_index) {
      weight_object <- readRDS(file.path(
        output_root,
        "03_权重",
        paste0(
          "weights_",
          trial,
          "_imp",
          sprintf("%02d", imputation_index),
          ".rds"
        )
      ))
      fit_object <- if (model_name == "cbps_ate") {
        weight_object$cbps_ate
      } else {
        weight_object$ato_brlogit
      }
      data.frame(
        n_analysis = nrow(weight_object$completed_data),
        ess_continuation = fit_object$ess_continuation,
        ess_deescalation = fit_object$ess_deescalation,
        stringsAsFactors = FALSE
      )
    })
    model_diagnostics <- do.call(rbind, model_diagnostics)
    point_rows[[length(point_rows) + 1L]] <- data.frame(
      trial = trial,
      specification = ifelse(
        model_name == "cbps_ate",
        "main_cbps_ate",
        "overlap_weight_ato_brlogit"
      ),
      current[, effect_fields, drop = FALSE],
      n_analysis = unique(model_diagnostics$n_analysis),
      truncated_fraction = 0,
      ess_continuation_min = min(model_diagnostics$ess_continuation),
      ess_deescalation_min = min(model_diagnostics$ess_deescalation),
      stringsAsFactors = FALSE
    )
  }
}

# 边界和截尾规格使用原始样本MIDS对象，而不是Bootstrap分布中位数。
for (trial in trial_registry$trial) {
  raw_data <- utils::read.csv(
    file.path(data_dir, paste0("analysis_dataset_", trial, ".csv")),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  raw_data <- augment_creatinine_only_outcome(raw_data, trial)
  mice_path <- file.path(
    output_root,
    "02_MICE",
    paste0("mids_", trial, "_m20_maxit50.rds")
  )
  if (!file.exists(mice_path)) {
    stop("缺少原始样本MICE对象：", mice_path)
  }
  saved_mice <- readRDS(mice_path)
  stopifnot(
    identical(saved_mice$trial, trial),
    saved_mice$imputation$m == mice_m,
    nrow(saved_mice$raw_data) == nrow(raw_data)
  )
  original_checkpoint <- list(
    sampled_row = seq_len(nrow(raw_data)),
    imputation = saved_mice$imputation,
    mice_seed = master_seed
  )
  for (specification in recalculated_specifications) {
    estimate <- estimate_specification(
      original_checkpoint,
      raw_data,
      trial,
      specification
    )
    point_rows[[length(point_rows) + 1L]] <- data.frame(
      trial = trial,
      specification = specification,
      as.list(estimate),
      stringsAsFactors = FALSE
    )
  }
}

# 合并并保存原始样本敏感性点估计，供ST-08和SF-09共同使用。
sensitivity_points <- do.call(rbind, point_rows)
atomic_write_csv(
  sensitivity_points,
  file.path(sensitivity_dir, "ST-08_sensitivity_point_estimates.csv")
)

# 初始化区间结果。
interval_rows <- list()

# 逐试验、规格和效应字段形成type 7百分位区间。
for (trial in trial_registry$trial) {
  for (specification in unique(sensitivity_replicates$specification)) {
    current <- sensitivity_replicates[
      sensitivity_replicates$trial == trial &
        sensitivity_replicates$specification == specification,
      ,
      drop = FALSE
    ]
    stopifnot(nrow(current) == bootstrap_B)
    for (effect_field in effect_fields) {
      interval <- stats::quantile(
        current[[effect_field]],
        probs = c(0.025, 0.5, 0.975),
        names = FALSE,
        type = 7
      )
      interval_rows[[length(interval_rows) + 1L]] <- data.frame(
        trial = trial,
        specification = specification,
        estimand_field = effect_field,
        bootstrap_B = bootstrap_B,
        ci_lower = interval[1],
        bootstrap_median = interval[2],
        ci_upper = interval[3],
        stringsAsFactors = FALSE
      )
    }
  }
}

# 合并区间表。
sensitivity_intervals <- do.call(rbind, interval_rows)

# 将每项区间与同一试验、同一规格的原始样本点估计连接。
point_long_rows <- list()
for (row_index in seq_len(nrow(sensitivity_points))) {
  for (effect_field in effect_fields) {
    point_long_rows[[length(point_long_rows) + 1L]] <- data.frame(
      trial = sensitivity_points$trial[row_index],
      specification = sensitivity_points$specification[row_index],
      estimand_field = effect_field,
      point_estimate = sensitivity_points[[effect_field]][row_index],
      stringsAsFactors = FALSE
    )
  }
}
point_long <- do.call(rbind, point_long_rows)
sensitivity_intervals <- merge(
  sensitivity_intervals,
  point_long,
  by = c("trial", "specification", "estimand_field"),
  all.x = TRUE,
  sort = FALSE
)
if (anyNA(sensitivity_intervals$point_estimate)) {
  stop("敏感性区间未能全部连接到原始样本点估计。")
}

# 保存区间表。
atomic_write_csv(
  sensitivity_intervals,
  file.path(sensitivity_dir, "ST-08_sensitivity_percentile_ci.csv")
)

# 连接原始样本权重诊断，形成正式ST-08长表。
diagnostic_columns <- c(
  "trial",
  "specification",
  "n_analysis",
  "truncated_fraction",
  "ess_continuation_min",
  "ess_deescalation_min"
)
st08_table <- merge(
  sensitivity_intervals,
  sensitivity_points[, diagnostic_columns, drop = FALSE],
  by = c("trial", "specification"),
  all.x = TRUE,
  sort = FALSE
)
st08_table$risk_set_rebuilt <- ifelse(
  st08_table$specification == "complete_case",
  "Yes: complete-case denominator",
  "No"
)
st08_table$estimate_95ci <- ifelse(
  grepl("^rr_", st08_table$estimand_field),
  sprintf(
    "%.3f (%.3f to %.3f)",
    st08_table$point_estimate,
    st08_table$ci_lower,
    st08_table$ci_upper
  ),
  sprintf(
    "%.3f (%.3f to %.3f)",
    st08_table$point_estimate,
    st08_table$ci_lower,
    st08_table$ci_upper
  )
)
atomic_write_csv(
  st08_table,
  file.path(sensitivity_dir, "ST-08_prespecified_sensitivity_results.csv")
)

# 生成可直接用于补充材料的单一三线表。
st08_display <- st08_table[, c(
  "trial",
  "specification",
  "estimand_field",
  "n_analysis",
  "risk_set_rebuilt",
  "truncated_fraction",
  "ess_continuation_min",
  "ess_deescalation_min",
  "estimate_95ci"
)]
names(st08_display) <- c(
  "Trial",
  "Specification",
  "Estimand",
  "N",
  "Risk set rebuilt",
  "Truncated fraction",
  "Minimum ESS: continuation",
  "Minimum ESS: de-escalation",
  "Estimate (95% CI)"
)
st08_flextable <- flextable::flextable(st08_display)
st08_flextable <- flextable::theme_booktabs(st08_flextable)
st08_flextable <- flextable::autofit(st08_flextable)
st08_document <- officer::read_docx()
st08_document <- officer::body_add_par(
  st08_document,
  "Supplementary Table 8. Prespecified sensitivity analyses",
  style = "heading 1"
)
st08_document <- flextable::body_add_flextable(
  st08_document,
  st08_flextable
)
print(
  st08_document,
  target = file.path(
    sensitivity_dir,
    "ST-08_prespecified_sensitivity_results.docx"
  )
)

# 登记无法仅凭当前640键数据重建的计划规格，禁止伪造结果。
unavailable_specifications <- data.frame(
  specification = "active_order_eligibility",
  status = "requires_separately_locked_input_or_model",
  reason = paste(
    "640键主数据不含活动医嘱资格敏感性扩展人群；",
    "该分析必须重新构建风险集、结局、协变量和权重，",
    "不得在主队列中伪造。"
  ),
  stringsAsFactors = FALSE
)

# 保存未执行规格登记，避免把缺失分析误写为阴性结果。
atomic_write_csv(
  unavailable_specifications,
  file.path(sensitivity_dir, "ST-08_unavailable_prespecified_specifications.csv")
)

# 输出完成信息。
message("ST-08可由当前运行包执行的6项规格已完成；另6项已登记所需额外输入。")
