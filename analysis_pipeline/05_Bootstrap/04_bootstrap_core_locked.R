# ==============================================================================
# Analysis: Phase 5.1B 单次Boot-MI重复执行引擎
# Date: 2026-07-29
# Random seed: 由trial_code与bootstrap重复号按锁定公式确定
# R: 4.4.2
#
# Method:
# 在一个平行试验内，从原始不完整患者—试验数据有放回抽样；随后重新执行
# m=20、maxit=50的定向尿量MICE，并在每份完成数据中重新拟合线性CBPS-ATE
# 和logistic ATO权重，最后估计7天AKI加权Aalen–Johansen风险及30天死亡风险。
# Purpose:
# 为正式B=1000患者级非参数bootstrap提供可复用、可审计、失败即停止的单次引擎。
# Source type:
# project_locked_protocol / project_locked_code / package_documentation
# Source name or URL:
# phase47K_outcome_mi_bootstrap_execution_protocol_v1.0_confirmed.md
# 94_phase44_execute_mice_and_diagnostics.R
# 107_phase47F_full_blinded_cbps_weights.R
# 112_phase47I_fit_ato_overlap_weights.R
# 122_phase51_primary_outcome_point_estimates.R
# https://amices.org/mice/reference/mice.html
# https://ngreifer.github.io/WeightIt/reference/method_cbps.html
# https://stat.ethz.ch/R-manual/R-devel/library/survival/html/survfit.formula.html
# Citation or commit/version:
# mice 3.19.10, commit 61083e667fa41cc4bc655492591b100bf8a7e443;
# WeightIt 1.7.0.9004, commit 595d2812c2f4d56ff7c14aaf6fcde76c59aa8c8f;
# survival 3.7-0; D-4.7K-OUTCOME-MI-BOOTSTRAP-EXECUTION-CONTRACT-V1.0
# Access date: 2026-07-29
# Adaptation notes:
# 仅执行NEW_MS01；保留原icu_stay_id追溯，同时用bootstrap_copy_id作为重复内唯一键。
# SMD、方差比和ESS只保存为诊断，不新增为单次bootstrap失败门槛。
# 对MICE内部remove.lindep()记录的两类已验证事件进行精确白名单审计：
# 24小时尿量PMM中的有序抗菌药.L/.Q对比列，以及目标观测子集中计数确为0的
# AKI前死亡竞争事件虚拟列；其他loggedEvents仍硬停止。
# anti-PSA的source_control_status仍使用锁定logreg统计步骤；仅在精确
# glm.fit不收敛警告下，以相同增广数据和起始规则提高内部迭代上限重试一次。
# Verification command:
# D:/software/R4.4.2/R-4.4.2/bin/x64/Rscript.exe --vanilla
#   02_code/125_phase51B_bootstrap_preflight.R
# ==============================================================================

# 逐行说明：将本行计算结果保存到对象 bootstrap_master_seed，供后续校验或估计使用。
bootstrap_master_seed <- 20260727L
# 逐行说明：将本行计算结果保存到对象 bootstrap_m，供后续校验或估计使用。
bootstrap_m <- 20L
# 逐行说明：将本行计算结果保存到对象 bootstrap_maxit，供后续校验或估计使用。
bootstrap_maxit <- 50L
# 逐行说明：将本行计算结果保存到对象 bootstrap_donors，供后续校验或估计使用。
bootstrap_donors <- 5L
# 逐行说明：将本行计算结果保存到对象 bootstrap_horizon_aki，供后续校验或估计使用。
bootstrap_horizon_aki <- 7
# 逐行说明：将本行计算结果保存到对象 bootstrap_tolerance，供后续校验或估计使用。
bootstrap_tolerance <- 1e-10

# 逐行说明：定义函数 derive_target_class_boot，封装本步骤的单一可复用逻辑。
derive_target_class_boot <- function(code, trial) {
  # 逐行说明：检查本行条件；条件成立时执行后续保护或分支逻辑。
  if (trial == "anti_mrsa") {
    # 逐行说明：立即返回括号内结果，结束当前函数。
    return(factor(
      # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
      ifelse(
        # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
        code %in% c(5, 8, 15),
        # 逐行说明：声明当前字符值、字段名或因子水平。
        "vancomycin_based",
        # 逐行说明：声明当前字符值、字段名或因子水平。
        "non_vancomycin_anti_mrsa"
      # 逐行说明：结束当前多行函数调用或表达式。
      ),
      # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
      levels = c("vancomycin_based", "non_vancomycin_anti_mrsa")
    # 逐行说明：结束当前多行函数调用或表达式。
    ))
  # 逐行说明：结束当前函数、循环或条件分支。
  }
  # 逐行说明：检查本行条件；条件成立时执行后续保护或分支逻辑。
  if (trial == "anti_psa") {
    # 逐行说明：立即返回括号内结果，结束当前函数。
    return(factor(
      # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
      ifelse(
        # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
        code %in% c(1, 2),
        # 逐行说明：声明当前字符值、字段名或因子水平。
        "antipseudomonal_cephalosporin",
        # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
        ifelse(
          # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
          code == 13,
          # 逐行说明：声明当前字符值、字段名或因子水平。
          "piperacillin_tazobactam",
          # 逐行说明：声明当前字符值、字段名或因子水平。
          "other_antipseudomonal"
        # 逐行说明：结束当前多行函数调用或表达式。
        )
      # 逐行说明：结束当前多行函数调用或表达式。
      ),
      # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
      levels = c(
        # 逐行说明：声明当前字符值、字段名或因子水平。
        "antipseudomonal_cephalosporin",
        # 逐行说明：声明当前字符值、字段名或因子水平。
        "piperacillin_tazobactam",
        # 逐行说明：声明当前字符值、字段名或因子水平。
        "other_antipseudomonal"
      # 逐行说明：结束当前多行函数调用或表达式。
      )
    # 逐行说明：结束当前多行函数调用或表达式。
    ))
  # 逐行说明：结束当前函数、循环或条件分支。
  }
  # 逐行说明：触发硬停止并报告原因，防止失败对象进入后续分析。
  stop("未知试验：", trial)
# 逐行说明：结束当前函数、循环或条件分支。
}

# 逐行说明：定义函数 build_bootstrap_working_data，封装本步骤的单一可复用逻辑。
build_bootstrap_working_data <- function(
  # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
  bootstrap_raw,
  # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
  reference_raw,
  # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
  trial,
  # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
  dictionary
# 逐行说明：结束当前多行函数调用或表达式。
) {
  # 逐行说明：将本行计算结果保存到对象 selected_dictionary，供后续校验或估计使用。
  selected_dictionary <- dictionary[dictionary$role != "new_ms03_only", ]
  # 逐行说明：将本行计算结果保存到对象 work，供后续校验或估计使用。
  work <- data.frame(matrix(
    # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
    nrow = nrow(bootstrap_raw),
    # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
    ncol = 0L,
    # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
    dimnames = list(NULL, character(0))
  # 逐行说明：结束当前多行函数调用或表达式。
  ))

  # 逐行说明：按本行定义的索引逐项循环，确保每个试验、插补或策略均被处理。
  for (field_index in seq_len(nrow(selected_dictionary))) {
    # 逐行说明：将本行计算结果保存到对象 field，供后续校验或估计使用。
    field <- selected_dictionary$field[field_index]
    # 逐行说明：将本行计算结果保存到对象 source_field，供后续校验或估计使用。
    source_field <- selected_dictionary$source_field[field_index]
    # 逐行说明：将本行计算结果保存到对象 storage_type，供后续校验或估计使用。
    storage_type <- selected_dictionary$storage_type[field_index]

    # 逐行说明：检查本行条件；条件成立时执行后续保护或分支逻辑。
    if (field == "target_class_grouped") {
      # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
      work[[field]] <- derive_target_class_boot(
        # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
        bootstrap_raw[[source_field]],
        # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
        trial
      # 逐行说明：结束当前多行函数调用或表达式。
      )
    # 逐行说明：结束上一分支并进入互斥的替代分支。
    } else if (field == "aki_competing_state") {
      # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
      work[[field]] <- factor(
        # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
        ifelse(
          # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
          bootstrap_raw[[source_field]] == 2,
          # 逐行说明：声明当前字符值、字段名或因子水平。
          "death_before_aki",
          # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
          ifelse(
            # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
            bootstrap_raw[[source_field]] == 3,
            # 逐行说明：声明当前字符值、字段名或因子水平。
            "discharge_before_aki",
            # 逐行说明：声明当前字符值、字段名或因子水平。
            "no_pre_aki_competing_event"
          # 逐行说明：结束当前多行函数调用或表达式。
          )
        # 逐行说明：结束当前多行函数调用或表达式。
        ),
        # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
        levels = c(
          # 逐行说明：声明当前字符值、字段名或因子水平。
          "no_pre_aki_competing_event",
          # 逐行说明：声明当前字符值、字段名或因子水平。
          "death_before_aki",
          # 逐行说明：声明当前字符值、字段名或因子水平。
          "discharge_before_aki"
        # 逐行说明：结束当前多行函数调用或表达式。
        )
      # 逐行说明：结束当前多行函数调用或表达式。
      )
    # 逐行说明：结束上一分支并进入互斥的替代分支。
    } else {
      # 逐行说明：将本行计算结果保存到对象 value，供后续校验或估计使用。
      value <- bootstrap_raw[[source_field]]
      # 逐行说明：检查本行条件；条件成立时执行后续保护或分支逻辑。
      if (field %in% c(
        # 逐行说明：声明当前字符值、字段名或因子水平。
        "other_broad_spectrum_antibiotics",
        # 逐行说明：声明当前字符值、字段名或因子水平。
        "source_control_status",
        # 逐行说明：声明当前字符值、字段名或因子水平。
        "nephrotoxic_drug_exposure"
      # 逐行说明：结束当前多行函数调用或表达式。
      )) {
        # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
        value[value == 9] <- NA
      # 逐行说明：结束当前函数、循环或条件分支。
      }

      # 逐行说明：检查本行条件；条件成立时执行后续保护或分支逻辑。
      if (storage_type == "factor") {
        # 逐行说明：将本行计算结果保存到对象 reference_value，供后续校验或估计使用。
        reference_value <- reference_raw[[source_field]]
        # 逐行说明：检查本行条件；条件成立时执行后续保护或分支逻辑。
        if (field %in% c(
          # 逐行说明：声明当前字符值、字段名或因子水平。
          "other_broad_spectrum_antibiotics",
          # 逐行说明：声明当前字符值、字段名或因子水平。
          "source_control_status",
          # 逐行说明：声明当前字符值、字段名或因子水平。
          "nephrotoxic_drug_exposure"
        # 逐行说明：结束当前多行函数调用或表达式。
        )) {
          # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
          reference_value[reference_value == 9] <- NA
        # 逐行说明：结束当前函数、循环或条件分支。
        }
        # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
        work[[field]] <- factor(
          # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
          value,
          # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
          levels = sort(unique(reference_value[!is.na(reference_value)]))
        # 逐行说明：结束当前多行函数调用或表达式。
        )
      # 逐行说明：结束上一分支并进入互斥的替代分支。
      } else if (storage_type == "ordered") {
        # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
        work[[field]] <- ordered(value, levels = c(0, 1, 2))
      # 逐行说明：结束上一分支并进入互斥的替代分支。
      } else if (storage_type == "numeric") {
        # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
        work[[field]] <- as.numeric(value)
      # 逐行说明：结束上一分支并进入互斥的替代分支。
      } else {
        # 逐行说明：触发硬停止并报告原因，防止失败对象进入后续分析。
        stop("未知存储类型：", storage_type)
      # 逐行说明：结束当前函数、循环或条件分支。
      }
    # 逐行说明：结束当前函数、循环或条件分支。
    }
  # 逐行说明：结束当前函数、循环或条件分支。
  }
  # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
  work
# 逐行说明：结束当前函数、循环或条件分支。
}

# 逐行说明：定义函数 rebuild_bootstrap_method，封装本步骤的单一可复用逻辑。
rebuild_bootstrap_method <- function(
  # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
  method_specification,
  # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
  trial,
  # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
  field_order
# 逐行说明：结束当前多行函数调用或表达式。
) {
  # 逐行说明：将本行计算结果保存到对象 rows，供后续校验或估计使用。
  rows <- method_specification[
    # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
    method_specification$trial == trial &
      # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
      method_specification$specification == "main_new_ms01",
  # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
  ]
  # 逐行说明：将本行计算结果保存到对象 rows，供后续校验或估计使用。
  rows <- rows[match(field_order, rows$field), ]
  # 逐行说明：检查本行条件；条件成立时执行后续保护或分支逻辑。
  if (anyNA(rows$field)) stop("MICE方法规格字段无法闭合。")
  # 逐行说明：将本行计算结果保存到对象 method，供后续校验或估计使用。
  method <- rows$active_method
  # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
  method[is.na(method)] <- ""
  # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
  names(method) <- rows$field
  # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
  method
# 逐行说明：结束当前函数、循环或条件分支。
}

# 逐行说明：定义函数 rebuild_bootstrap_predictor_matrix，封装本步骤的单一可复用逻辑。
rebuild_bootstrap_predictor_matrix <- function(
  # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
  predictor_specification,
  # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
  trial,
  # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
  field_order
# 逐行说明：结束当前多行函数调用或表达式。
) {
  # 逐行说明：将本行计算结果保存到对象 rows，供后续校验或估计使用。
  rows <- predictor_specification[
    # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
    predictor_specification$trial == trial &
      # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
      predictor_specification$specification == "main_new_ms01",
  # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
  ]
  # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
  rows$imputation_target <- factor(
    # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
    rows$imputation_target,
    # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
    levels = field_order
  # 逐行说明：结束当前多行函数调用或表达式。
  )
  # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
  rows$predictor <- factor(rows$predictor, levels = field_order)
  # 逐行说明：将本行计算结果保存到对象 matrix_object，供后续校验或估计使用。
  matrix_object <- stats::xtabs(
    # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
    include ~ imputation_target + predictor,
    # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
    data = rows,
    # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
    drop.unused.levels = FALSE
  # 逐行说明：结束当前多行函数调用或表达式。
  )
  # 逐行说明：将本行计算结果保存到对象 matrix_object，供后续校验或估计使用。
  matrix_object <- unclass(matrix_object)
  # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
  storage.mode(matrix_object) <- "integer"
  # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
  matrix_object[field_order, field_order, drop = FALSE]
# 逐行说明：结束当前函数、循环或条件分支。
}

# 逐行说明：定义函数 make_aki_event_factor_boot，封装本步骤的单一可复用逻辑。
make_aki_event_factor_boot <- function(aki_7d, competing_state) {
  # 逐行说明：将本行计算结果保存到对象 aki_7d，供后续校验或估计使用。
  aki_7d <- as.integer(as.character(aki_7d))
  # 逐行说明：将本行计算结果保存到对象 competing_state，供后续校验或估计使用。
  competing_state <- as.character(competing_state)
  # 逐行说明：将本行计算结果保存到对象 allowed_competing，供后续校验或估计使用。
  allowed_competing <- c(
    # 逐行说明：声明当前字符值、字段名或因子水平。
    "no_pre_aki_competing_event",
    # 逐行说明：声明当前字符值、字段名或因子水平。
    "death_before_aki",
    # 逐行说明：声明当前字符值、字段名或因子水平。
    "discharge_before_aki"
  # 逐行说明：结束当前多行函数调用或表达式。
  )
  # 逐行说明：检查本行条件；条件成立时执行后续保护或分支逻辑。
  if (
    # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
    length(aki_7d) != length(competing_state) ||
      # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
      anyNA(aki_7d) ||
      # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
      anyNA(competing_state) ||
      # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
      any(!aki_7d %in% c(0L, 1L)) ||
      # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
      any(!competing_state %in% allowed_competing) ||
      # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
      any(
        # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
        aki_7d == 1L &
          # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
          competing_state != "no_pre_aki_competing_event"
      # 逐行说明：结束当前多行函数调用或表达式。
      )
  # 逐行说明：结束当前多行函数调用或表达式。
  ) {
    # 逐行说明：触发硬停止并报告原因，防止失败对象进入后续分析。
    stop("AKI或竞争事件编码非法。")
  # 逐行说明：结束当前函数、循环或条件分支。
  }
  # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
  factor(
    # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
    ifelse(
      # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
      aki_7d == 1L,
      # 逐行说明：声明当前字符值、字段名或因子水平。
      "aki",
      # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
      ifelse(
        # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
        competing_state == "death_before_aki",
        # 逐行说明：声明当前字符值、字段名或因子水平。
        "death_before_aki",
        # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
        ifelse(
          # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
          competing_state == "discharge_before_aki",
          # 逐行说明：声明当前字符值、字段名或因子水平。
          "discharge_before_aki",
          # 逐行说明：声明当前字符值、字段名或因子水平。
          "censor"
        # 逐行说明：结束当前多行函数调用或表达式。
        )
      # 逐行说明：结束当前多行函数调用或表达式。
      )
    # 逐行说明：结束当前多行函数调用或表达式。
    ),
    # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
    levels = c(
      # 逐行说明：声明当前字符值、字段名或因子水平。
      "censor",
      # 逐行说明：声明当前字符值、字段名或因子水平。
      "aki",
      # 逐行说明：声明当前字符值、字段名或因子水平。
      "death_before_aki",
      # 逐行说明：声明当前字符值、字段名或因子水平。
      "discharge_before_aki"
    # 逐行说明：结束当前多行函数调用或表达式。
    )
  # 逐行说明：结束当前多行函数调用或表达式。
  )
# 逐行说明：结束当前函数、循环或条件分支。
}

# 逐行说明：定义函数 validate_aki_boot，封装本步骤的单一可复用逻辑。
validate_aki_boot <- function(time, event_factor, weight) {
  # 逐行说明：检查本行条件；条件成立时执行后续保护或分支逻辑。
  if (
    # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
    length(time) != length(event_factor) ||
      # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
      length(time) != length(weight) ||
      # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
      anyNA(time) ||
      # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
      anyNA(event_factor) ||
      # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
      anyNA(weight) ||
      # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
      any(!is.finite(time)) ||
      # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
      any(time < 0 | time > bootstrap_horizon_aki) ||
      # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
      any(!is.finite(weight)) ||
      # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
      any(weight <= 0) ||
      # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
      any(time == 0 & as.character(event_factor) != "censor")
  # 逐行说明：结束当前多行函数调用或表达式。
  ) {
    # 逐行说明：触发硬停止并报告原因，防止失败对象进入后续分析。
    stop("AKI时间、状态或权重QA失败。")
  # 逐行说明：结束当前函数、循环或条件分支。
  }
  # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
  invisible(TRUE)
# 逐行说明：结束当前函数、循环或条件分支。
}

# 逐行说明：定义函数 weighted_aj_boot，封装本步骤的单一可复用逻辑。
weighted_aj_boot <- function(time, event_factor, weight) {
  # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
  validate_aki_boot(time, event_factor, weight)
  # 逐行说明：将本行计算结果保存到对象 fit，供后续校验或估计使用。
  fit <- survival::survfit(
    # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
    survival::Surv(time, event_factor) ~ 1,
    # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
    data = data.frame(time, event_factor, weight),
    # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
    weights = weight,
    # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
    stype = 1,
    # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
    ctype = 1,
    # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
    se.fit = FALSE,
    # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
    time0 = TRUE,
    # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
    timefix = FALSE
  # 逐行说明：结束当前多行函数调用或表达式。
  )
  # 逐行说明：将本行计算结果保存到对象 index_at_horizon，供后续校验或估计使用。
  index_at_horizon <- max(which(fit$time <= bootstrap_horizon_aki))
  # 逐行说明：将本行计算结果保存到对象 probability，供后续校验或估计使用。
  probability <- fit$pstate[index_at_horizon, , drop = TRUE]
  # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
  names(probability) <- colnames(fit$pstate)
  # 逐行说明：检查本行条件；条件成立时执行后续保护或分支逻辑。
  if (
    # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
    any(!is.finite(probability)) ||
      # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
      any(probability < -bootstrap_tolerance) ||
      # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
      any(probability > 1 + bootstrap_tolerance) ||
      # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
      abs(sum(probability) - 1) > bootstrap_tolerance
  # 逐行说明：结束当前多行函数调用或表达式。
  ) {
    # 逐行说明：触发硬停止并报告原因，防止失败对象进入后续分析。
    stop("Aalen–Johansen状态概率QA失败。")
  # 逐行说明：结束当前函数、循环或条件分支。
  }
  # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
  probability
# 逐行说明：结束当前函数、循环或条件分支。
}

# 逐行说明：定义函数 weighted_aj_manual_boot，封装本步骤的单一可复用逻辑。
weighted_aj_manual_boot <- function(time, event_factor, weight) {
  # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
  validate_aki_boot(time, event_factor, weight)
  # 逐行说明：将本行计算结果保存到对象 event_character，供后续校验或估计使用。
  event_character <- as.character(event_factor)
  # 逐行说明：将本行计算结果保存到对象 terminal_states，供后续校验或估计使用。
  terminal_states <- c(
    # 逐行说明：声明当前字符值、字段名或因子水平。
    "aki",
    # 逐行说明：声明当前字符值、字段名或因子水平。
    "death_before_aki",
    # 逐行说明：声明当前字符值、字段名或因子水平。
    "discharge_before_aki"
  # 逐行说明：结束当前多行函数调用或表达式。
  )
  # 逐行说明：将本行计算结果保存到对象 event_times，供后续校验或估计使用。
  event_times <- sort(unique(time[
    # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
    event_character != "censor" &
      # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
      time <= bootstrap_horizon_aki
  # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
  ]))
  # 逐行说明：将本行计算结果保存到对象 probability，供后续校验或估计使用。
  probability <- c(
    # 逐行说明：声明当前字符值、字段名或因子水平。
    "(s0)" = 1,
    # 逐行说明：声明当前字符值、字段名或因子水平。
    "aki" = 0,
    # 逐行说明：声明当前字符值、字段名或因子水平。
    "death_before_aki" = 0,
    # 逐行说明：声明当前字符值、字段名或因子水平。
    "discharge_before_aki" = 0
  # 逐行说明：结束当前多行函数调用或表达式。
  )
  # 逐行说明：按本行定义的索引逐项循环，确保每个试验、插补或策略均被处理。
  for (current_time in event_times) {
    # 逐行说明：将本行计算结果保存到对象 risk_weight，供后续校验或估计使用。
    risk_weight <- sum(weight[time >= current_time])
    # 逐行说明：检查本行条件；条件成立时执行后续保护或分支逻辑。
    if (!is.finite(risk_weight) || risk_weight <= 0) {
      # 逐行说明：触发硬停止并报告原因，防止失败对象进入后续分析。
      stop("手工Aalen–Johansen风险集权重非法。")
    # 逐行说明：结束当前函数、循环或条件分支。
    }
    # 逐行说明：将本行计算结果保存到对象 event_weight，供后续校验或估计使用。
    event_weight <- vapply(
      # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
      terminal_states,
      # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
      function(state) {
        # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
        sum(weight[
          # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
          time == current_time &
            # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
            event_character == state
        # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
        ])
      # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
      },
      # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
      numeric(1)
    # 逐行说明：结束当前多行函数调用或表达式。
    )
    # 逐行说明：将本行计算结果保存到对象 increment，供后续校验或估计使用。
    increment <- event_weight / risk_weight
    # 逐行说明：将本行计算结果保存到对象 initial_before，供后续校验或估计使用。
    initial_before <- unname(probability["(s0)"])
    # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
    probability[terminal_states] <-
      # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
      probability[terminal_states] + initial_before * increment
    # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
    probability["(s0)"] <-
      # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
      initial_before * (1 - sum(increment))
  # 逐行说明：结束当前函数、循环或条件分支。
  }
  # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
  probability
# 逐行说明：结束当前函数、循环或条件分支。
}

# 逐行说明：定义函数 weighted_binary_risk_boot，封装本步骤的单一可复用逻辑。
weighted_binary_risk_boot <- function(outcome, weight) {
  # 逐行说明：将本行计算结果保存到对象 outcome，供后续校验或估计使用。
  outcome <- as.integer(as.character(outcome))
  # 逐行说明：检查本行条件；条件成立时执行后续保护或分支逻辑。
  if (
    # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
    length(outcome) != length(weight) ||
      # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
      anyNA(outcome) ||
      # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
      anyNA(weight) ||
      # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
      any(!outcome %in% c(0L, 1L)) ||
      # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
      any(!is.finite(weight)) ||
      # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
      any(weight <= 0)
  # 逐行说明：结束当前多行函数调用或表达式。
  ) {
    # 逐行说明：触发硬停止并报告原因，防止失败对象进入后续分析。
    stop("30天死亡或权重输入非法。")
  # 逐行说明：结束当前函数、循环或条件分支。
  }
  # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
  sum(weight * outcome) / sum(weight)
# 逐行说明：结束当前函数、循环或条件分支。
}

# 逐行说明：定义函数 extract_cbps_convergence_boot，封装本步骤的单一可复用逻辑。
extract_cbps_convergence_boot <- function(weightit_object) {
  # 逐行说明：将本行计算结果保存到对象 candidates，供后续校验或估计使用。
  candidates <- c(
    # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
    weightit_object$obj$convergence,
    # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
    weightit_object$obj$converge
  # 逐行说明：结束当前多行函数调用或表达式。
  )
  # 逐行说明：将本行计算结果保存到对象 candidates，供后续校验或估计使用。
  candidates <- candidates[!vapply(candidates, is.null, logical(1))]
  # 逐行说明：检查本行条件；条件成立时执行后续保护或分支逻辑。
  if (length(candidates) == 0L) return(NA_integer_)
  # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
  as.integer(candidates[[1]])
# 逐行说明：结束当前函数、循环或条件分支。
}

# 逐行说明：定义函数 effective_sample_size_boot，封装本步骤的单一可复用逻辑。
effective_sample_size_boot <- function(weight) {
  # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
  sum(weight)^2 / sum(weight^2)
# 逐行说明：结束当前函数、循环或条件分支。
}

# 逐行说明：定义函数 fit_bootstrap_weights，封装本步骤的单一可复用逻辑。
fit_bootstrap_weights <- function(data, formula, model) {
  # 逐行说明：将本行计算结果保存到对象 treatment，供后续校验或估计使用。
  treatment <- as.integer(as.character(data$deescalation))
  # 逐行说明：检查本行条件；条件成立时执行后续保护或分支逻辑。
  if (
    # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
    anyNA(treatment) ||
      # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
      !setequal(unique(treatment), c(0L, 1L))
  # 逐行说明：结束当前多行函数调用或表达式。
  ) {
    # 逐行说明：触发硬停止并报告原因，防止失败对象进入后续分析。
    stop("bootstrap样本或完成数据失去一个策略组。")
  # 逐行说明：结束当前函数、循环或条件分支。
  }

  # 逐行说明：将本行计算结果保存到对象 captured_warnings，供后续校验或估计使用。
  captured_warnings <- character(0)
  # 逐行说明：检查本行条件；条件成立时执行后续保护或分支逻辑。
  if (model == "cbps_ate") {
    # 逐行说明：将本行计算结果保存到对象 fit，供后续校验或估计使用。
    fit <- withCallingHandlers(
      # 逐行说明：调用WeightIt拟合锁定的倾向评分权重模型。
      WeightIt::weightit(
        # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
        formula = formula,
        # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
        data = data,
        # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
        method = "cbps",
        # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
        estimand = "ATE",
        # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
        stabilize = FALSE,
        # 固定使用恰好识别CBPS；over不代表重叠或阳性性诊断开关。
        over = FALSE,
        # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
        include.obj = TRUE
      # 逐行说明：结束当前多行函数调用或表达式。
      ),
      # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
      warning = function(condition) {
        # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
        captured_warnings <<- c(
          # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
          captured_warnings,
          # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
          conditionMessage(condition)
        # 逐行说明：结束当前多行函数调用或表达式。
        )
        # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
        invokeRestart("muffleWarning")
      # 逐行说明：结束当前函数、循环或条件分支。
      }
    # 逐行说明：结束当前多行函数调用或表达式。
    )
    # 逐行说明：将本行计算结果保存到对象 propensity_score，供后续校验或估计使用。
    propensity_score <- as.numeric(fit$ps)
    # 逐行说明：将本行计算结果保存到对象 raw_weight，供后续校验或估计使用。
    raw_weight <- as.numeric(fit$weights)
    # 逐行说明：将本行计算结果保存到对象 marginal_treated，供后续校验或估计使用。
    marginal_treated <- mean(treatment == 1L)
    # 逐行说明：将本行计算结果保存到对象 weight，供后续校验或估计使用。
    weight <- ifelse(
      # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
      treatment == 1L,
      # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
      raw_weight * marginal_treated,
      # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
      raw_weight * (1 - marginal_treated)
    # 逐行说明：结束当前多行函数调用或表达式。
    )
    # 逐行说明：将本行计算结果保存到对象 converged，供后续校验或估计使用。
    converged <- identical(extract_cbps_convergence_boot(fit), 0L)
  # 逐行说明：结束上一分支并进入互斥的替代分支。
  } else if (model == "ato") {
    # 逐行说明：将本行计算结果保存到对象 fit，供后续校验或估计使用。
    fit <- withCallingHandlers(
      # 逐行说明：调用WeightIt拟合锁定的倾向评分权重模型。
      WeightIt::weightit(
        # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
        formula = formula,
        # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
        data = data,
        # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
        method = "glm",
        # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
        estimand = "ATO",
        # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
        link = "logit",
        # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
        stabilize = FALSE,
        # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
        include.obj = TRUE
      # 逐行说明：结束当前多行函数调用或表达式。
      ),
      # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
      warning = function(condition) {
        # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
        captured_warnings <<- c(
          # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
          captured_warnings,
          # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
          conditionMessage(condition)
        # 逐行说明：结束当前多行函数调用或表达式。
        )
        # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
        invokeRestart("muffleWarning")
      # 逐行说明：结束当前函数、循环或条件分支。
      }
    # 逐行说明：结束当前多行函数调用或表达式。
    )
    # 逐行说明：将本行计算结果保存到对象 propensity_score，供后续校验或估计使用。
    propensity_score <- as.numeric(fit$ps)
    # 逐行说明：将本行计算结果保存到对象 weight，供后续校验或估计使用。
    weight <- as.numeric(fit$weights)
    # 逐行说明：将本行计算结果保存到对象 converged，供后续校验或估计使用。
    converged <- isTRUE(fit$obj$converged) &&
      # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
      all(is.finite(stats::coef(fit$obj)))
  # 逐行说明：结束上一分支并进入互斥的替代分支。
  } else {
    # 逐行说明：触发硬停止并报告原因，防止失败对象进入后续分析。
    stop("未知权重模型：", model)
  # 逐行说明：结束当前函数、循环或条件分支。
  }

  # 逐行说明：检查本行条件；条件成立时执行后续保护或分支逻辑。
  if (
    # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
    length(captured_warnings) > 0L ||
      # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
      !converged ||
      # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
      length(propensity_score) != nrow(data) ||
      # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
      any(!is.finite(propensity_score)) ||
      # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
      any(propensity_score <= 0 | propensity_score >= 1) ||
      # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
      length(weight) != nrow(data) ||
      # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
      any(!is.finite(weight)) ||
      # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
      any(weight <= 0)
  # 逐行说明：结束当前多行函数调用或表达式。
  ) {
    # 逐行说明：触发硬停止并报告原因，防止失败对象进入后续分析。
    stop(
      # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
      model,
      # 逐行说明：声明当前字符值、字段名或因子水平。
      "拟合或权重QA失败；warnings=",
      # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
      paste(captured_warnings, collapse = " | ")
    # 逐行说明：结束当前多行函数调用或表达式。
    )
  # 逐行说明：结束当前函数、循环或条件分支。
  }

  # 逐行说明：将本行计算结果保存到对象 balance，供后续校验或估计使用。
  balance <- cobalt::bal.tab(
    # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
    formula,
    # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
    data = data,
    # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
    weights = weight,
    # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
    estimand = if (model == "cbps_ate") "ATE" else "ATO",
    # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
    s.d.denom = "pooled",
    # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
    binary = "std",
    # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
    continuous = "std",
    # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
    un = TRUE
  # 逐行说明：结束当前多行函数调用或表达式。
  )$Balance
  # 逐行说明：将本行计算结果保存到对象 balance，供后续校验或估计使用。
  balance <- balance[balance$Type != "Distance", , drop = FALSE]
  # 逐行说明：将本行计算结果保存到对象 absolute_smd，供后续校验或估计使用。
  absolute_smd <- abs(balance$Diff.Adj)
  # 逐行说明：将本行计算结果保存到对象 maximum_smd_index，供后续校验或估计使用。
  maximum_smd_index <- which.max(absolute_smd)
  # 逐行说明：将本行计算结果保存到对象 ess，供后续校验或估计使用。
  ess <- vapply(
    # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
    c(0L, 1L),
    # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
    function(strategy) {
      # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
      effective_sample_size_boot(weight[treatment == strategy])
    # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
    },
    # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
    numeric(1)
  # 逐行说明：结束当前多行函数调用或表达式。
  )

  # 逐行说明：把本步骤的结果和诊断组织为命名列表。
  list(
    # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
    weight = weight,
    # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
    propensity_score = propensity_score,
    # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
    max_absolute_smd = absolute_smd[maximum_smd_index],
    # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
    max_absolute_smd_term = rownames(balance)[maximum_smd_index],
    # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
    ess_continuation = unname(ess[1]),
    # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
    ess_deescalation = unname(ess[2]),
    # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
    max_weight = max(weight),
    # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
    warning_n = length(captured_warnings)
  # 逐行说明：结束当前多行函数调用或表达式。
  )
# 逐行说明：结束当前函数、循环或条件分支。
}

# 逐行说明：定义函数 estimate_bootstrap_outcomes，封装本步骤的单一可复用逻辑。
estimate_bootstrap_outcomes <- function(data, weight) {
  # 逐行说明：将本行计算结果保存到对象 treatment，供后续校验或估计使用。
  treatment <- as.integer(as.character(data$deescalation))
  # 逐行说明：将本行计算结果保存到对象 event_factor，供后续校验或估计使用。
  event_factor <- make_aki_event_factor_boot(
    # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
    data$aki_7d,
    # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
    data$aki_competing_state
  # 逐行说明：结束当前多行函数调用或表达式。
  )
  # 逐行说明：将本行计算结果保存到对象 output，供后续校验或估计使用。
  output <- list()
  # 逐行说明：将本行计算结果保存到对象 manual_difference，供后续校验或估计使用。
  manual_difference <- 0
  # 逐行说明：按本行定义的索引逐项循环，确保每个试验、插补或策略均被处理。
  for (strategy in c(0L, 1L)) {
    # 逐行说明：将本行计算结果保存到对象 index，供后续校验或估计使用。
    index <- treatment == strategy
    # 逐行说明：检查本行条件；条件成立时执行后续保护或分支逻辑。
    if (!any(index)) stop("结局估计时缺少策略组。")
    # 逐行说明：将本行计算结果保存到对象 software，供后续校验或估计使用。
    software <- weighted_aj_boot(
      # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
      data$aki_followup_days[index],
      # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
      event_factor[index],
      # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
      weight[index]
    # 逐行说明：结束当前多行函数调用或表达式。
    )
    # 逐行说明：将本行计算结果保存到对象 manual，供后续校验或估计使用。
    manual <- weighted_aj_manual_boot(
      # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
      data$aki_followup_days[index],
      # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
      event_factor[index],
      # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
      weight[index]
    # 逐行说明：结束当前多行函数调用或表达式。
    )
    # 逐行说明：将本行计算结果保存到对象 manual_difference，供后续校验或估计使用。
    manual_difference <- max(
      # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
      manual_difference,
      # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
      abs(
        # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
        software[c(
          # 逐行说明：声明当前字符值、字段名或因子水平。
          "(s0)",
          # 逐行说明：声明当前字符值、字段名或因子水平。
          "aki",
          # 逐行说明：声明当前字符值、字段名或因子水平。
          "death_before_aki",
          # 逐行说明：声明当前字符值、字段名或因子水平。
          "discharge_before_aki"
        # 逐行说明：结束当前多行函数调用或表达式。
        )] -
          # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
          manual[c(
            # 逐行说明：声明当前字符值、字段名或因子水平。
            "(s0)",
            # 逐行说明：声明当前字符值、字段名或因子水平。
            "aki",
            # 逐行说明：声明当前字符值、字段名或因子水平。
            "death_before_aki",
            # 逐行说明：声明当前字符值、字段名或因子水平。
            "discharge_before_aki"
          # 逐行说明：结束当前多行函数调用或表达式。
          )]
      # 逐行说明：结束当前多行函数调用或表达式。
      )
    # 逐行说明：结束当前多行函数调用或表达式。
    )
    # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
    output[[paste0("risk_aki_", strategy)]] <- unname(software["aki"])
    # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
    output[[paste0("risk_death_before_aki_", strategy)]] <-
      # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
      unname(software["death_before_aki"])
    # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
    output[[paste0("risk_discharge_before_aki_", strategy)]] <-
      # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
      unname(software["discharge_before_aki"])
    # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
    output[[paste0("risk_death30_", strategy)]] <-
      # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
      weighted_binary_risk_boot(data$death_30d[index], weight[index])
  # 逐行说明：结束当前函数、循环或条件分支。
  }
  # 逐行说明：检查本行条件；条件成立时执行后续保护或分支逻辑。
  if (manual_difference > bootstrap_tolerance) {
    # 逐行说明：触发硬停止并报告原因，防止失败对象进入后续分析。
    stop("软件与手工Aalen–Johansen交叉核验失败。")
  # 逐行说明：结束当前函数、循环或条件分支。
  }
  # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
  output$aj_max_absolute_difference <- manual_difference
  # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
  output
# 逐行说明：结束当前函数、循环或条件分支。
}

# 逐行说明：定义函数 pool_bootstrap_imputations，封装本步骤的单一可复用逻辑。
pool_bootstrap_imputations <- function(imputation_rows) {
  # 逐行说明：将本行计算结果保存到对象 required，供后续校验或估计使用。
  required <- c(
    # 逐行说明：声明当前字符值、字段名或因子水平。
    "risk_aki_0",
    # 逐行说明：声明当前字符值、字段名或因子水平。
    "risk_aki_1",
    # 逐行说明：声明当前字符值、字段名或因子水平。
    "risk_death30_0",
    # 逐行说明：声明当前字符值、字段名或因子水平。
    "risk_death30_1"
  # 逐行说明：结束当前多行函数调用或表达式。
  )
  # 逐行说明：检查本行条件；条件成立时执行后续保护或分支逻辑。
  if (
    # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
    nrow(imputation_rows) != bootstrap_m ||
      # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
      anyNA(imputation_rows[, required]) ||
      # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
      any(!vapply(
        # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
        imputation_rows[, required],
        # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
        function(x) all(is.finite(x) & x >= 0 & x <= 1),
        # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
        logical(1)
      # 逐行说明：结束当前多行函数调用或表达式。
      ))
  # 逐行说明：结束当前多行函数调用或表达式。
  ) {
    # 逐行说明：触发硬停止并报告原因，防止失败对象进入后续分析。
    stop("20份插补风险合并输入非法。")
  # 逐行说明：结束当前函数、循环或条件分支。
  }
  # 逐行说明：将本行计算结果保存到对象 pooled，供后续校验或估计使用。
  pooled <- vapply(
    # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
    required,
    # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
    function(field) mean(imputation_rows[[field]]),
    # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
    numeric(1)
  # 逐行说明：结束当前多行函数调用或表达式。
  )
  # 逐行说明：检查本行条件；条件成立时执行后续保护或分支逻辑。
  if (pooled["risk_aki_0"] == 0 || pooled["risk_death30_0"] == 0) {
    # 逐行说明：触发硬停止并报告原因，防止失败对象进入后续分析。
    stop("继续组合并风险为0，风险比不可估计。")
  # 逐行说明：结束当前函数、循环或条件分支。
  }
  # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
  c(
    # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
    risk_aki_0 = unname(pooled["risk_aki_0"]),
    # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
    risk_aki_1 = unname(pooled["risk_aki_1"]),
    # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
    risk_death30_0 = unname(pooled["risk_death30_0"]),
    # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
    risk_death30_1 = unname(pooled["risk_death30_1"]),
    # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
    rd_aki = unname(pooled["risk_aki_1"] - pooled["risk_aki_0"]),
    # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
    rr_aki = unname(pooled["risk_aki_1"] / pooled["risk_aki_0"]),
    # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
    rd_death30 =
      # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
      unname(
        # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
        pooled["risk_death30_1"] - pooled["risk_death30_0"]
      # 逐行说明：结束当前多行函数调用或表达式。
      ),
    # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
    rr_death30 =
      # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
      unname(
        # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
        pooled["risk_death30_1"] / pooled["risk_death30_0"]
      # 逐行说明：结束当前多行函数调用或表达式。
      )
  # 逐行说明：结束当前多行函数调用或表达式。
  )
# 逐行说明：结束当前函数、循环或条件分支。
}

# 逐行说明：定义函数 run_single_bootstrap_repeat，封装本步骤的单一可复用逻辑。
run_single_bootstrap_repeat <- function(
  # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
  trial,
  # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
  trial_code,
  # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
  bootstrap_index,
  # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
  raw_data,
  # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
  dictionary,
  # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
  method_specification,
  # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
  predictor_specification,
  # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
  ps_formula
# 逐行说明：结束当前多行函数调用或表达式。
) {
  # 逐行说明：将本行计算结果保存到对象 sample_seed，供后续校验或估计使用。
  sample_seed <- bootstrap_master_seed +
    # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
    as.integer(trial_code) * 100000L +
    # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
    as.integer(bootstrap_index)
  # 逐行说明：将本行计算结果保存到对象 mice_seed，供后续校验或估计使用。
  mice_seed <- bootstrap_master_seed +
    # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
    as.integer(trial_code) * 100000L +
    # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
    50000L +
    # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
    as.integer(bootstrap_index)

  # 逐行说明：固定本步骤随机种子，保证抽样或插补可以复现。
  set.seed(sample_seed)
  # 逐行说明：将本行计算结果保存到对象 sampled_row，供后续校验或估计使用。
  sampled_row <- sample.int(
    # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
    nrow(raw_data),
    # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
    size = nrow(raw_data),
    # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
    replace = TRUE
  # 逐行说明：结束当前多行函数调用或表达式。
  )
  # 逐行说明：将本行计算结果保存到对象 bootstrap_raw，供后续校验或估计使用。
  bootstrap_raw <- raw_data[sampled_row, , drop = FALSE]
  # 逐行说明：将本行计算结果保存到对象 bootstrap_copy_id，供后续校验或估计使用。
  bootstrap_copy_id <- paste0(
    # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
    trial,
    # 逐行说明：声明当前字符值、字段名或因子水平。
    "_b",
    # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
    sprintf("%04d", bootstrap_index),
    # 逐行说明：声明当前字符值、字段名或因子水平。
    "_d",
    # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
    sprintf("%04d", seq_len(nrow(bootstrap_raw)))
  # 逐行说明：结束当前多行函数调用或表达式。
  )
  # 逐行说明：检查本行条件；条件成立时执行后续保护或分支逻辑。
  if (
    # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
    anyDuplicated(bootstrap_copy_id) > 0L ||
      # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
      length(bootstrap_copy_id) != nrow(bootstrap_raw)
  # 逐行说明：结束当前多行函数调用或表达式。
  ) {
    # 逐行说明：触发硬停止并报告原因，防止失败对象进入后续分析。
    stop("bootstrap_copy_id不唯一。")
  # 逐行说明：结束当前函数、循环或条件分支。
  }

  # 逐行说明：将本行计算结果保存到对象 work，供后续校验或估计使用。
  work <- build_bootstrap_working_data(
    # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
    bootstrap_raw,
    # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
    raw_data,
    # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
    trial,
    # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
    dictionary
  # 逐行说明：结束当前多行函数调用或表达式。
  )
  # 逐行说明：将本行计算结果保存到对象 method，供后续校验或估计使用。
  method <- rebuild_bootstrap_method(
    # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
    method_specification,
    # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
    trial,
    # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
    names(work)
  # 逐行说明：结束当前多行函数调用或表达式。
  )
  # anti-PSA source_control_status保持锁定logreg统计方法，仅绑定条件式重试包装器。
  if (identical(trial, "anti_psa")) {
    if (
      !identical(
        unname(method["source_control_status"]),
        "logreg"
      )
    ) {
      stop("anti_psa source_control_status锁定方法不是logreg。")
    }
    method["source_control_status"] <- "sourcecontrollogreg"
  }
  # 逐行说明：将本行计算结果保存到对象 predictor_matrix，供后续校验或估计使用。
  predictor_matrix <- rebuild_bootstrap_predictor_matrix(
    # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
    predictor_specification,
    # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
    trial,
    # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
    names(work)
  # 逐行说明：结束当前多行函数调用或表达式。
  )
  # 逐行说明：将本行计算结果保存到对象 urine_fields，供后续校验或估计使用。
  urine_fields <- c(
    # 逐行说明：声明当前字符值、字段名或因子水平。
    "urine_output_rate_6h_ml_kg_h",
    # 逐行说明：声明当前字符值、字段名或因子水平。
    "urine_output_rate_12h_ml_kg_h",
    # 逐行说明：声明当前字符值、字段名或因子水平。
    "urine_output_rate_24h_ml_kg_h"
  # 逐行说明：结束当前多行函数调用或表达式。
  )
  # 逐行说明：将本行计算结果保存到对象 active_targets，供后续校验或估计使用。
  active_targets <- names(method)[method != ""]
  # 逐行说明：检查本行条件；条件成立时执行后续保护或分支逻辑。
  if (!all(urine_fields %in% active_targets)) {
    # 逐行说明：触发硬停止并报告原因，防止失败对象进入后续分析。
    stop("三个尿量插补目标未全部处于活动状态。")
  # 逐行说明：结束当前函数、循环或条件分支。
  }
  # 逐行说明：将本行计算结果保存到对象 visit_sequence，供后续校验或估计使用。
  visit_sequence <- c(
    # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
    setdiff(active_targets, urine_fields),
    # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
    urine_fields
  # 逐行说明：结束当前多行函数调用或表达式。
  )

  # 记录MICE在返回mids对象时发出的loggedEvents汇总警告。
  captured_mice_warnings <- character(0)
  # 每个Bootstrap重复开始MICE前重置稀疏二分类内部拟合重试审计。
  if (
    !exists("reset_logreg100_retry_audit", mode = "function") ||
      !exists("get_logreg100_retry_audit", mode = "function")
  ) {
    stop("缺少logreg100条件式重试审计函数。")
  }
  reset_logreg100_retry_audit()
  # 执行锁定MICE；只静默“Number of logged events”汇总警告以保留事件明细。
  imputation <- withCallingHandlers(
    mice::mice(
      work,
      m = bootstrap_m,
      maxit = bootstrap_maxit,
      method = method,
      predictorMatrix = predictor_matrix,
      visitSequence = visit_sequence,
      donors = bootstrap_donors,
      seed = mice_seed,
      printFlag = FALSE,
      remove.collinear = FALSE,
      remove.constant = FALSE,
      polr.to.loggedEvents = TRUE
    ),
    warning = function(condition) {
      warning_message <- conditionMessage(condition)
      if (grepl(
        "^Number of logged events: [0-9]+$",
        warning_message
      )) {
        captured_mice_warnings <<- c(
          captured_mice_warnings,
          warning_message
        )
        invokeRestart("muffleWarning")
      }
    }
  )
  # 读取本次MICE中精确不收敛警告触发的条件式重试记录。
  mice_logreg_retry_summary <- get_logreg100_retry_audit()
  # 将无事件的NULL统一为固定列结构，便于后续检查点审计。
  mice_logged_events <- imputation$loggedEvents
  if (is.null(mice_logged_events)) {
    mice_logged_events <- data.frame(
      it = integer(0),
      im = integer(0),
      dep = character(0),
      meth = character(0),
      out = character(0),
      stringsAsFactors = FALSE
    )
  }
  # 第一类：anti_mrsa b=31验证的有序抗菌药.L/.Q对比列逐拟合剔除。
  allowed_ordered_contrast_event <- (
    mice_logged_events$dep == "urine_output_rate_24h_ml_kg_h" &
      mice_logged_events$meth == "pmm" &
      grepl(
        paste0(
          "^other_broad_spectrum_antibiotics\\.(L|Q)",
          "(, other_broad_spectrum_antibiotics\\.(L|Q))*$"
        ),
        mice_logged_events$out
      )
  )
  # 第二类：anti_mrsa b=33验证的、目标观测子集中计数确为0的竞争事件虚拟列。
  allowed_zero_count_death_event <- vapply(
    seq_len(nrow(mice_logged_events)),
    function(event_index) {
      dependent_field <- mice_logged_events$dep[event_index]
      if (
        !dependent_field %in% names(work) ||
          !dependent_field %in% names(method) ||
          mice_logged_events$out[event_index] !=
            "aki_competing_statedeath_before_aki" ||
          mice_logged_events$meth[event_index] !=
            unname(method[dependent_field])
      ) {
        return(FALSE)
      }
      observed_target <- !is.na(work[[dependent_field]])
      death_before_aki <- (
        as.character(work$aki_competing_state) == "death_before_aki"
      )
      sum(observed_target & death_before_aki, na.rm = TRUE) == 0L
    },
    logical(1)
  )
  # 合并两类方法不变且可从当前重抽样数据验证的内部剔除事件。
  allowed_lindep_event <- (
    allowed_ordered_contrast_event |
      allowed_zero_count_death_event
  )
  # 任一未知事件仍硬停止，禁止把未分类警告静默带入后续分析。
  if (
    imputation$m != bootstrap_m ||
      nrow(mice_logged_events) > 0L &&
        !all(allowed_lindep_event)
  ) {
    stop(
      "MICE未产生20份完成数据或出现未批准loggedEvents；",
      "事件数=",
      nrow(mice_logged_events)
    )
  }
  # 汇总已批准事件，随检查点永久保存。
  mice_logged_event_summary <- if (nrow(mice_logged_events) > 0L) {
    stats::aggregate(
      list(event_n = rep.int(1L, nrow(mice_logged_events))),
      by = mice_logged_events[c("dep", "meth", "out")],
      FUN = sum
    )
  } else {
    data.frame(
      dep = character(0),
      meth = character(0),
      out = character(0),
      event_n = integer(0),
      stringsAsFactors = FALSE
    )
  }

  # 逐行说明：将本行计算结果保存到对象 imputation_rows，供后续校验或估计使用。
  imputation_rows <- list()
  # 逐行说明：将本行计算结果保存到对象 diagnostic_rows，供后续校验或估计使用。
  diagnostic_rows <- list()
  # 逐行说明：按本行定义的索引逐项循环，确保每个试验、插补或策略均被处理。
  for (imputation_index in seq_len(bootstrap_m)) {
    # 逐行说明：将本行计算结果保存到对象 completed，供后续校验或估计使用。
    completed <- mice::complete(imputation, imputation_index)
    # 所有活动插补目标必须完整，不能只检查进入PS公式的变量。
    if (anyNA(completed[, active_targets, drop = FALSE])) {
      stop(
        "完成数据仍有活动插补目标缺失值；imputation=",
        imputation_index
      )
    }
    # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
    completed$bootstrap_copy_id <- bootstrap_copy_id
    # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
    completed$source_icu_stay_id <- bootstrap_raw$icu_stay_id
    # 逐行说明：检查本行条件；条件成立时执行后续保护或分支逻辑。
    if (
      # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
      nrow(completed) != nrow(bootstrap_raw) ||
        # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
        anyDuplicated(completed$bootstrap_copy_id) > 0L ||
        # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
        !identical(
          # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
          as.character(completed$bootstrap_copy_id),
          # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
          bootstrap_copy_id
        # 逐行说明：结束当前多行函数调用或表达式。
        )
    # 逐行说明：结束当前多行函数调用或表达式。
    ) {
      # 逐行说明：触发硬停止并报告原因，防止失败对象进入后续分析。
      stop("完成数据的复制键连接失败。")
    # 逐行说明：结束当前函数、循环或条件分支。
    }

    # 逐行说明：将本行计算结果保存到对象 model_data，供后续校验或估计使用。
    model_data <- completed
    # 逐行说明：检查本行条件；条件成立时执行后续保护或分支逻辑。
    if (
      # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
      anyNA(model_data[, all.vars(ps_formula), drop = FALSE])
    # 逐行说明：结束当前多行函数调用或表达式。
    ) {
      # 逐行说明：触发硬停止并报告原因，防止失败对象进入后续分析。
      stop("完成数据仍有PS模型缺失值。")
    # 逐行说明：结束当前函数、循环或条件分支。
    }
    # 逐行说明：按本行定义的索引逐项循环，确保每个试验、插补或策略均被处理。
    for (model in c("cbps_ate", "ato")) {
      # 逐行说明：将本行计算结果保存到对象 weight_fit，供后续校验或估计使用。
      weight_fit <- fit_bootstrap_weights(
        # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
        model_data,
        # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
        ps_formula,
        # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
        model
      # 逐行说明：结束当前多行函数调用或表达式。
      )
      # 逐行说明：将本行计算结果保存到对象 outcome_fit，供后续校验或估计使用。
      outcome_fit <- estimate_bootstrap_outcomes(
        # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
        model_data,
        # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
        weight_fit$weight
      # 逐行说明：结束当前多行函数调用或表达式。
      )
      # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
      imputation_rows[[length(imputation_rows) + 1L]] <- data.frame(
        # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
        trial = trial,
        # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
        bootstrap_index = bootstrap_index,
        # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
        imputation = imputation_index,
        # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
        model = model,
        # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
        risk_aki_0 = outcome_fit$risk_aki_0,
        # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
        risk_aki_1 = outcome_fit$risk_aki_1,
        # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
        risk_death30_0 = outcome_fit$risk_death30_0,
        # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
        risk_death30_1 = outcome_fit$risk_death30_1,
        # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
        aj_max_absolute_difference =
          # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
          outcome_fit$aj_max_absolute_difference,
        # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
        stringsAsFactors = FALSE
      # 逐行说明：结束当前多行函数调用或表达式。
      )
      # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
      diagnostic_rows[[length(diagnostic_rows) + 1L]] <- data.frame(
        # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
        trial = trial,
        # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
        bootstrap_index = bootstrap_index,
        # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
        imputation = imputation_index,
        # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
        model = model,
        # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
        max_absolute_smd = weight_fit$max_absolute_smd,
        # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
        max_absolute_smd_term = weight_fit$max_absolute_smd_term,
        # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
        ess_continuation = weight_fit$ess_continuation,
        # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
        ess_deescalation = weight_fit$ess_deescalation,
        # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
        max_weight = weight_fit$max_weight,
        # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
        warning_n = weight_fit$warning_n,
        ato_brglm_retry_used = weight_fit$ato_brglm_retry_used,
        ato_brglm_initial_warning_n =
          weight_fit$ato_brglm_initial_warning_n,
        ato_brglm_retry_maxit = weight_fit$ato_brglm_retry_maxit,
        # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
        stringsAsFactors = FALSE
      # 逐行说明：结束当前多行函数调用或表达式。
      )
    # 逐行说明：结束当前函数、循环或条件分支。
    }
  # 逐行说明：结束当前函数、循环或条件分支。
  }
  # 逐行说明：将本行计算结果保存到对象 imputation_rows，供后续校验或估计使用。
  imputation_rows <- do.call(rbind, imputation_rows)
  # 逐行说明：将本行计算结果保存到对象 diagnostic_rows，供后续校验或估计使用。
  diagnostic_rows <- do.call(rbind, diagnostic_rows)
  # 逐行说明：将本行计算结果保存到对象 pooled_rows，供后续校验或估计使用。
  pooled_rows <- list()
  # 逐行说明：按本行定义的索引逐项循环，确保每个试验、插补或策略均被处理。
  for (model in c("cbps_ate", "ato")) {
    # 逐行说明：将本行计算结果保存到对象 current，供后续校验或估计使用。
    current <- imputation_rows[imputation_rows$model == model, ]
    # 逐行说明：将本行计算结果保存到对象 pooled，供后续校验或估计使用。
    pooled <- pool_bootstrap_imputations(current)
    # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
    pooled_rows[[length(pooled_rows) + 1L]] <- data.frame(
      # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
      trial = trial,
      # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
      bootstrap_index = bootstrap_index,
      # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
      model = model,
      # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
      risk_aki_0 = unname(pooled["risk_aki_0"]),
      # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
      risk_aki_1 = unname(pooled["risk_aki_1"]),
      # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
      rd_aki = unname(pooled["rd_aki"]),
      # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
      rr_aki = unname(pooled["rr_aki"]),
      # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
      risk_death30_0 = unname(pooled["risk_death30_0"]),
      # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
      risk_death30_1 = unname(pooled["risk_death30_1"]),
      # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
      rd_death30 =
        # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
        unname(pooled["rd_death30"]),
      # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
      rr_death30 =
        # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
        unname(pooled["rr_death30"]),
      # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
      stringsAsFactors = FALSE
    # 逐行说明：结束当前多行函数调用或表达式。
    )
  # 逐行说明：结束当前函数、循环或条件分支。
  }

  # 逐行说明：把本步骤的结果和诊断组织为命名列表。
  list(
    # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
    status = "success",
    # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
    trial = trial,
    # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
    bootstrap_index = bootstrap_index,
    # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
    sample_seed = sample_seed,
    # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
    mice_seed = mice_seed,
    # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
    sampled_row = sampled_row,
    # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
    bootstrap_copy_id = bootstrap_copy_id,
    # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
    source_icu_stay_id = bootstrap_raw$icu_stay_id,
    # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
    sampled_unique_source_n = length(unique(bootstrap_raw$icu_stay_id)),
    # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
    sampled_strategy_0_n = sum(bootstrap_raw$deescalation == 0),
    # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
    sampled_strategy_1_n = sum(bootstrap_raw$deescalation == 1),
    # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
    mice_logged_event_n = nrow(mice_logged_events),
    # 保存逐类事件摘要，支持批次级自动审计。
    mice_logged_event_summary = mice_logged_event_summary,
    # 保存MICE汇总警告正文，不包含患者级数据。
    mice_warning_messages = captured_mice_warnings,
    # 保存稀疏二分类内部glm.fit条件式重试次数。
    mice_logreg_retry_n = nrow(mice_logreg_retry_summary),
    # 保存条件式重试的调用序号和迭代收敛审计。
    mice_logreg_retry_summary = mice_logreg_retry_summary,
    # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
    imputation_estimates = imputation_rows,
    # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
    weight_diagnostics = diagnostic_rows,
    # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
    pooled_estimates = do.call(rbind, pooled_rows),
    # 逐行说明：执行本行表达式；其输入和输出由相邻变量名及上方方法块限定。
    imputation = imputation
  # 逐行说明：结束当前多行函数调用或表达式。
  )
# 逐行说明：结束当前函数、循环或条件分支。
}
