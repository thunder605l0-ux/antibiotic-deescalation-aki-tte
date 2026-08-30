# ==============================================================================
# Analysis: 两个平行试验的锁定多重插补
# Date: 2026-07-29
# Random seed: 20260727
# R: 4.4.2
#
# Method:
# Purpose: 分别执行NEW_MS01的m=20、maxit=50定向尿量MICE。
# Source type: GitHub / local_knowledge_base / project_locked_protocol
# Source name or URL:
#   https://github.com/amices/mice
#   D:/source/obsidian
#   analysis_plan_locked.md
# Citation or commit/version:
#   amices/mice commit 61083e667fa41cc4bc655492591b100bf8a7e443
#   mice 3.19.10
#   D-4.7N-ANALYSIS-LOCK-V2.0
# Access date: 2026-07-29
# Adaptation notes:
#   尿量访问顺序固定为6h、12h、24h；PMM donors=5；
#   暴露和结局只作预测变量而不被插补；原始CSV保持只读。
# Verification command:
#   Rscript --vanilla 03_多重插补/02_run_multiple_imputation.R
# ==============================================================================

# 读取完整命令行以定位当前脚本。
# 使用用户锁定的绝对路径加载统一配置，不依赖当前工作目录。
source(
  file.path(Sys.getenv("TTE_RUNTIME_ROOT", unset = "."), "02_配置与变量字典", "analysis_config.R"),
  local = globalenv(),
  encoding = "UTF-8"
)
# 将工作目录固定到运行包根目录，保证交互式与命令行执行一致。
setwd(package_root)
# 核对MICE和哈希计算依赖。
require_packages(c("mice", "digest"))
# 固定检查mice必须从复制的GitHub隔离库加载。
expected_mice_version <- "3.19.10"
# 读取实际加载的mice版本。
actual_mice_version <- as.character(utils::packageVersion("mice"))
# 版本不符合锁定值时停止。
if (!identical(actual_mice_version, expected_mice_version)) {
  # 报告实际版本，便于研究者修复软件环境。
  stop(
    "mice版本不符合锁定值；实际=",
    actual_mice_version,
    "，期望=",
    expected_mice_version
  )
}

# 构造Boot-MI核心函数路径，以复用同一数据重编码和MICE矩阵逻辑。
core_path <- file.path(
  package_root,
  "05_Bootstrap",
  "04_bootstrap_core_locked.R"
)
# 加载经过逐行注释的锁定核心函数。
source(core_path, local = globalenv(), encoding = "UTF-8")
source(
  file.path(config_dir, "mice_imputation_extensions.R"),
  local = globalenv(),
  encoding = "UTF-8"
)
# 用统一配置覆盖核心文件中的MICE常量。
bootstrap_m <- mice_m
# 用统一配置覆盖核心文件中的迭代次数。
bootstrap_maxit <- mice_maxit
# 用统一配置覆盖核心文件中的PMM供体数。
bootstrap_donors <- mice_donors

# 读取工作变量字典。
dictionary <- utils::read.csv(
  file.path(config_dir, "mice_working_variable_dictionary.csv"),
  stringsAsFactors = FALSE,
  na.strings = ""
)
# 读取锁定MICE方法规格。
method_specification <- utils::read.csv(
  file.path(config_dir, "mice_method_specification.csv"),
  stringsAsFactors = FALSE,
  na.strings = ""
)
# 读取定向尿量预测矩阵的长表规格。
predictor_specification <- utils::read.csv(
  file.path(config_dir, "mice_predictor_matrix_directed_uo.csv"),
  stringsAsFactors = FALSE
)

# 声明两个试验的输入文件。
input_files <- c(
  anti_mrsa = "analysis_dataset_anti_mrsa.csv",
  anti_psa = "analysis_dataset_anti_psa.csv"
)
# 声明尿量变量的固定单向访问顺序。
urine_fields <- c(
  "urine_output_rate_6h_ml_kg_h",
  "urine_output_rate_12h_ml_kg_h",
  "urine_output_rate_24h_ml_kg_h"
)
# 初始化MICE运行摘要列表。
summary_rows <- list()

# 逐个平行试验执行独立插补。
for (trial_index in seq_len(nrow(trial_registry))) {
  # 读取当前试验注册信息。
  trial_row <- trial_registry[trial_index, , drop = FALSE]
  # 保存当前试验标识。
  trial <- trial_row$trial
  # 构造当前试验原始分析CSV路径。
  raw_path <- file.path(data_dir, input_files[[trial]])
  # 读取当前试验的原始不完整患者—试验数据。
  raw_data <- utils::read.csv(
    raw_path,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    na.strings = c("", "NA")
  )
  # 核对当前试验人数与锁定值一致。
  if (nrow(raw_data) != trial_row$expected_n) {
    # 人数漂移时停止，防止对错误队列插补。
    stop(trial, "原始数据人数与锁定值不一致。")
  }
  # 使用与正式Bootstrap相同的函数构建MICE工作数据。
  work_data <- build_bootstrap_working_data(
    raw_data,
    raw_data,
    trial,
    dictionary
  )
  # 根据当前试验和工作字段重建锁定插补方法向量。
  method <- rebuild_bootstrap_method(
    method_specification,
    trial,
    names(work_data)
  )
  # 根据当前试验和工作字段重建锁定预测矩阵。
  predictor_matrix <- rebuild_bootstrap_predictor_matrix(
    predictor_specification,
    trial,
    names(work_data)
  )
  # 提取所有需要主动插补的字段。
  active_targets <- names(method)[method != ""]
  # 确保6h、12h和24h尿量均处于主动插补状态。
  if (!all(urine_fields %in% active_targets)) {
    # 缺少任一尿量插补目标时停止。
    stop(trial, "三个尿量插补目标未全部启用。")
  }
  # 把其他插补变量放在尿量变量之前。
  visit_sequence <- c(
    setdiff(active_targets, urine_fields),
    urine_fields
  )
  # 原始样本MICE严格使用Analysis Lock登记的固定种子。
  # 只有Bootstrap内部才按试验代码和重复编号派生种子；
  # 原始样本若按试验派生种子将无法复现已锁定MIDS对象。
  trial_mice_seed <- master_seed
  # 执行20份、50次迭代的定向尿量MICE。
  imputation <- mice::mice(
    work_data,
    m = mice_m,
    maxit = mice_maxit,
    method = method,
    predictorMatrix = predictor_matrix,
    visitSequence = visit_sequence,
    donors = mice_donors,
    seed = trial_mice_seed,
    printFlag = FALSE,
    remove.collinear = FALSE,
    remove.constant = FALSE,
    polr.to.loggedEvents = TRUE
  )
  # 计算MICE记录的事件数；无事件时记为0。
  logged_event_n <- if (is.null(imputation$loggedEvents)) {
    # 返回0表示没有logged events。
    0L
  } else {
    # 返回实际logged events行数。
    nrow(imputation$loggedEvents)
  }
  # 要求MICE确实产生20份插补数据且无未解决事件。
  if (imputation$m != mice_m || logged_event_n > 0L) {
    # 插补合同失败时停止并报告试验。
    stop(
      trial,
      "MICE合同失败；m=",
      imputation$m,
      "，logged events=",
      logged_event_n
    )
  }
  # 初始化残余缺失计数。
  residual_missing_n <- 0L
  # 逐份完成数据核对主PS字段不存在残余缺失。
  for (imputation_index in seq_len(mice_m)) {
    # 提取当前完成数据。
    completed <- mice::complete(imputation, imputation_index)
    # 累加主PS公式字段中的残余缺失数。
    residual_missing_n <- residual_missing_n + sum(
      is.na(completed[, all.vars(ps_formula_main), drop = FALSE])
    )
    # 核对被动SOFA恒等式在当前完成数据中成立。
    passive_error <- max(abs(
      completed$nonrenal_sofa_change -
        (
          completed$recent_nonrenal_sofa -
            completed$early_nonrenal_sofa
        )
    ))
    # 被动恒等式误差超过数值容差时停止。
    if (!is.finite(passive_error) || passive_error > numeric_tolerance) {
      # 报告具体试验和插补编号。
      stop(trial, "第", imputation_index, "份插补的SOFA恒等式失败。")
    }
  }
  # 任何主PS字段残余缺失均不允许进入权重模型。
  if (residual_missing_n > 0L) {
    # 报告残余缺失单元格数量。
    stop(trial, "完成数据仍有主PS字段缺失：", residual_missing_n)
  }
  # 把原始数据、工作数据、MICE对象和配置共同保存，便于后续追溯。
  saved_object <- list(
    trial = trial,
    raw_data = raw_data,
    work_data = work_data,
    imputation = imputation,
    method = method,
    predictor_matrix = predictor_matrix,
    visit_sequence = visit_sequence,
    seed = trial_mice_seed,
    mice_version = actual_mice_version,
    ps_formula_main = ps_formula_main
  )
  # 构造当前试验MICE对象输出路径。
  object_path <- file.path(
    output_root,
    "02_MICE",
    paste0("mids_", trial, "_m20_maxit50.rds")
  )
  # 原子保存当前试验MICE对象。
  atomic_save_rds(saved_object, object_path)
  # 把当前试验MICE运行摘要追加到列表。
  summary_rows[[length(summary_rows) + 1L]] <- data.frame(
    trial = trial,
    n = nrow(raw_data),
    m = imputation$m,
    maxit = mice_maxit,
    donors = mice_donors,
    seed = trial_mice_seed,
    logged_event_n = logged_event_n,
    residual_missing_n = residual_missing_n,
    object_path = file.path(
      "09_输出结果",
      "02_MICE",
      basename(object_path)
    ),
    stringsAsFactors = FALSE
  )
}

# 合并两个试验的插补摘要。
summary_table <- do.call(rbind, summary_rows)
# 构造MICE摘要CSV路径。
summary_path <- file.path(
  output_root,
  "02_MICE",
  "mice_run_summary.csv"
)
# 原子写出MICE运行摘要。
atomic_write_csv(summary_table, summary_path)
# 向控制台报告MICE完成状态。
message("MICE完成：两个试验均为m=20、maxit=50，logged events=0。")
