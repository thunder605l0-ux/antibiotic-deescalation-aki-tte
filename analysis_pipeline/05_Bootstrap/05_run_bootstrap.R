# ==============================================================================
# Analysis: 两个平行试验的可恢复正式Boot-MI
# Date: 2026-07-29
# Random seed: 20260727及确定性派生种子
# R: 4.4.2
#
# Method:
# Purpose: 从b=1开始运行患者级非参数Bootstrap，并在每次重复中重新MICE和加权。
# Source type: project_locked_code / R_official_documentation / self_written
# Source name or URL:
#   04_bootstrap_single_repeat.R
#   https://stat.ethz.ch/R-manual/R-devel/library/base/html/conditions.html
# Citation or commit/version:
#   D-4.7N-ANALYSIS-LOCK-V2.0
#   D-4.7K-BOOT-MI-BRLOGIT-AMENDMENT-V1.1
# Access date: 2026-07-29
# Adaptation notes:
#   不复制旧检查点；新目录从b=1运行；成功任务支持安全跳过和续跑；
#   任何警告或错误均保存任务、种子、阶段和调用摘要后停止。
# Verification command:
#   Rscript --vanilla 05_Bootstrap/05_run_bootstrap.R --trial=anti_mrsa --b_start=1 --b_end=10
# ==============================================================================

# 定义读取--name=value命令行参数的函数。
get_argument <- function(name, default) {
  # 读取传给当前脚本的尾随参数。
  arguments <- commandArgs(trailingOnly = TRUE)
  # 构造当前参数的固定前缀。
  prefix <- paste0("--", name, "=")
  # 选择匹配当前前缀的参数。
  matched <- arguments[startsWith(arguments, prefix)]
  # 未提供参数时返回默认值。
  if (length(matched) == 0L) {
    # 返回调用者指定的默认值。
    return(default)
  }
  # 同一个参数重复出现时停止。
  if (length(matched) != 1L) {
    # 报告重复的参数名。
    stop("参数重复：", name)
  }
  # 删除参数名前缀，仅返回参数值。
  sub(prefix, "", matched, fixed = TRUE)
}

# 读取完整命令行以定位当前脚本。
# 使用用户锁定的绝对路径定位br.logit单次Boot-MI引擎。
engine_path <- normalizePath(
  file.path(Sys.getenv("TTE_RUNTIME_ROOT", unset = "."), "05_Bootstrap", "04_bootstrap_single_repeat.R"),
  winslash = "/",
  mustWork = TRUE
)
# 加载统一配置、锁定核心和br.logit适配函数。
source(engine_path, local = globalenv(), encoding = "UTF-8")
# 将工作目录固定到运行包根目录，保证交互式与命令行执行一致。
setwd(package_root)
# 核对任务运行需要的包。
require_packages(c("mice", "WeightIt", "brglm2", "cobalt", "survival"))

# 读取试验参数；默认按注册表顺序运行两个试验。
trial_argument <- get_argument("trial", "all")
# 读取Bootstrap起始编号；默认从b=1开始。
b_start <- as.integer(get_argument("b_start", "1"))
# 读取Bootstrap结束编号；默认运行至锁定B=1000。
b_end <- as.integer(get_argument("b_end", as.character(bootstrap_B)))
# 检查起止编号均为有效范围内整数。
if (
  is.na(b_start) ||
    is.na(b_end) ||
    b_start < 1L ||
    b_end > bootstrap_B ||
    b_start > b_end
) {
  # 报告允许的编号范围。
  stop("b_start和b_end必须满足1≤b_start≤b_end≤1000。")
}
# 根据trial参数选择一个或两个平行试验。
selected_trials <- if (identical(trial_argument, "all")) {
  # 返回两个锁定试验。
  trial_registry$trial
} else if (trial_argument %in% trial_registry$trial) {
  # 返回用户明确指定的单个锁定试验。
  trial_argument
} else {
  # 拒绝未知试验名称。
  stop("trial必须是anti_mrsa、anti_psa或all。")
}

# 读取MICE工作变量字典。
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
# 读取锁定定向尿量预测矩阵。
predictor_specification <- utils::read.csv(
  file.path(config_dir, "mice_predictor_matrix_directed_uo.csv"),
  stringsAsFactors = FALSE
)
# 声明两个试验的输入文件。
input_files <- c(
  anti_mrsa = "analysis_dataset_anti_mrsa.csv",
  anti_psa = "analysis_dataset_anti_psa.csv"
)

# 构造Bootstrap输出根目录。
bootstrap_output <- file.path(output_root, "04_Bootstrap")
# 构造成功检查点目录。
checkpoint_dir <- file.path(bootstrap_output, "checkpoints")
# 构造任务状态目录。
status_dir <- file.path(bootstrap_output, "status")
# 构造失败对象目录。
failure_dir <- file.path(bootstrap_output, "failures")
# 确保三个输出目录存在。
dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
# 确保状态目录存在。
dir.create(status_dir, recursive = TRUE, showWarnings = FALSE)
# 确保失败目录存在。
dir.create(failure_dir, recursive = TRUE, showWarnings = FALSE)

# 逐个选定试验执行连续Bootstrap任务。
for (trial in selected_trials) {
  # 读取当前试验注册行。
  trial_row <- trial_registry[
    trial_registry$trial == trial,
    ,
    drop = FALSE
  ]
  # 构造当前试验原始分析CSV路径。
  raw_path <- file.path(data_dir, input_files[[trial]])
  # 读取当前试验原始不完整患者—试验数据。
  raw_data <- utils::read.csv(
    raw_path,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    na.strings = c("", "NA")
  )
  # 核对当前试验人数和暴露组人数。
  if (
    nrow(raw_data) != trial_row$expected_n ||
      sum(raw_data$deescalation == 1L) !=
        trial_row$expected_deescalation_n ||
      sum(raw_data$deescalation == 0L) !=
        trial_row$expected_continuation_n
  ) {
    # 队列身份不一致时停止。
    stop(trial, "原始数据人数或策略人数与锁定值不一致。")
  }
  # 按用户指定范围逐个运行Bootstrap重复。
  for (bootstrap_index in seq.int(b_start, b_end)) {
    # 构造稳定任务标识。
    task_id <- paste0(
      trial,
      "_b",
      sprintf("%04d", bootstrap_index)
    )
    # 构造当前任务成功检查点路径。
    checkpoint_path <- file.path(
      checkpoint_dir,
      paste0(task_id, "_checkpoint.rds")
    )
    # 构造当前任务状态CSV路径。
    status_path <- file.path(
      status_dir,
      paste0(task_id, "_status.csv")
    )
    # 构造当前任务失败对象路径。
    failure_path <- file.path(
      failure_dir,
      paste0(task_id, "_failure.rds")
    )
    # 已存在成功检查点时先验证身份再安全跳过。
    if (file.exists(checkpoint_path)) {
      # 读取已有检查点。
      existing <- readRDS(checkpoint_path)
      # 验证检查点属于当前试验、重复且状态成功。
      checkpoint_valid <- (
        identical(existing$status, "success") &&
          identical(existing$trial, trial) &&
          identical(
            as.integer(existing$bootstrap_index),
            as.integer(bootstrap_index)
          )
      )
      # 检查点身份异常时停止，不覆盖可疑对象。
      if (!checkpoint_valid) {
        # 报告异常任务标识。
        stop(task_id, "已有检查点身份异常。")
      }
      # 输出安全跳过信息。
      message(task_id, "已有有效检查点，跳过。")
      # 进入下一个Bootstrap任务。
      next
    }
    # 记录任务开始时间。
    started_at <- format(Sys.time(), "%Y-%m-%d %H:%M:%S%z")
    # 初始化外层警告正文。
    outer_warning <- character(0)
    # 在统一警告处理器和错误捕获器中运行单次完整Boot-MI。
    result <- tryCatch(
      # 把任何未被内部处理的警告升级为可持久化错误。
      withCallingHandlers(
        # 调用锁定的单次Boot-MI函数。
        run_single_bootstrap_repeat(
          trial = trial,
          trial_code = trial_row$trial_code,
          bootstrap_index = bootstrap_index,
          raw_data = raw_data,
          dictionary = dictionary,
          method_specification = method_specification,
          predictor_specification = predictor_specification,
          ps_formula = ps_formula_main
        ),
        # 定义全任务外层警告处理函数。
        warning = function(condition) {
          # 保存当前警告正文。
          outer_warning <<- c(
            outer_warning,
            conditionMessage(condition)
          )
          # 把警告立即升级为错误，使外层tryCatch进入失败分支。
          stop(
            "OUTER_WARNING_CAPTURED[",
            task_id,
            "]: ",
            conditionMessage(condition),
            call. = FALSE
          )
        }
      ),
      # 返回错误条件对象而不是丢失上下文。
      error = function(condition) condition
    )
    # 判断单次任务是否失败。
    if (inherits(result, "error")) {
      # 构造包含任务、种子、阶段和调用摘要的失败对象。
      failure_object <- list(
        task_id = task_id,
        trial = trial,
        bootstrap_index = bootstrap_index,
        stage = "run_single_bootstrap_repeat",
        sample_seed =
          master_seed + trial_row$trial_code * 100000L + bootstrap_index,
        mice_seed =
          master_seed + trial_row$trial_code * 100000L + 50000L +
            bootstrap_index,
        error_message = conditionMessage(result),
        outer_warnings = outer_warning,
        call = deparse(conditionCall(result)),
        started_at = started_at,
        failed_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S%z")
      )
      # 原子保存失败对象，确保警告和错误不会仅留在控制台。
      atomic_save_rds(failure_object, failure_path)
      # 构造失败状态表。
      failure_status <- data.frame(
        task_id = task_id,
        trial = trial,
        bootstrap_index = bootstrap_index,
        status = "failed",
        error_message = conditionMessage(result),
        failure_file = file.path(
          "09_输出结果",
          "04_Bootstrap",
          "failures",
          basename(failure_path)
        ),
        stringsAsFactors = FALSE
      )
      # 原子写出失败状态。
      atomic_write_csv(failure_status, status_path)
      # 保存失败证据后停止，不跳过、不换种子。
      stop(
        task_id,
        "失败并已保存证据：",
        conditionMessage(result)
      )
    }
    # 为成功对象补充明确的当前引擎身份。
    result$engine <- "boot_mi_brlogit_independent_package_v1"
    # 为成功对象补充Analysis Lock身份。
    result$analysis_lock <- "D-4.7N-ANALYSIS-LOCK-V2.0"
    # 为成功对象登记完成时间。
    result$completed_at <- format(Sys.time(), "%Y-%m-%d %H:%M:%S%z")
    # 原子保存成功检查点。
    atomic_save_rds(result, checkpoint_path)
    # 构造成功状态表。
    success_status <- data.frame(
      task_id = task_id,
      trial = trial,
      bootstrap_index = bootstrap_index,
      status = "success",
      sample_seed = result$sample_seed,
      mice_seed = result$mice_seed,
      checkpoint_file = file.path(
        "09_输出结果",
        "04_Bootstrap",
        "checkpoints",
        basename(checkpoint_path)
      ),
      stringsAsFactors = FALSE
    )
    # 原子写出成功状态。
    atomic_write_csv(success_status, status_path)
    # 向控制台报告当前任务完成。
    message(task_id, "完成。")
  }
}

# 报告本次请求范围全部完成。
message(
  "本次Bootstrap运行完成：trial=",
  trial_argument,
  "，b=",
  b_start,
  "–",
  b_end,
  "。"
)
