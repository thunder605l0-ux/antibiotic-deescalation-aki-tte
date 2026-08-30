# ==============================================================================
# Analysis: 独立运行包原始分析数据结构核验
# Date: 2026-07-29
# Random seed: not applicable
# R: 4.4.2
#
# Method:
# Purpose: 在任何插补或模型拟合前核对数据身份、键、人数、暴露和结局编码。
# Source type: local_knowledge_base / project_locked_protocol / self_written
# Source name or URL:
#   analysis_plan_locked.md
# Citation or commit/version: D-4.7N-ANALYSIS-LOCK-V2.0
# Access date: 2026-07-29
# Adaptation notes: 只读检查两份52列CSV，不修改原始数据。
# Verification command:
#   Rscript --vanilla 00_运行说明/01_validate_input_data.R
# ==============================================================================

# 读取完整命令行以定位本脚本。
# 使用用户锁定的绝对路径加载统一配置，不依赖当前工作目录。
source(
  file.path(Sys.getenv("TTE_RUNTIME_ROOT", unset = "."), "02_配置与变量字典", "analysis_config.R"),
  local = globalenv(),
  encoding = "UTF-8"
)
# 将工作目录固定到运行包根目录，保证交互式与命令行执行一致。
setwd(package_root)
# 核对哈希计算所需的digest包。
require_packages("digest")

# 声明两个只读输入文件及其锁定SHA-256。
input_registry <- data.frame(
  # 使用试验标识连接trial_registry。
  trial = c("anti_mrsa", "anti_psa"),
  # 登记运行包内相对文件名。
  file = c(
    "analysis_dataset_anti_mrsa.csv",
    "analysis_dataset_anti_psa.csv"
  ),
  # 登记复制前原项目最终输入的SHA-256。
  expected_sha256 = c(
    "F04D40E660EC8F1FE849FD7F433CD66DD6BE64E01C3772E25F50D112F68E5DD6",
    "A366F8928BEE8142FD34C0B0DAF85EDCCF8A6F0669F79EC4082E535FE22E034F"
  ),
  # 保持字符列为字符类型。
  stringsAsFactors = FALSE
)

# 声明必须存在的关键字段。
required_fields <- c(
  "patient_id",
  "hospital_admission_id",
  "icu_stay_id",
  "target_coverage_trial",
  "deescalation",
  "index_time",
  "aki_7d",
  "aki_followup_days",
  "aki_7d_event_code",
  "death_30d"
)

# 初始化逐试验QA结果列表。
qa_rows <- list()
# 逐个平行试验核验输入数据。
for (row_index in seq_len(nrow(input_registry))) {
  # 读取当前输入注册行。
  registry_row <- input_registry[row_index, , drop = FALSE]
  # 根据试验标识取得锁定人数。
  trial_row <- trial_registry[
    trial_registry$trial == registry_row$trial,
    ,
    drop = FALSE
  ]
  # 构造当前CSV的绝对路径。
  data_path <- file.path(data_dir, registry_row$file)
  # 输入文件缺失时立即停止。
  if (!file.exists(data_path)) {
    # 报告具体缺失路径。
    stop("缺少原始分析数据：", data_path)
  }
  # 计算当前CSV的SHA-256并转换为大写。
  actual_sha256 <- toupper(digest::digest(
    file = data_path,
    algo = "sha256",
    serialize = FALSE
  ))
  # 读取CSV但不转换字段名或字符列。
  data <- utils::read.csv(
    data_path,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    na.strings = c("", "NA")
  )
  # 检查关键字段是否全部存在。
  fields_ok <- all(required_fields %in% names(data))
  # 字段缺失时报告完整差集并停止。
  if (!fields_ok) {
    # 输出所有缺失字段名称。
    stop(
      registry_row$trial,
      "缺少字段：",
      paste(setdiff(required_fields, names(data)), collapse = ", ")
    )
  }
  # 检查患者—试验主键icu_stay_id是否唯一。
  key_unique <- anyDuplicated(data$icu_stay_id) == 0L
  # 检查试验字段是否只包含当前试验。
  observed_trial_values <- unique(data$target_coverage_trial)
  # 要求试验字段只有一个值且与注册表数字试验编码一致。
  trial_value_ok <- (
    length(observed_trial_values) == 1L &&
      identical(
        as.integer(observed_trial_values),
        as.integer(trial_row$trial_code)
      )
  )
  # 检查暴露是否完整且仅含0和1。
  exposure_ok <- (
    !anyNA(data$deescalation) &&
      setequal(unique(data$deescalation), c(0L, 1L))
  )
  # 检查主要AKI结局是否完整且仅含0和1。
  aki_ok <- (
    !anyNA(data$aki_7d) &&
      all(data$aki_7d %in% c(0L, 1L))
  )
  # 检查30天死亡结局是否完整且仅含0和1。
  death_ok <- (
    !anyNA(data$death_30d) &&
      all(data$death_30d %in% c(0L, 1L))
  )
  # 把当前试验的全部核验结果组织为一行。
  qa_rows[[length(qa_rows) + 1L]] <- data.frame(
    trial = registry_row$trial,
    file = registry_row$file,
    n = nrow(data),
    p = ncol(data),
    deescalation_n = sum(data$deescalation == 1L),
    continuation_n = sum(data$deescalation == 0L),
    sha256 = actual_sha256,
    sha256_ok = identical(actual_sha256, registry_row$expected_sha256),
    n_ok = nrow(data) == trial_row$expected_n,
    p_ok = ncol(data) == 52L,
    key_unique = key_unique,
    trial_value_ok = trial_value_ok,
    exposure_ok = exposure_ok,
    aki_ok = aki_ok,
    death_ok = death_ok,
    stringsAsFactors = FALSE
  )
}

# 合并两个试验的QA行。
qa_table <- do.call(rbind, qa_rows)
# 检查所有逻辑型QA字段是否全部通过。
all_passed <- all(unlist(qa_table[
  ,
  c(
    "sha256_ok",
    "n_ok",
    "p_ok",
    "key_unique",
    "trial_value_ok",
    "exposure_ok",
    "aki_ok",
    "death_ok"
  )
]))
# 构造输入QA输出路径。
qa_path <- file.path(
  output_root,
  "01_输入QA",
  "input_data_qa.csv"
)
# 原子写出输入QA表。
atomic_write_csv(qa_table, qa_path)
# 任一QA失败时停止并保留已写出的审计表。
if (!all_passed) {
  # 报告QA文件路径供研究者定位失败项。
  stop("输入数据QA未全部通过；详见：", qa_path)
}
# 向控制台报告两个试验及总行数。
message(
  "输入QA通过：抗MRSA ",
  qa_table$n[qa_table$trial == "anti_mrsa"],
  " 行；抗PSA ",
  qa_table$n[qa_table$trial == "anti_psa"],
  " 行；合计 ",
  sum(qa_table$n),
  " 行。"
)
