# ==============================================================================
# Analysis: 独立R运行包统一配置
# Date: 2026-07-29
# Random seed: 20260727
# R: 4.4.2
#
# Method:
# Purpose: 集中声明路径、试验、随机种子、MICE、CBPS、结局和输出参数。
# Source type: local_knowledge_base / project_locked_protocol / self_written
# Source name or URL:
#   D:/source/obsidian/method_wiki/causal_inference/g_methods_ipw_standardization.md
#   analysis_plan_locked.md
# Citation or commit/version:
#   D-4.7N-ANALYSIS-LOCK-V2.0
#   D-4.6-ANALYSIS-OUTPUT-PLAN-V1.0
# Access date: 2026-07-29
# Adaptation notes:
#   本文件集中管理固定绝对根路径及全部派生路径，不实施统计估计，不改变锁定方法。
# Verification command:
#   parse(file="02_配置与变量字典/analysis_config.R")
# ==============================================================================

# 禁止R自动把字符列转换成因子，避免不同R版本产生隐式类型差异。
options(stringsAsFactors = FALSE)

# 从环境变量读取独立运行包根目录，避免公开代码包含本机绝对路径。
package_root_value <- Sys.getenv("TTE_RUNTIME_ROOT", unset = "")
if (!nzchar(package_root_value)) {
  stop("Set the TTE_RUNTIME_ROOT environment variable before running the analysis.")
}
package_root <- normalizePath(package_root_value, winslash = "/", mustWork = TRUE)

# 构造只读原始分析数据目录。
data_dir <- file.path(package_root, "01_原始分析数据")
# 构造配置、字典和锁定软件目录。
config_dir <- file.path(package_root, "02_配置与变量字典")
# 构造所有派生结果的唯一根目录。
output_root <- file.path(package_root, "09_输出结果")

# 声明两个平行试验及其固定内部代码。
trial_registry <- data.frame(
  # 使用稳定英文试验标识作为文件名和对象键。
  trial = c("anti_mrsa", "anti_psa"),
  # 使用1和2参与确定性随机种子派生。
  trial_code = c(1L, 2L),
  # 登记锁定风险集人数，供输入QA使用而不参与估计。
  expected_n = c(268L, 372L),
  # 登记降阶梯组人数，防止暴露编码或数据版本漂移。
  expected_deescalation_n = c(83L, 57L),
  # 登记继续组人数，防止暴露编码或数据版本漂移。
  expected_continuation_n = c(185L, 315L),
  # 保持字符列为字符类型。
  stringsAsFactors = FALSE
)

# 固定整个项目的主随机种子。
master_seed <- 20260727L
# 固定每次MICE生成20份插补数据。
mice_m <- 20L
# 固定每条MICE链运行50次迭代。
mice_maxit <- 50L
# 固定PMM候选供体数为5。
mice_donors <- 5L
# 固定每个平行试验执行1000次患者级Bootstrap。
bootstrap_B <- 1000L
# 固定主要AKI随访时点为T0后7天。
aki_horizon_days <- 7
# 固定数值交叉核验容差。
numeric_tolerance <- 1e-10

# 定义锁定CBPS-ATE补充估计对象使用的NEW_MS01倾向评分公式；当前投稿主分析为ATO。
ps_formula_main <- stats::as.formula(
  paste(
    "deescalation ~ age_years + male + target_class_grouped +",
    "diabetes_mellitus + congestive_heart_failure +",
    "chronic_obstructive_pulmonary_disease + chronic_kidney_disease +",
    "hematologic_malignancy + solid_malignancy + liver_disease_severity +",
    "immunosuppression + sepsis_onset_setting + infection_site_group5 +",
    "other_broad_spectrum_antibiotics + source_control_status +",
    "nephrotoxic_drug_exposure + baseline_creatinine_mg_dl +",
    "target_antibiotic_administration_count_72h + recent_nonrenal_sofa +",
    "nonrenal_sofa_change + creatinine_change_mg_dl +",
    "urine_output_rate_6h_ml_kg_h"
  )
)

# 定义将R对象原子化写入磁盘的函数，避免中断留下半个RDS文件。
atomic_save_rds <- function(object, path) {
  # 确保目标目录存在。
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  # 在同一目录构造临时文件名，以便后续原子重命名。
  temporary_path <- paste0(path, ".tmp")
  # 先把完整对象序列化到临时文件。
  saveRDS(object, temporary_path)
  # 若目标文件已存在则先删除，确保Windows下重命名行为确定。
  if (file.exists(path)) {
    # 仅删除调用者明确指定的单个派生文件。
    unlink(path)
  }
  # 将完整临时文件重命名为正式文件。
  renamed <- file.rename(temporary_path, path)
  # 重命名失败时停止，避免误报保存成功。
  if (!isTRUE(renamed)) {
    # 抛出包含目标路径的错误信息。
    stop("无法原子写入RDS：", path)
  }
  # 返回不可见TRUE供自动化测试判断。
  invisible(TRUE)
}

# 定义将数据框原子化写为CSV的函数。
atomic_write_csv <- function(data, path) {
  # 确保目标目录存在。
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  # 在同一目录构造临时CSV路径。
  temporary_path <- paste0(path, ".tmp")
  # 以UTF-8兼容方式写出数据，不保存行名。
  utils::write.csv(
    data,
    temporary_path,
    row.names = FALSE,
    na = ""
  )
  # 若目标文件已存在则仅删除该单个派生文件。
  if (file.exists(path)) {
    # 删除旧派生版本，为原子重命名腾出目标名。
    unlink(path)
  }
  # 将写完的临时文件重命名为正式CSV。
  renamed <- file.rename(temporary_path, path)
  # 重命名失败时立即停止。
  if (!isTRUE(renamed)) {
    # 抛出包含目标路径的错误信息。
    stop("无法原子写入CSV：", path)
  }
  # 返回不可见TRUE供调用脚本继续执行。
  invisible(TRUE)
}

# 定义依赖包检查函数，但不在分析脚本中自动联网安装软件。
require_packages <- function(packages) {
  # 对每个包执行静默可用性检查。
  available <- vapply(
    packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
  # 提取所有缺失包的名称。
  missing_packages <- packages[!available]
  # 存在缺失包时一次性报告完整清单。
  if (length(missing_packages) > 0L) {
    # 停止运行，避免使用未锁定的替代实现。
    stop("缺少R包：", paste(missing_packages, collapse = ", "))
  }
  # 返回不可见TRUE表示依赖检查通过。
  invisible(TRUE)
}

# 把项目隔离的mice与WeightIt/cobalt软件库放到R库搜索路径最前端。
.libPaths(c(
  file.path(
    config_dir,
    "R_library",
    "mice-github-61083e667fa41cc4bc655492591b100bf8a7e443"
  ),
  file.path(
    config_dir,
    "R_library",
    "phase47-weightit-cobalt-github"
  ),
  .libPaths()
))

# 固定全局随机种子；每个Bootstrap重复仍在引擎内使用确定性派生种子。
set.seed(master_seed)
