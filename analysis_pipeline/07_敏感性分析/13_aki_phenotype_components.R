# ==============================================================================
# ST-10：AKI电子表型组成
# Method: 将肌酐、尿量和新启动RRT路径按试验与策略交叉汇总
# Purpose: 透明报告三路径表型及其重叠，不改变主要复合AKI定义
# Source type: local_knowledge_base + self_written
# Source name or URL: D:/source/obsidian；本项目锁定三路径KDIGO结局对象
# Citation or commit/version: D-3.5-THREE-PATH-KDIGO-OUTCOME-V1.1
# Access date: 2026-07-29
# Adaptation notes: 旧中间表仅按最终640键内连接后使用，未进入最终风险集的记录全部剔除。
# Verification command: Rscript 07_敏感性分析/13_aki_phenotype_components.R
# ==============================================================================

# 读取统一配置。
source(
  file.path(Sys.getenv("TTE_RUNTIME_ROOT", unset = "."), "02_配置与变量字典", "analysis_config.R"),
  encoding = "UTF-8"
)
setwd(package_root)

# 检查表格依赖。
require_packages(c("flextable", "officer"))

# 读取三路径患者级中间表。
components <- utils::read.csv(
  file.path(data_dir, "aki_three_path_patient_summary.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# 读取两个最终试验数据并建立最终键。
final_keys <- do.call(
  rbind,
  lapply(trial_registry$trial, function(trial_name) {
    current <- utils::read.csv(
      file.path(data_dir, paste0("analysis_dataset_", trial_name, ".csv")),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    data.frame(
      stay_id = current$icu_stay_id,
      trial = trial_name,
      strategy_final = current$treatment_strategy,
      aki_7d_final = current$aki_7d,
      stringsAsFactors = FALSE
    )
  })
)

# 用stay_id和试验双键限制至最终640条记录。
components_final <- merge(
  final_keys,
  components,
  by = c("stay_id", "trial"),
  all.x = TRUE,
  sort = FALSE
)

# 要求每个最终键恰好连接一次。
stopifnot(nrow(components_final) == sum(trial_registry$expected_n))

# 要求关键路径字段均成功连接。
required_component_fields <- c(
  "creatinine_aki",
  "uo_aki_fully_post_t0",
  "rrt_event",
  "aki7_primary"
)

# 关键路径缺失时停止。
if (anyNA(components_final[, required_component_fields])) {
  stop("最终640键与AKI三路径中间表连接不完整。")
}

# 要求主要AKI与最终分析CSV一致。
stopifnot(
  all(
    as.integer(components_final$aki7_primary) ==
      as.integer(components_final$aki_7d_final)
  )
)

# 生成互斥路径组合标签。
components_final$path_combination <- with(
  components_final,
  paste0(
    ifelse(creatinine_aki, "C", ""),
    ifelse(uo_aki_fully_post_t0, "U", ""),
    ifelse(rrt_event, "R", "")
  )
)

# 将无路径事件标为None。
components_final$path_combination[
  components_final$path_combination == ""
] <- "None"

# 定义逐组汇总函数。
summarize_group <- function(data) {
  # 返回路径事件数与总AKI数。
  data.frame(
    n = nrow(data),
    aki_any_n = sum(data$aki7_primary),
    creatinine_path_n = sum(data$creatinine_aki),
    urine_output_path_n = sum(data$uo_aki_fully_post_t0),
    rrt_path_n = sum(data$rrt_event),
    creatinine_only_n = sum(data$path_combination == "C"),
    urine_only_n = sum(data$path_combination == "U"),
    rrt_only_n = sum(data$path_combination == "R"),
    creatinine_urine_n = sum(data$path_combination == "CU"),
    creatinine_rrt_n = sum(data$path_combination == "CR"),
    urine_rrt_n = sum(data$path_combination == "UR"),
    all_three_n = sum(data$path_combination == "CUR"),
    stringsAsFactors = FALSE
  )
}

# 按试验和策略拆分。
split_groups <- split(
  components_final,
  list(components_final$trial, components_final$strategy_final),
  drop = TRUE
)

# 汇总全部组。
summary_rows <- lapply(names(split_groups), function(group_name) {
  current <- split_groups[[group_name]]
  cbind(
    trial = unique(current$trial),
    treatment_strategy = unique(current$strategy_final),
    summarize_group(current),
    stringsAsFactors = FALSE
  )
})

# 合并汇总结果。
phenotype_table <- do.call(rbind, summary_rows)

# 定义输出目录。
sensitivity_dir <- file.path(output_root, "06_敏感性分析")

# 创建输出目录。
dir.create(sensitivity_dir, recursive = TRUE, showWarnings = FALSE)

# 保存CSV。
atomic_write_csv(
  phenotype_table,
  file.path(sensitivity_dir, "ST-10_aki_phenotype_components.csv")
)

# 创建Word三线表。
phenotype_flex <- flextable::flextable(phenotype_table)

# 应用三线表格式。
phenotype_flex <- flextable::theme_booktabs(phenotype_flex)

# 自动调整列宽。
phenotype_flex <- flextable::autofit(phenotype_flex)

# 保存DOCX。
flextable::save_as_docx(
  "Supplementary Table ST-10" = phenotype_flex,
  path = file.path(
    sensitivity_dir,
    "ST-10_aki_phenotype_components.docx"
  )
)

# 保存用于独立核对的患者级路径组合。
atomic_write_csv(
  components_final[, c(
    "stay_id",
    "trial",
    "strategy_final",
    "aki7_primary",
    "creatinine_aki",
    "uo_aki_fully_post_t0",
    "rrt_event",
    "path_combination"
  )],
  file.path(sensitivity_dir, "ST-10_aki_phenotype_components_patient_QA.csv")
)

# 输出完成信息。
message("ST-10已生成并限制于最终640键。")
