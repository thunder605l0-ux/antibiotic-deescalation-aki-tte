# ==============================================================================
# ST-02：抗菌药清单及MIMIC-IV映射
# Method: 合并锁定目标药词表与系统性给药途径词表
# Purpose: 生成抗MRSA/抗PSA药物、数据库映射和途径规则补充表
# Source type: article_supplement + self_written
# Source name or URL: 文献a抗菌药定义；本项目锁定药物和途径词表
# Citation or commit/version: D-3.5-TIME-EXPOSURE-TARGET-ABX与D-3.5-DAY3-SYSTEMIC-DICTIONARY-ROUTE
# Access date: 2026-07-29
# Adaptation notes: 不从分析结果反向增删药物；只格式化既有锁定词表。
# Verification command: Rscript 08_未完成图表/15_antibacterial_dictionary_table.R
# ==============================================================================

# 读取统一配置。
source(
  file.path(Sys.getenv("TTE_RUNTIME_ROOT", unset = "."), "02_配置与变量字典", "analysis_config.R"),
  encoding = "UTF-8"
)
setwd(package_root)

# 检查表格依赖。
require_packages(c("flextable", "officer"))

# 读取目标抗菌药词表。
target_dictionary <- utils::read.csv(
  file.path(config_dir, "target_antibiotic_dictionary.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# 读取系统性抗菌药原始映射词表。
systemic_dictionary <- utils::read.csv(
  file.path(config_dir, "systemic_antibacterial_mimic_raw_dictionary.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# 读取给药途径词表。
route_dictionary <- utils::read.csv(
  file.path(config_dir, "systemic_antibiotic_route_dictionary.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# 为三个词表增加来源标签后纵向保存，避免猜测不一致的连接键。
target_dictionary$table_section <- "Trial-specific target antibiotics"
systemic_dictionary$table_section <- "MIMIC medication-name mapping"
route_dictionary$table_section <- "Route inclusion and exclusion rules"

# 统一三个数据框的全部字段。
all_columns <- unique(c(
  names(target_dictionary),
  names(systemic_dictionary),
  names(route_dictionary)
))

# 补齐目标药词表缺失字段。
for (column_name in setdiff(all_columns, names(target_dictionary))) {
  target_dictionary[[column_name]] <- NA_character_
}

# 补齐系统性药物词表缺失字段。
for (column_name in setdiff(all_columns, names(systemic_dictionary))) {
  systemic_dictionary[[column_name]] <- NA_character_
}

# 补齐途径词表缺失字段。
for (column_name in setdiff(all_columns, names(route_dictionary))) {
  route_dictionary[[column_name]] <- NA_character_
}

# 按统一列顺序合并。
dictionary_table <- rbind(
  target_dictionary[, all_columns, drop = FALSE],
  systemic_dictionary[, all_columns, drop = FALSE],
  route_dictionary[, all_columns, drop = FALSE]
)

# 定义输出目录。
table_dir <- file.path(output_root, "09_tables")

# 创建输出目录。
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

# 保存完整可审计CSV。
atomic_write_csv(
  dictionary_table,
  file.path(table_dir, "ST-02_antibacterial_dictionary.csv")
)

# 创建Word文档；395行药名映射使用officer原生表以避免大型flextable渲染过慢。
dictionary_doc <- officer::read_docx()

# 写入目标药分表标题。
dictionary_doc <- officer::body_add_par(
  dictionary_doc,
  "ST-02A Trial-specific target antibiotics",
  style = "heading 1"
)

# 写入目标药分表。
dictionary_doc <- officer::body_add_table(
  dictionary_doc,
  value = target_dictionary
)

# 写入数据库药名映射分表标题。
dictionary_doc <- officer::body_add_par(
  dictionary_doc,
  "ST-02B MIMIC medication-name mapping",
  style = "heading 1"
)

# 写入数据库药名映射分表。
dictionary_doc <- officer::body_add_table(
  dictionary_doc,
  value = systemic_dictionary
)

# 写入给药途径分表标题。
dictionary_doc <- officer::body_add_par(
  dictionary_doc,
  "ST-02C Route inclusion and exclusion rules",
  style = "heading 1"
)

# 写入给药途径分表。
dictionary_doc <- officer::body_add_table(
  dictionary_doc,
  value = route_dictionary
)

# 保存Word版本。
print(
  dictionary_doc,
  target = file.path(table_dir, "ST-02_antibacterial_dictionary.docx")
)

# 输出完成信息。
message("ST-02已生成。")
