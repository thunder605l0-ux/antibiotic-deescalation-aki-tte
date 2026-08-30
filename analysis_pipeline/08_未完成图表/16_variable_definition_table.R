# ==============================================================================
# ST-03：变量编码、时间窗、DAG角色和缺失处理
# Method: 合并变量编码表与MICE工作变量字典
# Purpose: 形成统计分析使用的完整变量定义补充表
# Source type: self_written
# Source name or URL: 本项目锁定变量编码表、DAG映射和MICE规格
# Citation or commit/version: D-3.4C-PRIMARY-ADJUSTMENT-SET-V1.0与D-4.7N-ANALYSIS-LOCK-V2.0
# Access date: 2026-07-29
# Adaptation notes: 仅连接和格式化已有定义；连接失败字段保留为空并在QA表中报告。
# Verification command: Rscript 08_未完成图表/16_variable_definition_table.R
# ==============================================================================

# 读取统一配置。
source(
  file.path(Sys.getenv("TTE_RUNTIME_ROOT", unset = "."), "02_配置与变量字典", "analysis_config.R"),
  encoding = "UTF-8"
)
setwd(package_root)

# 检查表格依赖。
require_packages(c("flextable", "officer"))

# 读取英文变量编码表。
coding_en <- utils::read.csv(
  file.path(config_dir, "variable_coding_table_en.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# 读取中文变量编码表，作为审计对照。
coding_cn <- utils::read.csv(
  file.path(config_dir, "variable_coding_table_cn.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# 读取MICE工作变量字典。
mice_dictionary <- utils::read.csv(
  file.path(config_dir, "mice_working_variable_dictionary.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# 明确把编码表clean_field映射到MICE字典field，避免按位置连接。
mice_dictionary$clean_field <- mice_dictionary$field

# 锁定共同连接键。
join_key <- "clean_field"

# 合并编码表和MICE字典。
variable_table <- merge(
  coding_en,
  mice_dictionary,
  by = join_key,
  all.x = TRUE,
  sort = FALSE,
  suffixes = c("_coding", "_mice")
)

# 恢复英文编码表原始顺序。
variable_table$.order <- match(variable_table[[join_key]], coding_en[[join_key]])
variable_table <- variable_table[order(variable_table$.order), , drop = FALSE]
variable_table$.order <- NULL

# 为全部52列建立可审计的时间窗、DAG角色、模型状态和数据来源。
variable_table$measurement_window <- "Not applicable or fixed at cohort construction"
variable_table$dag_construct <- "Not a covariate node"
variable_table$main_model_status <- "Not included in propensity-score model"
variable_table$missing_data_handling <- "Not imputed"
variable_table$data_source <- "Locked MIMIC-IV extraction and derivation chain"

# 定义按字段批量更新元数据的辅助函数。
set_metadata <- function(fields, column, value) {
  index <- variable_table$clean_field %in% fields
  variable_table[index, column] <<- value
}

# 患者、住院和ICU键以及人口学来源。
set_metadata(
  c("patient_id", "hospital_admission_id", "icu_stay_id"),
  "data_source",
  "MIMIC-IV patients/admissions/icustays identifiers"
)
set_metadata(
  c("age_years", "male", "race_white"),
  "data_source",
  "MIMIC-IV patients and admissions"
)
set_metadata(
  c("age_years", "male", "race_white"),
  "measurement_window",
  "Known by T0; demographic baseline"
)

# 暴露、时间锚点和抗菌治疗字段。
antibiotic_fields <- c(
  "target_coverage_trial",
  "treatment_strategy",
  "deescalation",
  "target_antibiotic_start_time",
  "sepsis_to_target_start_hours",
  "initial_target_antibiotic",
  "other_broad_spectrum_antibiotics",
  "target_antibiotic_administration_count_72h"
)
set_metadata(
  antibiotic_fields,
  "data_source",
  "MIMIC-IV eMAR, eMAR detail, pharmacy and POE-derived antibiotic records"
)
set_metadata(
  c("sepsis_onset_time", "sepsis_onset_setting"),
  "data_source",
  "MIMIC-IV Sepsis-3 derived cohort and admission/ICU times"
)
set_metadata(
  c("target_antibiotic_start_time", "initial_target_antibiotic"),
  "measurement_window",
  "First eligible target-antibiotic administration"
)
set_metadata(
  "index_time",
  "measurement_window",
  "First eligible target-antibiotic administration plus 72 hours"
)
set_metadata(
  "sepsis_to_target_start_hours",
  "measurement_window",
  "Sepsis time to first eligible target-antibiotic administration"
)
set_metadata(
  "target_antibiotic_administration_count_72h",
  "measurement_window",
  "From first target administration through T0"
)
set_metadata(
  c("treatment_strategy", "deescalation"),
  "measurement_window",
  "Strategy known and applicable at T0; no post-T0 grace period"
)

# 结局和竞争事件来源及随访窗口。
aki_fields <- c(
  "aki_followup_end_time",
  "aki_7d",
  "aki_7d_onset_time",
  "aki_7d_max_stage",
  "aki_7d_cross_t0_uo_sensitivity",
  "aki_7d_cross_t0_uo_onset_time",
  "aki_7d_cross_t0_uo_max_stage",
  "aki_7d_multistate_status",
  "aki_followup_days",
  "aki_7d_event_code"
)
set_metadata(
  aki_fields,
  "data_source",
  "MIMIC-IV labevents, chartevents, outputevents, procedureevents and admissions"
)
set_metadata(
  aki_fields,
  "measurement_window",
  "From T0 through T0+7 days, death, or live discharge, whichever occurs first"
)
death_fields <- c(
  "death_30d",
  "death_30d_boundary_sensitivity",
  "death_time_from_index_days"
)
set_metadata(
  death_fields,
  "data_source",
  "MIMIC-IV admissions deathtime and patients date of death"
)
set_metadata(
  death_fields,
  "measurement_window",
  "From T0 through T0+30 days"
)

# 共病和基线感染特征。
comorbidity_fields <- c(
  "diabetes_mellitus",
  "congestive_heart_failure",
  "chronic_obstructive_pulmonary_disease",
  "chronic_kidney_disease",
  "hematologic_malignancy",
  "solid_malignancy",
  "liver_disease_severity",
  "immunosuppression"
)
set_metadata(
  comorbidity_fields,
  "data_source",
  "MIMIC-IV diagnoses_icd and prespecified clinical code lists"
)
set_metadata(
  comorbidity_fields,
  "measurement_window",
  "Pre-existing condition known by T0"
)
set_metadata(
  c("infection_site_group5", "microbiology_state_group4"),
  "data_source",
  "T0-available structured infection information and microbiologyevents"
)
set_metadata(
  c("infection_site_group5", "microbiology_state_group4"),
  "measurement_window",
  "Specimens collected by T0; result status restricted to information available by T0"
)
set_metadata(
  "source_control_status",
  "data_source",
  "MIMIC-IV procedureevents and prespecified source-control procedure mapping"
)
set_metadata(
  "source_control_status",
  "measurement_window",
  "Source-control procedure completed before T0"
)
set_metadata(
  "nephrotoxic_drug_exposure",
  "data_source",
  "MIMIC-IV eMAR/pharmacy nephrotoxic-drug exposure mapping"
)
set_metadata(
  "nephrotoxic_drug_exposure",
  "measurement_window",
  "Exposure documented before T0"
)

# 肾功能、尿量和非肾SOFA轨迹。
set_metadata(
  c("baseline_creatinine_mg_dl", "creatinine_change_mg_dl"),
  "data_source",
  "MIMIC-IV labevents creatinine measurements"
)
set_metadata(
  "baseline_creatinine_mg_dl",
  "measurement_window",
  "Prespecified baseline value available before T0"
)
set_metadata(
  "creatinine_change_mg_dl",
  "measurement_window",
  "T0-preceding creatinine trajectory; no post-T0 measurement"
)
set_metadata(
  c(
    "early_nonrenal_sofa",
    "recent_nonrenal_sofa",
    "nonrenal_sofa_change"
  ),
  "data_source",
  "MIMIC-IV chart/laboratory/treatment components of non-renal SOFA"
)
set_metadata(
  "early_nonrenal_sofa",
  "measurement_window",
  "Early pretreatment/T0-preceding severity window"
)
set_metadata(
  "recent_nonrenal_sofa",
  "measurement_window",
  "Most recent fully observed non-renal SOFA window ending no later than T0"
)
set_metadata(
  "nonrenal_sofa_change",
  "measurement_window",
  "Recent minus early non-renal SOFA, both windows ending no later than T0"
)
urine_fields <- c(
  "urine_output_rate_6h_ml_kg_h",
  "urine_output_rate_12h_ml_kg_h",
  "urine_output_rate_24h_ml_kg_h"
)
set_metadata(
  urine_fields,
  "data_source",
  "MIMIC-IV outputevents/chartevents and prespecified body-weight denominator"
)
set_metadata(
  "urine_output_rate_6h_ml_kg_h",
  "measurement_window",
  "Complete 6-hour support window ending no later than T0"
)
set_metadata(
  "urine_output_rate_12h_ml_kg_h",
  "measurement_window",
  "Complete 12-hour support window ending no later than T0"
)
set_metadata(
  "urine_output_rate_24h_ml_kg_h",
  "measurement_window",
  "Complete 24-hour support window ending no later than T0"
)

# 锁定NEW_MS01操作变量与DAG构念一一映射。
dag_map <- c(
  age_years = "AGE",
  male = "SEX",
  initial_target_antibiotic = "TARGET_CLASS",
  diabetes_mellitus = "DM",
  congestive_heart_failure = "CHF",
  chronic_obstructive_pulmonary_disease = "COPD",
  chronic_kidney_disease = "CKD",
  hematologic_malignancy = "HEM_MAL",
  solid_malignancy = "SOLID_MAL",
  liver_disease_severity = "LIVER_DIS",
  immunosuppression = "IMMUNOSUPP",
  sepsis_onset_setting = "ONSET_SETTING",
  infection_site_group5 = "INF_SITE",
  other_broad_spectrum_antibiotics = "OTHER_BSA",
  source_control_status = "SOURCE_CTRL",
  nephrotoxic_drug_exposure = "NEPHRO_DRUG",
  baseline_creatinine_mg_dl = "SCR_BASE",
  target_antibiotic_administration_count_72h = "TARGET_EXPOSURE",
  recent_nonrenal_sofa = "CLIN_TRAJ",
  nonrenal_sofa_change = "CLIN_TRAJ",
  creatinine_change_mg_dl = "SCR_DELTA",
  urine_output_rate_6h_ml_kg_h = "UO_RATE"
)
for (field in names(dag_map)) {
  set_metadata(field, "dag_construct", unname(dag_map[field]))
  set_metadata(field, "main_model_status", "Included in NEW_MS01")
}
set_metadata(
  "sepsis_to_target_start_hours",
  "dag_construct",
  "ABX_TIME"
)
set_metadata(
  "sepsis_to_target_start_hours",
  "main_model_status",
  "Sensitivity only: NEW_MS01 plus ABX_TIME"
)
set_metadata("race_white", "dag_construct", "RACE_W")
set_metadata(
  "race_white",
  "main_model_status",
  "Sensitivity only: NEW_MS03 replaces SEX"
)
set_metadata(
  "microbiology_state_group4",
  "dag_construct",
  "MICRO_STATE"
)
set_metadata(
  "microbiology_state_group4",
  "main_model_status",
  "Proxy-disturbance sensitivity only"
)
set_metadata(
  c(
    "early_nonrenal_sofa",
    "urine_output_rate_12h_ml_kg_h",
    "urine_output_rate_24h_ml_kg_h"
  ),
  "main_model_status",
  "MICE auxiliary variable; not entered jointly in the PS formula"
)

# 从锁定MICE字典填入实际插补方法；非MICE字段保留Not imputed。
for (row_index in seq_len(nrow(variable_table))) {
  main_method <- variable_table$planned_method_main[row_index]
  if (!is.na(main_method) && nzchar(main_method)) {
    variable_table$missing_data_handling[row_index] <- paste0(
      "MICE: ",
      main_method
    )
  }
}
set_metadata(
  "nonrenal_sofa_change",
  "missing_data_handling",
  "Passive derivation: recent non-renal SOFA minus early non-renal SOFA"
)
set_metadata(
  "race_white",
  "missing_data_handling",
  "Not imputed in NEW_MS01; logistic MICE in NEW_MS03"
)
set_metadata(
  c("infection_site_group5", "microbiology_state_group4"),
  "missing_data_handling",
  "Prespecified information-state category; not treated as ordinary missingness"
)

# 定义输出目录。
table_dir <- file.path(output_root, "09_tables")

# 创建输出目录。
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

# 保存主补充表CSV。
atomic_write_csv(
  variable_table,
  file.path(table_dir, "ST-03_variable_definitions.csv")
)

# 保存中英文键一致性QA。
key_qa <- data.frame(
  item = c(
    "English coding rows",
    "Chinese coding rows",
    "English unique keys",
    "Chinese unique keys",
    "Keys present only in English",
    "Keys present only in Chinese"
  ),
  value = c(
    nrow(coding_en),
    nrow(coding_cn),
    length(unique(coding_en[[join_key]])),
    length(unique(coding_cn[[join_key]])),
    length(setdiff(coding_en[[join_key]], coding_cn[[join_key]])),
    length(setdiff(coding_cn[[join_key]], coding_en[[join_key]]))
  ),
  stringsAsFactors = FALSE
)

# 输出键一致性QA。
atomic_write_csv(
  key_qa,
  file.path(table_dir, "ST-03_variable_definition_key_QA.csv")
)

# 创建Word三线表。
variable_flex <- flextable::flextable(variable_table)

# 应用三线表风格。
variable_flex <- flextable::theme_booktabs(variable_flex)

# 自动调整版式。
variable_flex <- flextable::autofit(variable_flex)

# 保存DOCX。
flextable::save_as_docx(
  "Supplementary Table ST-03" = variable_flex,
  path = file.path(table_dir, "ST-03_variable_definitions.docx")
)

# 输出完成信息。
message("ST-03已生成。")
