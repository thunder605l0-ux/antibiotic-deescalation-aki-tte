# Analysis: Build the locked English publication display dictionary
# Date: 2026-08-08
# Random seed: not applicable (deterministic dictionary build)
# R: 4.4.2
#
# Method: deterministic joins and exact semantic checks
# Purpose: separate machine-readable fields from reader-facing English labels
# Source type: project-locked dictionaries and controlled reporting tables
# Source name or URL: variable_coding_table_en/cn.csv; analysis_operation_variable_inventory_v1.0.csv;
#   main_iptw_operation_variable_dictionary_v1.0.csv; sf06_balance_term_labels_v1.0.csv;
#   ST-03_variable_definitions.csv
# Citation or commit/version: RD-4.6-ANALYSIS-OUTPUT-PLAN-V1.0; publication display dictionary v1.0
# Access date: 2026-08-08
# Adaptation notes: reader-facing race_white label is the researcher-approved "Race";
#   the machine field, binary coding, and missingness definition are unchanged.
# Verification command: D:/software/R4.4.2/R-4.4.2/bin/x64/Rscript.exe 02_配置与变量字典/21_build_publication_display_dictionary.R

options(stringsAsFactors = FALSE)

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(script_arg)) {
  normalizePath(sub("^--file=", "", script_arg[[1L]]), winslash = "/", mustWork = TRUE)
} else if (!is.null(sys.frame(1)$ofile)) {
  normalizePath(sys.frame(1)$ofile, winslash = "/", mustWork = TRUE)
} else {
  stop("Cannot determine script path.")
}
project_root <- dirname(dirname(script_path))
config_dir <- file.path(project_root, "02_配置与变量字典")
controlled_table_dir <- file.path(project_root, "09_输出结果", "09_tables", "controlled_reporting")

read_locked_csv <- function(path) {
  if (!file.exists(path)) stop("Missing locked source: ", path)
  utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE, fileEncoding = "UTF-8")
}

coding_en <- read_locked_csv(file.path(config_dir, "variable_coding_table_en.csv"))
coding_cn <- read_locked_csv(file.path(config_dir, "variable_coding_table_cn.csv"))
inventory <- read_locked_csv(file.path(config_dir, "analysis_operation_variable_inventory_v1.0.csv"))
dag_dictionary <- read_locked_csv(file.path(config_dir, "main_iptw_operation_variable_dictionary_v1.0.csv"))
term_dictionary <- read_locked_csv(file.path(config_dir, "sf06_balance_term_labels_v1.0.csv"))
st03 <- read_locked_csv(file.path(project_root, "09_输出结果", "09_tables", "ST-03_variable_definitions.csv"))

stopifnot(nrow(coding_en) == 52L, nrow(coding_cn) == 52L)
stopifnot(identical(coding_en$clean_field, coding_cn$clean_field))
stopifnot(!anyDuplicated(coding_en$clean_field))

field_labels <- c(
  patient_id = "Patient identifier",
  hospital_admission_id = "Hospital admission identifier",
  icu_stay_id = "ICU stay identifier",
  target_coverage_trial = "Target-coverage trial",
  treatment_strategy = "Treatment strategy",
  deescalation = "De-escalation strategy",
  sepsis_onset_time = "Sepsis onset time",
  target_antibiotic_start_time = "Target antibiotic initiation time",
  index_time = "Index time",
  sepsis_to_target_start_hours = "Time from sepsis onset to target antibiotic initiation",
  aki_followup_end_time = "End of AKI follow-up",
  aki_7d = "7-day incident AKI",
  aki_7d_onset_time = "Onset time of 7-day AKI",
  aki_7d_max_stage = "Maximum 7-day AKI stage",
  aki_7d_cross_t0_uo_sensitivity = "7-day AKI allowing urine-output windows across index time",
  aki_7d_cross_t0_uo_onset_time = "Onset time of 7-day AKI allowing urine-output windows across index time",
  aki_7d_cross_t0_uo_max_stage = "Maximum 7-day AKI stage allowing urine-output windows across index time",
  aki_7d_multistate_status = "7-day AKI multistate outcome",
  death_30d = "30-day all-cause mortality",
  death_30d_boundary_sensitivity = "30-day mortality boundary sensitivity outcome",
  death_time_from_index_days = "Time from index to death",
  aki_followup_days = "Observed AKI follow-up time",
  aki_7d_event_code = "AKI competing-event code",
  age_years = "Age",
  male = "Male sex",
  initial_target_antibiotic = "Initial target antibiotic class",
  diabetes_mellitus = "Diabetes mellitus",
  congestive_heart_failure = "Congestive heart failure",
  chronic_obstructive_pulmonary_disease = "Chronic obstructive pulmonary disease",
  chronic_kidney_disease = "Chronic kidney disease",
  hematologic_malignancy = "Hematologic malignancy",
  solid_malignancy = "Solid malignancy",
  liver_disease_severity = "Liver disease severity",
  immunosuppression = "Immunosuppression",
  sepsis_onset_setting = "Sepsis onset setting",
  infection_site_group5 = "Infection site",
  other_broad_spectrum_antibiotics = "Other broad-spectrum antibiotic exposure",
  source_control_status = "Source-control procedure",
  nephrotoxic_drug_exposure = "Nephrotoxic drug exposure",
  baseline_creatinine_mg_dl = "Baseline serum creatinine",
  target_antibiotic_administration_count_72h = "Target antibiotic administrations during the first 72 hours",
  early_nonrenal_sofa = "Early nonrenal SOFA score",
  recent_nonrenal_sofa = "Recent nonrenal SOFA score",
  nonrenal_sofa_change = "Change in nonrenal SOFA score",
  creatinine_change_mg_dl = "Pre-index change in serum creatinine",
  urine_output_rate_6h_ml_kg_h = "6-hour urine-output rate",
  urine_output_rate_12h_ml_kg_h = "12-hour urine-output rate",
  urine_output_rate_24h_ml_kg_h = "24-hour urine-output rate",
  race_white = "Race",
  skilled_nursing_facility_source = "Skilled nursing facility admission source",
  acute_viral_infection_status = "Acute viral infection status",
  microbiology_state_group4 = "Microbiology status at index time"
)

internal_labels_cn <- c(
  patient_id = "患者标识符", hospital_admission_id = "住院标识符", icu_stay_id = "ICU住院标识符",
  target_coverage_trial = "目标覆盖试验", treatment_strategy = "治疗策略", deescalation = "降阶梯策略",
  sepsis_onset_time = "脓毒症起始时间", target_antibiotic_start_time = "目标抗菌药物起始时间",
  index_time = "索引时间", sepsis_to_target_start_hours = "脓毒症至目标抗菌药物起始时间",
  aki_followup_end_time = "AKI随访结束时间", aki_7d = "7天新发AKI", aki_7d_onset_time = "7天AKI起始时间",
  aki_7d_max_stage = "7天AKI最高分期", aki_7d_cross_t0_uo_sensitivity = "允许尿量窗跨越索引时间的7天AKI",
  aki_7d_cross_t0_uo_onset_time = "允许尿量窗跨越索引时间的7天AKI起始时间",
  aki_7d_cross_t0_uo_max_stage = "允许尿量窗跨越索引时间的7天AKI最高分期",
  aki_7d_multistate_status = "7天AKI多状态结局", death_30d = "30天全因死亡",
  death_30d_boundary_sensitivity = "30天死亡边界敏感性结局", death_time_from_index_days = "索引至死亡时间",
  aki_followup_days = "AKI观察随访时间", aki_7d_event_code = "AKI竞争事件代码", age_years = "年龄",
  male = "男性", initial_target_antibiotic = "初始目标抗菌药物类别", diabetes_mellitus = "糖尿病",
  congestive_heart_failure = "充血性心力衰竭", chronic_obstructive_pulmonary_disease = "慢性阻塞性肺疾病",
  chronic_kidney_disease = "慢性肾病", hematologic_malignancy = "血液系统恶性肿瘤",
  solid_malignancy = "实体恶性肿瘤", liver_disease_severity = "肝病严重程度",
  immunosuppression = "免疫抑制", sepsis_onset_setting = "脓毒症起病场景", infection_site_group5 = "感染部位",
  other_broad_spectrum_antibiotics = "其他广谱抗菌药物暴露", source_control_status = "感染源控制操作",
  nephrotoxic_drug_exposure = "肾毒性药物暴露", baseline_creatinine_mg_dl = "基线血清肌酐",
  target_antibiotic_administration_count_72h = "最初72小时目标抗菌药物给药次数",
  early_nonrenal_sofa = "早期非肾脏SOFA评分", recent_nonrenal_sofa = "近期非肾脏SOFA评分",
  nonrenal_sofa_change = "非肾脏SOFA评分变化", creatinine_change_mg_dl = "索引前血清肌酐变化",
  urine_output_rate_6h_ml_kg_h = "6小时尿量率", urine_output_rate_12h_ml_kg_h = "12小时尿量率",
  urine_output_rate_24h_ml_kg_h = "24小时尿量率", race_white = "种族",
  skilled_nursing_facility_source = "专业护理机构入院来源", acute_viral_infection_status = "急性病毒感染状态",
  microbiology_state_group4 = "索引时间微生物学状态"
)

if (!setequal(names(field_labels), coding_en$clean_field)) stop("English field-label coverage is not exactly 52 fields.")
if (!setequal(names(internal_labels_cn), coding_en$clean_field)) stop("Chinese internal-label coverage is not exactly 52 fields.")

inventory_label <- setNames(inventory$variable, inventory$field)
dag_from_st03 <- setNames(st03$dag_construct, st03$clean_field)
dag_from_st03[dag_from_st03 %in% c("", "Not a covariate node")] <- NA_character_
dag_primary <- setNames(dag_dictionary$dag_concept, dag_dictionary$field)

publication_unit <- coding_en$unit
publication_unit[publication_unit == "None"] <- ""

dictionary <- data.frame(
  clean_field = coding_en$clean_field,
  display_order = coding_en$number,
  display_label = unname(field_labels[coding_en$clean_field]),
  table_label = unname(field_labels[coding_en$clean_field]),
  figure_label = unname(field_labels[coding_en$clean_field]),
  unit = publication_unit,
  variable_type = coding_en$variable_type,
  factor_levels_en = coding_en$numeric_code_definition,
  internal_label_cn = unname(internal_labels_cn[coding_en$clean_field]),
  factor_levels_cn = coding_cn$numeric_code_definition,
  dag_code = unname(dag_from_st03[coding_en$clean_field]),
  inventory_label_original = unname(inventory_label[coding_en$clean_field]),
  source_provenance = ifelse(
    coding_en$clean_field %in% inventory$field,
    "variable_coding_table_en.csv; analysis_operation_variable_inventory_v1.0.csv; ST-03_variable_definitions.csv",
    "variable_coding_table_en.csv; variable_coding_table_cn.csv; ST-03_variable_definitions.csv"
  ),
  active = TRUE,
  stringsAsFactors = FALSE
)

dictionary$dag_code[is.na(dictionary$dag_code)] <- unname(dag_primary[dictionary$clean_field[is.na(dictionary$dag_code)]])
dictionary$dag_code[is.na(dictionary$dag_code)] <- ""

# Preserve exact coding semantics and explicitly document the approved presentation-only race label.
race_row <- dictionary[dictionary$clean_field == "race_white", , drop = FALSE]
stopifnot(nrow(race_row) == 1L)
stopifnot(race_row$display_label == "Race")
stopifnot(race_row$factor_levels_en == "0=explicitly non-White; 1=White; unknown race=missing")

qa <- data.frame(
  check_id = c(
    "DICT01", "DICT02", "DICT03", "DICT04", "DICT05", "DICT06", "DICT07", "DICT08",
    "DICT09", "DICT10", "DICT11", "DICT12"
  ),
  description = c(
    "English coding table has exactly 52 rows",
    "Chinese coding table has exactly 52 rows",
    "English and Chinese clean fields match in order",
    "All clean fields are unique",
    "Every clean field has exactly one nonempty English display label",
    "Every clean field has exactly one internal Chinese label",
    "No reader-facing field label contains an underscore",
    "No reader-facing field label exposes factor syntax",
    "race_white is displayed as Race",
    "race_white coding and missingness definition are unchanged",
    "Publication units exactly preserve the English coding table units",
    "All 33 locked balance terms retain a nonempty English display label"
  ),
  pass = c(
    nrow(coding_en) == 52L,
    nrow(coding_cn) == 52L,
    identical(coding_en$clean_field, coding_cn$clean_field),
    !anyDuplicated(dictionary$clean_field),
    nrow(dictionary) == 52L && all(nzchar(dictionary$display_label)),
    all(nzchar(dictionary$internal_label_cn)),
    !any(grepl("_", dictionary$display_label, fixed = TRUE)),
    !any(grepl("factor(", dictionary$display_label, fixed = TRUE)),
    identical(race_row$display_label, "Race"),
    identical(race_row$factor_levels_en, "0=explicitly non-White; 1=White; unknown race=missing"),
    identical(dictionary$unit, publication_unit),
    nrow(term_dictionary) == 33L && !anyDuplicated(term_dictionary$term) && all(nzchar(term_dictionary$display_label))
  ),
  stringsAsFactors = FALSE
)

term_dictionary$display_label[term_dictionary$term == "race_white"] <- "Race"
direct_term_index <- match(term_dictionary$term, dictionary$clean_field)
direct_term_rows <- !is.na(direct_term_index)
term_dictionary$display_label[direct_term_rows] <- dictionary$display_label[direct_term_index[direct_term_rows]]
direct_units <- dictionary$unit[direct_term_index[direct_term_rows]]
has_unit <- nzchar(direct_units)
term_dictionary$display_label[which(direct_term_rows)[has_unit]] <- paste0(
  term_dictionary$display_label[which(direct_term_rows)[has_unit]], ", ", direct_units[has_unit]
)
term_dictionary$display_label <- gsub("Initial target class:", "Initial target antibiotic class:", term_dictionary$display_label, fixed = TRUE)
term_dictionary$source_provenance <- "sf06_balance_term_labels_v1.0.csv; publication display dictionary v1.0"
term_dictionary$active <- TRUE

if (!all(qa$pass)) {
  print(qa[!qa$pass, , drop = FALSE])
  stop("Publication display dictionary QA failed.")
}

utils::write.csv(
  dictionary,
  file.path(config_dir, "publication_display_labels_v1.0.csv"),
  row.names = FALSE,
  na = "",
  fileEncoding = "UTF-8"
)
utils::write.csv(
  term_dictionary,
  file.path(config_dir, "publication_model_term_labels_v1.0.csv"),
  row.names = FALSE,
  na = "",
  fileEncoding = "UTF-8"
)
utils::write.csv(
  qa,
  file.path(config_dir, "publication_display_labels_v1.0_QA.csv"),
  row.names = FALSE,
  na = "",
  fileEncoding = "UTF-8"
)

cat("PUBLICATION_DISPLAY_DICTIONARY_ROWS=", nrow(dictionary), "\n", sep = "")
cat("PUBLICATION_MODEL_TERM_ROWS=", nrow(term_dictionary), "\n", sep = "")
cat("PUBLICATION_DISPLAY_DICTIONARY_QA=PASS\n")
