# ==============================================================================
# 脚本用途：对同时进入 anti-MRSA 与 antipseudomonal 两个目标试验的患者，
#   进行两项 T0 决策的 2×2 事后描述性交叉分类，并报告各格 7 天 AKI 风险。
# 方法来源：项目锁定的两个试验级分析数据集；Clopper-Pearson 精确二项区间
#   使用 R stats::binom.test（R 4.4.2）。
# 引用 / commit：RD-4.6-ANALYSIS-OUTPUT-PLAN-V1.0；两个锁定试验数据集的
#   SHA-256 由项目 manifest 追踪。
# 访问日期：2026-08-28。
# 适配说明：本分析为事后描述性分析。两个试验的 T0 不强制相同，因此分别
#   报告 anti-MRSA 锚定与 antipseudomonal 锚定的 AKI 风险；不拟合联合倾向
#   评分、不使用既有单试验 ATO 权重、不运行 bootstrap、不计算 P 值或正式交互。
# 验证命令：Rscript --vanilla 23_joint_strategy_crossclassification_2026-08-28.R
# ==============================================================================

options(stringsAsFactors = FALSE, warn = 1)

project_root <- normalizePath(Sys.getenv("TTE_PROJECT_ROOT", unset = "."), winslash = "/", mustWork = TRUE)
runtime_root <- normalizePath(Sys.getenv("TTE_RUNTIME_ROOT", unset = "."), winslash = "/", mustWork = TRUE)
raw_dir <- file.path(runtime_root, "01_原始分析数据")
output_dir <- file.path(project_root, "data")

mrsa_path <- file.path(raw_dir, "analysis_dataset_anti_mrsa.csv")
psa_path <- file.path(raw_dir, "analysis_dataset_anti_psa.csv")
table_path <- file.path(output_dir, "Supplementary Table 7_Joint Strategy Cross-Classification.csv")
audit_path <- file.path(output_dir, "Supplementary Table 7_Joint Strategy Cross-Classification_audit.csv")

mrsa <- utils::read.csv(mrsa_path, check.names = FALSE, stringsAsFactors = FALSE)
psa <- utils::read.csv(psa_path, check.names = FALSE, stringsAsFactors = FALSE)

key_fields <- c("patient_id", "hospital_admission_id", "icu_stay_id")
required_fields <- c(key_fields, "index_time", "deescalation", "aki_7d")
if (!all(required_fields %in% names(mrsa)) || !all(required_fields %in% names(psa))) {
  stop("The locked trial datasets do not contain the required joint-analysis fields.")
}
if (anyDuplicated(mrsa[key_fields]) || anyDuplicated(psa[key_fields])) {
  stop("Trial datasets contain duplicated participant keys.")
}

shared <- merge(
  mrsa[, required_fields],
  psa[, required_fields],
  by = key_fields,
  suffixes = c("_mrsa", "_psa"),
  all = FALSE,
  sort = FALSE
)
if (nrow(shared) != 171L) stop("Expected 171 shared participants; found ", nrow(shared), ".")

parse_time <- function(x) {
  as.POSIXct(x, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
}
shared$t0_mrsa <- parse_time(shared$index_time_mrsa)
shared$t0_psa <- parse_time(shared$index_time_psa)
if (anyNA(shared$t0_mrsa) || anyNA(shared$t0_psa)) stop("A shared-participant T0 could not be parsed.")
shared$t0_difference_hours <- as.numeric(difftime(shared$t0_psa, shared$t0_mrsa, units = "hours"))
shared$absolute_t0_difference_hours <- abs(shared$t0_difference_hours)
shared$t0_within_24h <- shared$absolute_t0_difference_hours <= 24

shared$deescalation_mrsa <- as.integer(shared$deescalation_mrsa)
shared$deescalation_psa <- as.integer(shared$deescalation_psa)
shared$aki_7d_mrsa <- as.integer(shared$aki_7d_mrsa)
shared$aki_7d_psa <- as.integer(shared$aki_7d_psa)
binary_fields <- c("deescalation_mrsa", "deescalation_psa", "aki_7d_mrsa", "aki_7d_psa")
if (any(!vapply(shared[binary_fields], function(x) all(x %in% 0:1), logical(1)))) {
  stop("A strategy or AKI field is not binary in the shared cohort.")
}

shared$mrsa_strategy <- ifelse(shared$deescalation_mrsa == 1, "De-escalation", "Continuation")
shared$psa_strategy <- ifelse(shared$deescalation_psa == 1, "De-escalation", "Continuation")
shared$cell_order <- 1L + 2L * shared$deescalation_mrsa + shared$deescalation_psa

exact_interval <- function(events, total) {
  interval <- stats::binom.test(events, total, conf.level = 0.95)$conf.int
  c(risk = events / total, lower = interval[1], upper = interval[2])
}

format_count_percent <- function(count, total) {
  sprintf("%d (%.1f)", count, 100 * count / total)
}

format_events_n <- function(events, total) sprintf("%d/%d", events, total)

format_risk_ci <- function(estimate) {
  sprintf("%.1f (%.1f to %.1f)", 100 * estimate["risk"], 100 * estimate["lower"], 100 * estimate["upper"])
}

audit_rows <- list()
display_rows <- list()

for (cell in 1:4) {
  selected <- shared$cell_order == cell
  cell_data <- shared[selected, , drop = FALSE]
  total <- nrow(cell_data)
  if (total == 0L) stop("A joint-strategy cell is empty.")
  events_mrsa <- sum(cell_data$aki_7d_mrsa)
  events_psa <- sum(cell_data$aki_7d_psa)
  estimate_mrsa <- exact_interval(events_mrsa, total)
  estimate_psa <- exact_interval(events_psa, total)
  within_24h <- sum(cell_data$t0_within_24h)

  audit_rows[[cell]] <- data.frame(
    order = cell,
    anti_mrsa_strategy = unique(cell_data$mrsa_strategy),
    antipseudomonal_strategy = unique(cell_data$psa_strategy),
    n = total,
    t0_within_24h_n = within_24h,
    t0_within_24h_proportion = within_24h / total,
    anti_mrsa_anchored_aki_events = events_mrsa,
    anti_mrsa_anchored_risk = unname(estimate_mrsa["risk"]),
    anti_mrsa_anchored_ci_lower = unname(estimate_mrsa["lower"]),
    anti_mrsa_anchored_ci_upper = unname(estimate_mrsa["upper"]),
    antipseudomonal_anchored_aki_events = events_psa,
    antipseudomonal_anchored_risk = unname(estimate_psa["risk"]),
    antipseudomonal_anchored_ci_lower = unname(estimate_psa["lower"]),
    antipseudomonal_anchored_ci_upper = unname(estimate_psa["upper"]),
    stringsAsFactors = FALSE
  )

  display_rows[[cell]] <- data.frame(
    `Anti-MRSA strategy` = unique(cell_data$mrsa_strategy),
    `Antipseudomonal strategy` = unique(cell_data$psa_strategy),
    N = total,
    `T0 within 24 h, n (%)` = format_count_percent(within_24h, total),
    `Anti-MRSA-anchored AKI events/N` = format_events_n(events_mrsa, total),
    `Anti-MRSA-anchored 7-day AKI risk (95% CI), %` = format_risk_ci(estimate_mrsa),
    `Antipseudomonal-anchored AKI events/N` = format_events_n(events_psa, total),
    `Antipseudomonal-anchored 7-day AKI risk (95% CI), %` = format_risk_ci(estimate_psa),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

audit <- do.call(rbind, audit_rows)
display <- do.call(rbind, display_rows)

overall <- data.frame(
  shared_n = nrow(shared),
  exact_same_t0_n = sum(shared$absolute_t0_difference_hours == 0),
  t0_within_24h_n = sum(shared$t0_within_24h),
  t0_within_24h_proportion = mean(shared$t0_within_24h),
  absolute_t0_difference_median_hours = stats::median(shared$absolute_t0_difference_hours),
  absolute_t0_difference_q1_hours = unname(stats::quantile(shared$absolute_t0_difference_hours, 0.25, type = 7)),
  absolute_t0_difference_q3_hours = unname(stats::quantile(shared$absolute_t0_difference_hours, 0.75, type = 7)),
  individual_aki_concordant_n = sum(shared$aki_7d_mrsa == shared$aki_7d_psa),
  individual_aki_concordant_proportion = mean(shared$aki_7d_mrsa == shared$aki_7d_psa),
  stringsAsFactors = FALSE
)

audit$shared_n <- overall$shared_n
audit$exact_same_t0_n <- overall$exact_same_t0_n
audit$overall_t0_within_24h_n <- overall$t0_within_24h_n
audit$absolute_t0_difference_median_hours <- overall$absolute_t0_difference_median_hours
audit$absolute_t0_difference_q1_hours <- overall$absolute_t0_difference_q1_hours
audit$absolute_t0_difference_q3_hours <- overall$absolute_t0_difference_q3_hours
audit$individual_aki_concordant_n <- overall$individual_aki_concordant_n

utils::write.csv(display, table_path, row.names = FALSE, fileEncoding = "UTF-8")
utils::write.csv(audit, audit_path, row.names = FALSE, fileEncoding = "UTF-8")

cat("SHARED_N=", overall$shared_n, "\n", sep = "")
cat("EXACT_SAME_T0_N=", overall$exact_same_t0_n, "\n", sep = "")
cat("T0_WITHIN_24H_N=", overall$t0_within_24h_n, "\n", sep = "")
cat(
  "ABSOLUTE_T0_DIFFERENCE_HOURS_MEDIAN_IQR=",
  sprintf(
    "%.1f (%.1f to %.1f)",
    overall$absolute_t0_difference_median_hours,
    overall$absolute_t0_difference_q1_hours,
    overall$absolute_t0_difference_q3_hours
  ),
  "\n",
  sep = ""
)
cat("INDIVIDUAL_AKI_CONCORDANT_N=", overall$individual_aki_concordant_n, "\n", sep = "")
print(display, row.names = FALSE)
