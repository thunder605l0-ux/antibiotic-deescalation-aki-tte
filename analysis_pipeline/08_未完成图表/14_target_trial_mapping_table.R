# ==============================================================================
# ST-01：目标试验方案与MIMIC-IV映射
# Method: 将锁定method.md中的目标试验要素结构化为补充表
# Purpose: 生成可审计的目标试验方案映射CSV和DOCX
# Source type: article_supplement + self_written
# Source name or URL: 文献a补充材料的target trial specification table结构
# Citation or commit/version: 本项目锁定决定D-1.5B-DAY4-LANDMARK-IPTW与D-3.11-DESIGN-LOCK-V1.0
# Access date: 2026-07-29
# Adaptation notes: 内容严格按本项目锁定方法，不机械复制文献a的人群资格差异。
# Verification command: Rscript 08_未完成图表/14_target_trial_mapping_table.R
# ==============================================================================

# 读取统一配置。
source(
  file.path(Sys.getenv("TTE_RUNTIME_ROOT", unset = "."), "02_配置与变量字典", "analysis_config.R"),
  encoding = "UTF-8"
)
setwd(package_root)

# 检查表格输出依赖。
require_packages(c("flextable", "officer"))

# 逐项定义目标试验与观察性模拟。
target_trial_table <- data.frame(
  Element = c(
    "Target population",
    "Parallel trials",
    "Treatment strategies",
    "Treatment assignment",
    "Time zero",
    "Grace period",
    "Follow-up",
    "Primary outcome",
    "Competing events",
    "Secondary outcome",
    "Causal contrast",
    "Analysis"
  ),
  Hypothetical_target_trial = c(
    "Adults with ICU sepsis who remain eligible at the 72-hour decision point after initiation of trial-specific broad-spectrum coverage.",
    "Anti-MRSA coverage trial and anti-pseudomonal/resistant Gram-negative coverage trial.",
    "Immediate discontinuation of trial-specific target coverage versus continuation at time zero.",
    "Random assignment at time zero.",
    "72 hours after the first eligible administration of the trial-specific target antibiotic.",
    "None.",
    "From time zero through 7 days for AKI and through 30 days for all-cause mortality.",
    "New-onset in-hospital KDIGO AKI within 7 days.",
    "Death before AKI and live discharge before AKI.",
    "30-day all-cause mortality.",
    "Average treatment effect in the overlap population (ATO): effect of immediate discontinuation versus continuation among eligible patients for whom both strategies were clinically plausible at the 72-hour decision point.",
    "Marginal risks, risk difference as the primary effect measure, and risk ratio as a supplementary measure."
  ),
  MIMIC_IV_emulation = c(
    "First identifiable hospitalization and first ICU stay; sepsis used to identify the index population; trial-specific target coverage documented during hours 48–72; alive, observable, AKI-free, and without target-MDRO positivity at time zero.",
    "Constructed and analyzed separately; one hospitalization may contribute to both trials when separately eligible.",
    "Strategy determined only from orders known by time zero and applicable after time zero; no post-time-zero order is used to classify treatment.",
    "Observed clinical assignment adjusted using the locked prespecified adjustment-set covariates.",
    "First eligible eMAR administration plus exactly 72 hours.",
    "None; no CCW, cloning, or artificial censoring.",
    "Starts synchronously with eligibility and assignment at time zero.",
    "Three-path phenotype using creatinine, urine output, and new exact-time RRT.",
    "Handled with weighted Aalen–Johansen estimation.",
    "Exact in-hospital death time preferred; date-only death used under the locked boundary rule.",
    "Primary analysis: overlap weighting from bias-reduced logistic propensity-score models (de-escalated: 1−PS[X]; continued: PS[X]). Complementary analysis: stabilised CBPS ATE-IPTW for the full eligible population.",
    "Twenty imputations; primary ATO propensity scores and overlap weights; nested patient bootstrap with imputation, weighting, and outcome estimation repeated within each resample. Complementary CBPS-ATE results were reported separately."
  ),
  stringsAsFactors = FALSE
)

# 定义输出目录。
table_dir <- file.path(output_root, "09_tables")

# 创建输出目录。
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

# 保存统计数据层CSV。
atomic_write_csv(
  target_trial_table,
  file.path(table_dir, "ST-01_target_trial_mapping.csv")
)

# 创建三线表对象。
target_trial_flex <- flextable::flextable(target_trial_table)

# 设置表头。
target_trial_flex <- flextable::set_header_labels(
  target_trial_flex,
  Element = "Target trial element",
  Hypothetical_target_trial = "Hypothetical target trial",
  MIMIC_IV_emulation = "MIMIC-IV emulation"
)

# 应用简洁学术表格格式。
target_trial_flex <- flextable::theme_booktabs(target_trial_flex)

# 自动调整列宽。
target_trial_flex <- flextable::autofit(target_trial_flex)

# 保存DOCX。
flextable::save_as_docx(
  "Supplementary Table ST-01" = target_trial_flex,
  path = file.path(table_dir, "ST-01_target_trial_mapping.docx")
)

# 输出完成信息。
message("ST-01已生成。")
