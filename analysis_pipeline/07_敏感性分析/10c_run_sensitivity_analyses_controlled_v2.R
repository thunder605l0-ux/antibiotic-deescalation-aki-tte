# ==============================================================================
# ST-08：预设敏感性分析受控恢复 V2（trial worker）
# Method: 复用锁定敏感性函数；仅对NEW_MS03的race_white绑定官方logreg等价
#         条件重试，并复用mice_remove_lindep_verified_v2精确事件审计。
# Source type: project_locked_code + upstream_package_source + method_paper
# Source name or URL:
#   1) 07_敏感性分析/10_run_sensitivity_analyses.R
#   2) 07_敏感性分析/10b_run_sensitivity_analyses_controlled.R
#   3) https://github.com/amices/mice/blob/
#      61083e667fa41cc4bc655492591b100bf8a7e443/R/mice.impute.logreg.R
# Citation or commit/version:
#   1) locked sensitivity script MD5 9c20d440972610c41e193cc7bf31eab1
#   2) mice commit 61083e667fa41cc4bc655492591b100bf8a7e443
#   3) White IR, Daniel R, Royston P. Stat Med. 2010;29:2920-2931.
# Access date: 2026-08-02
# Adaptation notes:
#   - race_white首次内部glm.fit保留官方默认maxit=25。
#   - 仅对精确“glm.fit: algorithm did not converge”且
#     fit$converged=FALSE，以相同增广数据、quasibinomial-logit和起始规则
#     将maxit提高至2000重试一次。
#   - 其他警告、重试警告、重试失败或未知loggedEvents均硬停止。
#   - 不改变m=20、MICE链maxit=50、donors=5、预测矩阵、访问顺序、
#     bootstrap编号、checkpoint保存种子、权重或结局估计量。
# Verification command:
#   D:/software/R4.4.2/R-4.4.2/bin/x64/Rscript.exe --vanilla
#   07_敏感性分析/10c_run_sensitivity_analyses_controlled_v2.R
#   --trial=anti_mrsa --b_start=31 --b_end=31
#
# CBPS 条件式重试扩展（2026-08-03）：
#   - 方法依据：Imai K, Ratkovic M. J R Stat Soc Series B.
#     2014;76(1):243-263. doi:10.1111/rssb.12027.
#   - 上游代码：WeightIt 1.7.0.9004, commit
#     595d2812c2f4d56ff7c14aaf6fcde76c59aa8c8f。
#   - 仅对精确 CBPS 非收敛警告且 convergence code != 0，
#     在同数据、公式、方法、起始规则和收敛阈值下将
#     maxit 从默认 5000 提高至 20000 重试一次。
#   - 适配边界、审计字段和验证命令见：
#     docs/superpowers/specs/
#     2026-08-03-st08-cbps-single-retry-controlled-v2-design.md。
#   - 访问日期：2026-08-03。
#
# nephrotoxic_drug_exposure logreg100 最终层救援（2026-08-04）：
#   - 仅在初始maxit=100和第一层maxit=2000均出现唯一、精确的
#     “glm.fit: algorithm did not converge”且fit$converged=FALSE时，
#     以同一增广数据、quasibinomial-logit和零起始规则进行一次且仅一次
#     maxit=20000最终拟合。
#   - 该规则仅存在于ST-08 NEW_MS03，不改写全局MICE扩展或主分析。
#   - 诊断、门槛、审计和验证顺序见：docs/superpowers/specs/
#     2026-08-04-st08-b0079-logreg100-diagnostic-and-rescue-design.md。
#   - 方法来源：mice.impute.logreg，commit
#     61083e667fa41cc4bc655492591b100bf8a7e443；R stats::glm.control；
#     White IR, Daniel R, Royston P. Stat Med. 2010;29:2920-2931。
#   - 访问日期：2026-08-04。
#
# 零风险比边界验证扩展（2026-08-04）：
#   - 仅当risk_1恰为0、risk_0严格大于0、RR与risk_1/risk_0一致，
#     且全部风险、RD、RR和既有QA通过时，允许有限边界RR=0。
#   - 负RR、零分母、非有限值、风险越界、恒等式不一致或未知警告
#     仍然硬停止；不加入连续性校正或伪计数，不修改估计量。
#   - 方法与实施边界见：docs/superpowers/specs/
#     2026-08-04-st08-b0373-zero-risk-ratio-boundary-design.md。
#   - 代码依据：项目锁定RR公式（10_run_sensitivity_analyses.R，
#     MD5=9c20d440972610c41e193cc7bf31eab1）；R 4.4.2 base
#     Arithmetic Operators；访问日期：2026-08-04。
#
# race_white有序对比列线性依赖审计扩展（2026-08-05）：
#   - 仅允许dep=race_white、meth=racewhitelogreg，且out严格为
#     other_broad_spectrum_antibiotics.L/Q或其逗号标准组合。
#   - 运行时同时核验race_white锁定方法为logreg、二分类水平为0/1，
#     有序抗菌药因子水平为0/1/2，且预测矩阵明确包含该预测变量。
#   - 方法来源：mice 3.19.0内部remove.lindep()与updateLog()；
#     amices/mice commit 61083e667fa41cc4bc655492591b100bf8a7e443；
#     访问日期：2026-08-05。
#   - 适配边界、保护证据和验证命令见：docs/superpowers/specs/
#     2026-08-05-st08-b0584-race-white-ordered-contrast-lindep-design.md。
#
# b=635事件路由纠正（2026-08-05）：
#   - 仅将dep=race_white、meth=racewhitelogreg且out精确为有序抗菌药
#     L/Q列的事件交给race_white线性依赖规则。
#   - aki_competing_statedeath_before_aki事件继续且仅由既有零计数规则
#     逐事件审核；不新增白名单、不改变插补或估计方法。
#   - 诊断、路由字段和验证顺序见：docs/superpowers/specs/
#     2026-08-05-st08-b0635-event-routing-correction-design.md。
#   - 代码来源：R 4.4.2逻辑运算；mice 3.19.0 remove.lindep()/
#     updateLog()；访问日期：2026-08-05。
# ==============================================================================

# 定位当前脚本；同时支持Rscript和source()。
get_current_script_path_v2 <- function() {
  command_line <- commandArgs(trailingOnly = FALSE)
  file_argument <- grep("^--file=", command_line, value = TRUE)
  if (length(file_argument) == 1L) {
    return(normalizePath(
      sub("^--file=", "", file_argument),
      winslash = "/",
      mustWork = TRUE
    ))
  }
  frame_files <- vapply(
    sys.frames(),
    function(frame) {
      if (is.null(frame$ofile)) "" else as.character(frame$ofile)
    },
    character(1)
  )
  frame_files <- frame_files[nzchar(frame_files)]
  if (length(frame_files) > 0L) {
    return(normalizePath(
      frame_files[length(frame_files)],
      winslash = "/",
      mustWork = TRUE
    ))
  }
  stop("无法定位10c脚本；请使用Rscript或source()运行。")
}

script_path_v2 <- get_current_script_path_v2()
package_root_v2 <- normalizePath(
  file.path(dirname(script_path_v2), ".."),
  winslash = "/",
  mustWork = TRUE
)
controlled_v2_worker_path <- file.path(
  package_root_v2,
  "07_敏感性分析",
  "10c_run_sensitivity_analyses_controlled_v2.R"
)
if (!file.exists(controlled_v2_worker_path)) {
  stop("缺少10c受控V2 worker脚本：", controlled_v2_worker_path)
}
controlled_v1_path <- file.path(
  package_root_v2,
  "07_敏感性分析",
  "10b_run_sensitivity_analyses_controlled.R"
)
if (!file.exists(controlled_v1_path)) {
  stop("缺少10b受控脚本：", controlled_v1_path)
}

# 从10b读取初始化及函数定义，不执行其worker或finalizer。
controlled_v1_lines <- readLines(
  controlled_v1_path,
  warn = FALSE,
  encoding = "UTF-8"
)
controlled_v1_driver_marker <- which(
  trimws(controlled_v1_lines) ==
    "# 逐试验、逐正式编号计算；已验证RDS直接续用。"
)
if (length(controlled_v1_driver_marker) != 1L) {
  stop("10b执行分界标记不唯一；停止以避免错误加载。")
}
eval(
  parse(text = controlled_v1_lines[
    seq_len(controlled_v1_driver_marker - 1L)
  ]),
  envir = .GlobalEnv
)

# 固定并核验本恢复方案批准时的关键文件。
expected_md5_v2 <- c(
  locked_sensitivity = "9c20d440972610c41e193cc7bf31eab1",
  controlled_v1 = "91245d6f22b53e99b968c1f33b161e42",
  mice_extension = "a1443a8a91a547d6b8b53ae30d71bc45",
  anti_mrsa_b0031_failure = "c42b8cd40e4290b0a2f69c6734220792",
  anti_psa_b0934_failure = "2e9b244d16b9020fd94fb4aa4824d619",
  anti_psa_b0052_v2_failure = "6a79fce5ecfd1636601960d37bf5752f",
  anti_mrsa_b0140_v2_failure = "c14bf320ab4eb617bc4eed647991f4b5",
  anti_psa_b0079_v2_failure = "b5b0bedda197620decb0e4ec2c1bad68",
  anti_psa_b0079_checkpoint = "4e0f1283896c7c52006c9f8a7e6f94ec",
  anti_psa_b0079_diagnostic = "23205ba97392aa992495230882f72c49",
  anti_mrsa_b0373_v2_failure = "b6d67a8d123692fb9fb56ef2a8284b24",
  anti_mrsa_b0373_checkpoint = "7c2f68af03006e4b6414aa69c7a59717",
  anti_mrsa_b0373_diagnostic = "605aaed07047465238967ae0a5178561",
  anti_mrsa_b0584_v2_failure = "fb53fe63c9dcd698f04236301478e1de",
  anti_mrsa_b0584_checkpoint = "e97ab0bd9d517ce9441960f9abf83233",
  anti_mrsa_b0584_diagnostic = "c1f386a47d9d67fac5effe6e009c0bff",
  anti_mrsa_b0584_preflight = "94d30cbde4d2400f9a48a0ce842b8628",
  anti_mrsa_b0584_worker_before = "e0f4966cada1d67155861fccc91dd1d9",
  anti_mrsa_b0635_v2_failure = "29e131366953ee45952d9ab72980065b",
  anti_mrsa_b0635_checkpoint = "508599b143e2f61a06fb54d3d20477f0",
  anti_mrsa_b0635_diagnostic = "6e10395c04ebb97b73fa31094baa5a0a",
  anti_mrsa_b0635_diagnostic_script = "25edf613e2e4452cefb034875307d4d5",
  anti_mrsa_b0635_worker_before = "fc9707b7f3ede605d87b7674f986a981"
)
mice_extension_path_v2 <- file.path(
  config_dir,
  "mice_imputation_extensions.R"
)
reviewed_b0031_failure_path_v2 <- file.path(
  failure_dir,
  "anti_mrsa_b0031_sensitivity_failure.rds"
)

assert_md5_v2 <- function(path, expected, label) {
  if (!file.exists(path)) {
    stop("缺少", label, "：", path)
  }
  observed <- tolower(unname(tools::md5sum(path)))
  if (!identical(observed, tolower(expected))) {
    stop(
      label,
      " MD5不匹配；expected=",
      expected,
      "；observed=",
      observed
    )
  }
  invisible(observed)
}

assert_md5_v2(
  locked_script_path,
  expected_md5_v2[["locked_sensitivity"]],
  "锁定敏感性脚本"
)
assert_md5_v2(
  controlled_v1_path,
  expected_md5_v2[["controlled_v1"]],
  "10b受控脚本"
)
assert_md5_v2(
  mice_extension_path_v2,
  expected_md5_v2[["mice_extension"]],
  "MICE扩展"
)
assert_md5_v2(
  reviewed_b0031_failure_path_v2,
  expected_md5_v2[["anti_mrsa_b0031_failure"]],
  "anti-MRSA b=31敏感性失败证据"
)
assert_md5_v2(
  preserved_failure_path,
  expected_md5_v2[["anti_psa_b0934_failure"]],
  "anti-PSA b=934失败证据"
)

repeat_audit_dir_v2 <- file.path(sensitivity_dir, "repeat_audits_v2")
failure_dir_v2 <- file.path(sensitivity_dir, "failures_v2")
dir.create(repeat_audit_dir_v2, recursive = TRUE, showWarnings = FALSE)
dir.create(failure_dir_v2, recursive = TRUE, showWarnings = FALSE)

reviewed_cbps_failure_registry_v2 <- data.frame(
  trial = c("anti_psa", "anti_mrsa"),
  bootstrap_index = c(52L, 140L),
  file_name = c(
    "anti_psa_b0052_sensitivity_controlled_v2_failure.rds",
    "anti_mrsa_b0140_sensitivity_controlled_v2_failure.rds"
  ),
  failure_md5 = c(
    expected_md5_v2[["anti_psa_b0052_v2_failure"]],
    expected_md5_v2[["anti_mrsa_b0140_v2_failure"]]
  ),
  checkpoint_md5 = c(
    "fe45638c507447b7723f747432a3c55d",
    "5c18b0b0c7598528a02fdac1a285939b"
  ),
  stringsAsFactors = FALSE
)
reviewed_cbps_failure_registry_v2$path <- file.path(
  failure_dir_v2,
  reviewed_cbps_failure_registry_v2$file_name
)
reviewed_cbps_failure_message_v2 <- paste0(
  "cbps_ate拟合或权重QA失败；warnings=",
  "The optimization failed to converge; try again with a higher value of `maxit`."
)

is_reviewed_cbps_failure_v2 <- function(path, trial, bootstrap_index) {
  registry_row <- reviewed_cbps_failure_registry_v2[
    reviewed_cbps_failure_registry_v2$trial == trial &
      reviewed_cbps_failure_registry_v2$bootstrap_index == bootstrap_index,
    ,
    drop = FALSE
  ]
  if (nrow(registry_row) != 1L || !file.exists(path)) return(FALSE)
  if (!identical(basename(path), registry_row$file_name[[1]])) return(FALSE)
  if (!identical(
    tolower(unname(tools::md5sum(path))),
    tolower(registry_row$failure_md5[[1]])
  )) return(FALSE)

  evidence <- readRDS(path)
  identical(evidence$status, "failed") &&
    identical(evidence$execution_version, "controlled_v2") &&
    identical(evidence$trial, trial) &&
    identical(as.integer(evidence$bootstrap_index), bootstrap_index) &&
    identical(
      tolower(evidence$checkpoint_md5),
      tolower(registry_row$checkpoint_md5[[1]])
    ) &&
    identical(evidence$message, reviewed_cbps_failure_message_v2) &&
    identical(
      tolower(evidence$locked_script_md5),
      expected_md5_v2[["locked_sensitivity"]]
    )
}

for (registry_index in seq_len(nrow(reviewed_cbps_failure_registry_v2))) {
  registry_row <- reviewed_cbps_failure_registry_v2[
    registry_index,
    ,
    drop = FALSE
  ]
  if (!is_reviewed_cbps_failure_v2(
    registry_row$path[[1]],
    registry_row$trial[[1]],
    registry_row$bootstrap_index[[1]]
  )) {
    stop(
      "已审核 CBPS V2 failure 证据不匹配：",
      registry_row$path[[1]]
    )
  }
}

# 登记 b=79 原始 logreg100 失败和只读诊断证据；两者必须永久保留。
reviewed_b0079_failure_path_v2 <- file.path(
  failure_dir_v2,
  "anti_psa_b0079_sensitivity_controlled_v2_failure.rds"
)
reviewed_b0079_diagnostic_path_v2 <- file.path(
  sensitivity_dir,
  "diagnostic_history",
  "st08_logreg100_anti_psa_b0079",
  "diagnostic_evidence_20260804_084312.rds"
)
reviewed_b0079_failure_message_v2 <- paste0(
  "logreg100内部glm.fit条件式重试失败；警告=",
  "glm.fit: algorithm did not converge；converged=FALSE"
)

is_reviewed_b0079_failure_v2 <- function(path) {
  if (!file.exists(path)) return(FALSE)
  if (!identical(
    tolower(unname(tools::md5sum(path))),
    expected_md5_v2[["anti_psa_b0079_v2_failure"]]
  )) return(FALSE)
  evidence <- readRDS(path)
  identical(evidence$status, "failed") &&
    identical(evidence$execution_version, "controlled_v2") &&
    identical(evidence$trial, "anti_psa") &&
    identical(as.integer(evidence$bootstrap_index), 79L) &&
    identical(
      tolower(evidence$checkpoint_md5),
      expected_md5_v2[["anti_psa_b0079_checkpoint"]]
    ) &&
    identical(evidence$message, reviewed_b0079_failure_message_v2) &&
    identical(
      tolower(evidence$locked_script_md5),
      expected_md5_v2[["locked_sensitivity"]]
    ) &&
    identical(
      tolower(evidence$controlled_v2_script_md5),
      "219d51bd947baef03641ea89d23b0f6a"
    )
}

is_verified_b0079_diagnostic_v2 <- function(path) {
  if (!file.exists(path)) return(FALSE)
  if (!identical(
    tolower(unname(tools::md5sum(path))),
    expected_md5_v2[["anti_psa_b0079_diagnostic"]]
  )) return(FALSE)
  evidence <- readRDS(path)
  metadata <- evidence$metadata
  candidate <- evidence$candidate_fits
  comparison <- evidence$numerical_comparison
  required_maxit <- c(100L, 2000L, 5000L, 10000L, 20000L)
  is.list(evidence) &&
    identical(evidence$status, "diagnostic_complete") &&
    is.data.frame(metadata) &&
    nrow(metadata) == 1L &&
    identical(metadata$trial, "anti_psa") &&
    identical(as.integer(metadata$bootstrap_index), 79L) &&
    identical(metadata$specification, "NEW_MS03") &&
    identical(metadata$target_field, "nephrotoxic_drug_exposure") &&
    isTRUE(metadata$eligible_for_controlled_rescue) &&
    identical(
      tolower(metadata$failure_md5_before),
      expected_md5_v2[["anti_psa_b0079_v2_failure"]]
    ) &&
    identical(
      tolower(metadata$failure_md5_after),
      expected_md5_v2[["anti_psa_b0079_v2_failure"]]
    ) &&
    is.data.frame(candidate) &&
    all(required_maxit %in% candidate$internal_maxit) &&
    all(candidate$converged[candidate$internal_maxit >= 5000L]) &&
    all(candidate$warning_n[candidate$internal_maxit >= 5000L] == 0L) &&
    is.data.frame(comparison) &&
    nrow(comparison) == 1L &&
    isTRUE(comparison$reference_converged) &&
    isTRUE(comparison$target_converged) &&
    comparison$max_absolute_coefficient_difference <= 1e-6 &&
    comparison$max_absolute_probability_difference <= 1e-8 &&
    identical(evidence$formal_files_written, FALSE)
}

if (!is_reviewed_b0079_failure_v2(reviewed_b0079_failure_path_v2)) {
  stop("anti-PSA b=79原始V2失败证据不匹配。")
}
if (!is_verified_b0079_diagnostic_v2(reviewed_b0079_diagnostic_path_v2)) {
  stop("anti-PSA b=79只读诊断证据未通过预设门槛。")
}

# 登记 b=373 原始零RR失败和只读诊断证据；两者必须永久保留。
reviewed_b0373_failure_path_v2 <- file.path(
  failure_dir_v2,
  "anti_mrsa_b0373_sensitivity_controlled_v2_failure.rds"
)
reviewed_b0373_diagnostic_path_v2 <- file.path(
  sensitivity_dir,
  "diagnostic_history",
  "st08_zero_rr_anti_mrsa_b0373",
  "diagnostic_evidence_20260804_225520.rds"
)
reviewed_b0373_failure_message_v2 <-
  "单重复敏感性结果包含非正RR。"

is_reviewed_b0373_failure_v2 <- function(path) {
  if (!file.exists(path)) return(FALSE)
  if (!identical(
    tolower(unname(tools::md5sum(path))),
    expected_md5_v2[["anti_mrsa_b0373_v2_failure"]]
  )) return(FALSE)
  evidence <- readRDS(path)
  identical(evidence$status, "failed") &&
    identical(evidence$execution_version, "controlled_v2") &&
    identical(evidence$trial, "anti_mrsa") &&
    identical(as.integer(evidence$bootstrap_index), 373L) &&
    identical(
      tolower(evidence$checkpoint_md5),
      expected_md5_v2[["anti_mrsa_b0373_checkpoint"]]
    ) &&
    identical(evidence$message, reviewed_b0373_failure_message_v2) &&
    identical(
      tolower(evidence$locked_script_md5),
      expected_md5_v2[["locked_sensitivity"]]
    ) &&
    identical(
      tolower(evidence$controlled_v2_script_md5),
      "fba5c3861dca833f63a2111c949774b1"
    )
}

is_verified_b0373_diagnostic_v2 <- function(path) {
  if (!file.exists(path)) return(FALSE)
  if (!identical(
    tolower(unname(tools::md5sum(path))),
    expected_md5_v2[["anti_mrsa_b0373_diagnostic"]]
  )) return(FALSE)
  evidence <- readRDS(path)
  metadata <- evidence$metadata
  boundary <- evidence$zero_rr_boundary_summary
  event_count <- evidence$complete_case_event_count
  is.list(evidence) &&
    identical(evidence$status, "diagnostic_complete") &&
    is.data.frame(metadata) &&
    nrow(metadata) == 1L &&
    identical(metadata$trial, "anti_mrsa") &&
    identical(as.integer(metadata$bootstrap_index), 373L) &&
    identical(metadata$specification, "complete_case") &&
    identical(metadata$outcome, "death_30d") &&
    identical(as.integer(metadata$warning_n), 0L) &&
    identical(as.integer(metadata$complete_case_n), 230L) &&
    identical(as.integer(metadata$continuation_death_n), 19L) &&
    identical(as.integer(metadata$deescalation_death_n), 0L) &&
    identical(as.integer(metadata$zero_rr_boundary_n), 1L) &&
    isTRUE(metadata$all_key_values_finite) &&
    isTRUE(metadata$eligible_for_zero_rr_boundary) &&
    identical(
      tolower(metadata$failure_md5_before),
      expected_md5_v2[["anti_mrsa_b0373_v2_failure"]]
    ) &&
    identical(
      tolower(metadata$failure_md5_after),
      expected_md5_v2[["anti_mrsa_b0373_v2_failure"]]
    ) &&
    identical(
      tolower(metadata$checkpoint_md5_before),
      expected_md5_v2[["anti_mrsa_b0373_checkpoint"]]
    ) &&
    identical(
      tolower(metadata$checkpoint_md5_after),
      expected_md5_v2[["anti_mrsa_b0373_checkpoint"]]
    ) &&
    identical(
      tolower(metadata$controlled_v2_script_md5),
      "7bb47fb770e468aeca987fc51fc5dff4"
    ) &&
    identical(
      tolower(metadata$diagnostic_script_md5),
      "c9e1450c3472dbcf6a30d8ea64df9d4c"
    ) &&
    identical(metadata$formal_files_written, FALSE) &&
    identical(evidence$formal_files_written, FALSE) &&
    is.data.frame(boundary) &&
    nrow(boundary) == 1L &&
    identical(boundary$trial, "anti_mrsa") &&
    identical(as.integer(boundary$bootstrap_index), 373L) &&
    identical(boundary$specification, "complete_case") &&
    identical(boundary$outcome, "death_30d") &&
    boundary$risk_0 > 0 &&
    identical(boundary$risk_1, 0) &&
    identical(boundary$rr, 0) &&
    isTRUE(boundary$denominator_positive) &&
    isTRUE(boundary$risk_in_range) &&
    boundary$rd_identity_error <= 1e-12 &&
    boundary$rr_identity_error <= 1e-12 &&
    isTRUE(boundary$boundary_verified) &&
    is.data.frame(event_count) &&
    sum(event_count$n) == 230L &&
    sum(event_count$n[
      event_count$deescalation == 1L & event_count$death_30d == 1L
    ]) == 0L
}

if (!is_reviewed_b0373_failure_v2(reviewed_b0373_failure_path_v2)) {
  stop("anti-MRSA b=373原始V2失败证据不匹配。")
}
if (!is_verified_b0373_diagnostic_v2(reviewed_b0373_diagnostic_path_v2)) {
  stop("anti-MRSA b=373零RR只读诊断证据未通过预设门槛。")
}

# 登记 b=584 原始MICE事件失败、只读复现和候选预检；全部永久保留。
reviewed_b0584_failure_path_v2 <- file.path(
  failure_dir_v2,
  "anti_mrsa_b0584_sensitivity_controlled_v2_failure.rds"
)
reviewed_b0584_diagnostic_dir_v2 <- file.path(
  sensitivity_dir,
  "diagnostic_history",
  "st08_new_ms03_anti_mrsa_b0584"
)
reviewed_b0584_diagnostic_path_v2 <- file.path(
  reviewed_b0584_diagnostic_dir_v2,
  "diagnostic_evidence_20260805_065811.rds"
)
reviewed_b0584_preflight_path_v2 <- file.path(
  reviewed_b0584_diagnostic_dir_v2,
  "candidate_exact_whitelist_preflight_20260805_070402.rds"
)
reviewed_b0584_failure_message_v2 <- paste0(
  "NEW_MS03 MICE未产生完整MIDS或出现未批准loggedEvents；事件数=46"
)
race_white_lindep_strategy_v3 <-
  "race_white_ordered_contrast_lindep_verified_v1"

is_exact_b0584_event_summary_v3 <- function(summary) {
  if (!is.data.frame(summary) || nrow(summary) != 2L) return(FALSE)
  summary <- summary[order(summary$out), , drop = FALSE]
  identical(summary$dep, rep("race_white", 2L)) &&
    identical(summary$meth, rep("racewhitelogreg", 2L)) &&
    identical(
      summary$out,
      c(
        "other_broad_spectrum_antibiotics.L",
        "other_broad_spectrum_antibiotics.Q"
      )
    ) &&
    identical(as.integer(summary$event_n), c(20L, 26L))
}

is_reviewed_b0584_failure_v2 <- function(path) {
  if (!file.exists(path)) return(FALSE)
  if (!identical(
    tolower(unname(tools::md5sum(path))),
    expected_md5_v2[["anti_mrsa_b0584_v2_failure"]]
  )) return(FALSE)
  evidence <- readRDS(path)
  identical(evidence$status, "failed") &&
    identical(evidence$execution_version, "controlled_v2") &&
    identical(evidence$trial, "anti_mrsa") &&
    identical(as.integer(evidence$bootstrap_index), 584L) &&
    identical(
      tolower(evidence$checkpoint_md5),
      expected_md5_v2[["anti_mrsa_b0584_checkpoint"]]
    ) &&
    identical(evidence$message, reviewed_b0584_failure_message_v2) &&
    identical(
      tolower(evidence$locked_script_md5),
      expected_md5_v2[["locked_sensitivity"]]
    ) &&
    identical(
      tolower(evidence$controlled_v2_script_md5),
      expected_md5_v2[["anti_mrsa_b0584_worker_before"]]
    )
}

is_verified_b0584_diagnostic_v2 <- function(path) {
  if (!file.exists(path)) return(FALSE)
  if (!identical(
    tolower(unname(tools::md5sum(path))),
    expected_md5_v2[["anti_mrsa_b0584_diagnostic"]]
  )) return(FALSE)
  evidence <- readRDS(path)
  is.list(evidence) &&
    identical(evidence$status, "diagnostic_completed") &&
    identical(evidence$diagnostic_only, TRUE) &&
    identical(evidence$trial, "anti_mrsa") &&
    identical(as.integer(evidence$bootstrap_index), 584L) &&
    identical(
      tolower(evidence$checkpoint_md5),
      expected_md5_v2[["anti_mrsa_b0584_checkpoint"]]
    ) &&
    identical(as.integer(evidence$mice_seed), 20411311L) &&
    identical(evidence$error_message, reviewed_b0584_failure_message_v2) &&
    identical(as.integer(evidence$event_n), 46L) &&
    is.data.frame(evidence$events) &&
    nrow(evidence$events) == 46L &&
    all(evidence$events$dep == "race_white") &&
    all(evidence$events$meth == "racewhitelogreg") &&
    all(!evidence$allowed) &&
    is_exact_b0584_event_summary_v3(evidence$summary) &&
    identical(
      tolower(evidence$worker_md5),
      expected_md5_v2[["anti_mrsa_b0584_worker_before"]]
    )
}

is_verified_b0584_preflight_v2 <- function(path) {
  if (!file.exists(path)) return(FALSE)
  if (!identical(
    tolower(unname(tools::md5sum(path))),
    expected_md5_v2[["anti_mrsa_b0584_preflight"]]
  )) return(FALSE)
  evidence <- readRDS(path)
  audit <- evidence$audit
  structural <- evidence$structural
  is.list(evidence) &&
    identical(evidence$status, "success") &&
    identical(evidence$diagnostic_only, TRUE) &&
    identical(evidence$candidate_strategy, race_white_lindep_strategy_v3) &&
    identical(evidence$trial, "anti_mrsa") &&
    identical(as.integer(evidence$bootstrap_index), 584L) &&
    identical(
      tolower(evidence$checkpoint_md5),
      expected_md5_v2[["anti_mrsa_b0584_checkpoint"]]
    ) &&
    is.data.frame(evidence$result) &&
    nrow(evidence$result) == length(expected_specifications) &&
    all(vapply(
      evidence$result[c(
        "risk_aki_0", "risk_aki_1", "rd_aki", "rr_aki",
        "risk_death30_0", "risk_death30_1", "rd_death30", "rr_death30"
      )],
      function(value) all(is.finite(value)),
      logical(1)
    )) &&
    is.list(audit) &&
    identical(audit$status, "success") &&
    isTRUE(audit$all_logged_events_allowed) &&
    nrow(audit$mice_logged_events) == 46L &&
    is_exact_b0584_event_summary_v3(audit$mice_logged_event_summary) &&
    is.list(structural) &&
    all(structural$race_allowed) &&
    all(!structural$existing_allowed) &&
    identical(structural$race_class, "factor") &&
    identical(structural$race_levels, c("0", "1")) &&
    identical(structural$antibiotic_class, c("ordered", "factor")) &&
    identical(structural$antibiotic_levels, c("0", "1", "2")) &&
    identical(
      tolower(evidence$worker_md5),
      expected_md5_v2[["anti_mrsa_b0584_worker_before"]]
    )
}

if (!is_reviewed_b0584_failure_v2(reviewed_b0584_failure_path_v2)) {
  stop("anti-MRSA b=584原始V2失败证据不匹配。")
}
if (!is_verified_b0584_diagnostic_v2(reviewed_b0584_diagnostic_path_v2)) {
  stop("anti-MRSA b=584只读MICE事件诊断证据不匹配。")
}
if (!is_verified_b0584_preflight_v2(reviewed_b0584_preflight_path_v2)) {
  stop("anti-MRSA b=584候选精确白名单预检未通过。")
}

# 登记 b=635 原始失败和同编号同种子事件路由诊断；全部永久保留。
reviewed_b0635_failure_path_v2 <- file.path(
  failure_dir_v2,
  "anti_mrsa_b0635_sensitivity_controlled_v2_failure.rds"
)
reviewed_b0635_diagnostic_script_path_v2 <- file.path(
  package_root_v2,
  "07_敏感性分析",
  "10j_diagnose_st08_b0635_event_routing_controlled_v2.R"
)
reviewed_b0635_diagnostic_path_v2 <- file.path(
  sensitivity_dir,
  "diagnostic_history",
  "st08_event_routing_anti_mrsa_b0635",
  "diagnostic_evidence_20260805_211027.rds"
)
reviewed_b0635_failure_message_v2 <- paste0(
  "NEW_MS03 MICE未产生完整MIDS或出现未批准loggedEvents；事件数=5000"
)

is_reviewed_b0635_failure_v2 <- function(path) {
  if (!file.exists(path)) return(FALSE)
  if (!identical(
    tolower(unname(tools::md5sum(path))),
    expected_md5_v2[["anti_mrsa_b0635_v2_failure"]]
  )) return(FALSE)
  evidence <- readRDS(path)
  identical(evidence$status, "failed") &&
    identical(evidence$execution_version, "controlled_v2") &&
    identical(evidence$trial, "anti_mrsa") &&
    identical(as.integer(evidence$bootstrap_index), 635L) &&
    identical(
      tolower(evidence$checkpoint_md5),
      expected_md5_v2[["anti_mrsa_b0635_checkpoint"]]
    ) &&
    identical(evidence$message, reviewed_b0635_failure_message_v2) &&
    identical(
      tolower(evidence$locked_script_md5),
      expected_md5_v2[["locked_sensitivity"]]
    ) &&
    identical(
      tolower(evidence$controlled_v2_script_md5),
      expected_md5_v2[["anti_mrsa_b0635_worker_before"]]
    )
}

is_exact_b0635_route_summary_v4 <- function(summary) {
  if (!is.data.frame(summary) || nrow(summary) != 5L) return(FALSE)
  summary <- summary[order(summary$dep, summary$meth), , drop = FALSE]
  identical(
    sort(summary$dep),
    sort(c(
      "urine_output_rate_6h_ml_kg_h",
      "urine_output_rate_12h_ml_kg_h",
      "urine_output_rate_24h_ml_kg_h",
      "other_broad_spectrum_antibiotics",
      "race_white"
    ))
  ) &&
    identical(
      sort(paste(summary$dep, summary$meth, sep = "|")),
      sort(c(
        "urine_output_rate_6h_ml_kg_h|pmm",
        "urine_output_rate_12h_ml_kg_h|pmm",
        "urine_output_rate_24h_ml_kg_h|pmm",
        "other_broad_spectrum_antibiotics|polr",
        "race_white|racewhitelogreg"
      ))
    ) &&
    all(summary$out == "aki_competing_statedeath_before_aki") &&
    all(summary$route == "existing_v2") &&
    all(as.integer(summary$event_n) == 1000L) &&
    all(as.integer(summary$existing_allowed_n) == 1000L) &&
    all(as.integer(summary$race_lindep_allowed_n) == 0L) &&
    all(as.integer(summary$corrected_allowed_n) == 1000L)
}

is_verified_b0635_diagnostic_v2 <- function(path) {
  if (!file.exists(path)) return(FALSE)
  if (!identical(
    tolower(unname(tools::md5sum(path))),
    expected_md5_v2[["anti_mrsa_b0635_diagnostic"]]
  )) return(FALSE)
  if (!identical(
    tolower(unname(tools::md5sum(
      reviewed_b0635_diagnostic_script_path_v2
    ))),
    expected_md5_v2[["anti_mrsa_b0635_diagnostic_script"]]
  )) return(FALSE)
  evidence <- readRDS(path)
  is.list(evidence) &&
    identical(evidence$status, "diagnostic_complete") &&
    identical(evidence$diagnostic_only, TRUE) &&
    identical(evidence$trial, "anti_mrsa") &&
    identical(as.integer(evidence$bootstrap_index), 635L) &&
    identical(
      tolower(evidence$checkpoint_md5_before),
      expected_md5_v2[["anti_mrsa_b0635_checkpoint"]]
    ) &&
    identical(
      tolower(evidence$checkpoint_md5_after),
      expected_md5_v2[["anti_mrsa_b0635_checkpoint"]]
    ) &&
    identical(
      tolower(evidence$failure_md5_before),
      expected_md5_v2[["anti_mrsa_b0635_v2_failure"]]
    ) &&
    identical(
      tolower(evidence$failure_md5_after),
      expected_md5_v2[["anti_mrsa_b0635_v2_failure"]]
    ) &&
    identical(as.integer(evidence$event_n), 5000L) &&
    is.data.frame(evidence$events) &&
    nrow(evidence$events) == 5000L &&
    all(evidence$original_existing_allowed) &&
    identical(as.integer(evidence$corrected_race_lindep_event_n), 0L) &&
    all(evidence$corrected_allowed) &&
    all(evidence$route == "existing_v2") &&
    is_exact_b0635_route_summary_v4(evidence$route_summary) &&
    isTRUE(evidence$shared_structure_verified) &&
    is.data.frame(evidence$result) &&
    nrow(evidence$result) == length(expected_specifications) &&
    isTRUE(evidence$all_key_values_finite) &&
    identical(evidence$formal_repeat_written, FALSE) &&
    identical(evidence$formal_audit_written, FALSE) &&
    identical(
      tolower(evidence$worker_md5),
      expected_md5_v2[["anti_mrsa_b0635_worker_before"]]
    ) &&
    identical(
      tolower(evidence$diagnostic_script_md5),
      expected_md5_v2[["anti_mrsa_b0635_diagnostic_script"]]
    )
}

if (!is_reviewed_b0635_failure_v2(reviewed_b0635_failure_path_v2)) {
  stop("anti-MRSA b=635原始V2失败证据不匹配。")
}
if (!is_verified_b0635_diagnostic_v2(reviewed_b0635_diagnostic_path_v2)) {
  stop("anti-MRSA b=635事件路由只读诊断证据未通过。")
}

# 仅用于ST-08 NEW_MS03的nephrotoxic logreg100最终层救援审计。
.logreg100_escalation_audit_v2 <- new.env(parent = emptyenv())
.logreg100_escalation_audit_v2$rows <- list()

reset_logreg100_escalation_audit_v2 <- function() {
  .logreg100_escalation_audit_v2$rows <- list()
  invisible(NULL)
}

empty_logreg100_escalation_audit_v2 <- function() {
  data.frame(
    retry_strategy = character(0),
    target_field = character(0),
    locked_method = character(0),
    call_index = integer(0),
    initial_maxit = integer(0),
    initial_iterations = integer(0),
    initial_warning = character(0),
    initial_converged = logical(0),
    first_retry_maxit = integer(0),
    first_retry_iterations = integer(0),
    first_retry_warning = character(0),
    first_retry_converged = logical(0),
    final_retry_maxit = integer(0),
    final_retry_iterations = integer(0),
    final_retry_warning = character(0),
    final_retry_converged = logical(0),
    observed_n = integer(0),
    missing_n = integer(0),
    observed_event_n = integer(0),
    coefficient_all_finite = logical(0),
    covariance_all_finite = logical(0),
    covariance_chol_ok = logical(0),
    probability_all_finite = logical(0),
    probability_min = numeric(0),
    probability_max = numeric(0),
    rescue_verified = logical(0),
    stringsAsFactors = FALSE
  )
}

get_logreg100_escalation_audit_v2 <- function() {
  if (length(.logreg100_escalation_audit_v2$rows) == 0L) {
    return(empty_logreg100_escalation_audit_v2())
  }
  do.call(rbind, .logreg100_escalation_audit_v2$rows)
}

# 覆盖仅限当前ST-08进程；全局MICE扩展文件及主Bootstrap均保持原样。
mice.impute.logreg100 <- function(y, ry, x, wy = NULL, ...) {
  .logreg100_retry_audit_environment$call_n <-
    .logreg100_retry_audit_environment$call_n + 1L
  current_call_index <- .logreg100_retry_audit_environment$call_n
  if (is.null(wy)) wy <- !ry

  # 完全复用mice官方数据增广和锁定的quasibinomial-logit模型。
  augmented <- mice:::augment(y, ry, x, wy)
  x <- cbind(1, as.matrix(augmented$x))
  y <- augmented$y
  ry <- augmented$ry
  wy <- augmented$wy
  weight <- augmented$w

  fit_with_warning_capture_v2 <- function(internal_maxit) {
    captured_warnings <- character(0)
    fit <- withCallingHandlers(
      stats::glm.fit(
        x = x[ry, , drop = FALSE],
        y = y[ry],
        family = stats::quasibinomial(link = "logit"),
        weights = weight[ry],
        control = stats::glm.control(maxit = internal_maxit)
      ),
      warning = function(condition) {
        captured_warnings <<- c(
          captured_warnings,
          conditionMessage(condition)
        )
        invokeRestart("muffleWarning")
      }
    )
    list(fit = fit, warnings = captured_warnings)
  }

  initial <- fit_with_warning_capture_v2(100L)
  fit <- initial$fit
  exact_initial_nonconvergence <- (
    length(initial$warnings) == 1L &&
      identical(
        initial$warnings,
        "glm.fit: algorithm did not converge"
      ) &&
      !isTRUE(fit$converged)
  )

  if (exact_initial_nonconvergence) {
    initial_iterations <- as.integer(fit$iter)
    first_retry <- fit_with_warning_capture_v2(2000L)
    first_retry_ok <- (
      length(first_retry$warnings) == 0L &&
        isTRUE(first_retry$fit$converged) &&
        all(is.finite(stats::coef(first_retry$fit)))
    )
    if (first_retry_ok) {
      # 保持既有V1审计结构，历史和新增普通成功调用完全兼容。
      .logreg100_retry_audit_environment$rows[[
        length(.logreg100_retry_audit_environment$rows) + 1L
      ]] <- data.frame(
        target_field = "nephrotoxic_drug_exposure",
        locked_method = "logreg",
        call_index = current_call_index,
        initial_maxit = 100L,
        initial_iterations = initial_iterations,
        retry_maxit = 2000L,
        retry_iterations = as.integer(first_retry$fit$iter),
        observed_n = sum(ry),
        missing_n = sum(wy),
        observed_event_n = sum(as.numeric(y[ry]) == 2L),
        stringsAsFactors = FALSE
      )
      fit <- first_retry$fit
    } else {
      exact_first_retry_nonconvergence <- (
        length(first_retry$warnings) == 1L &&
          identical(
            first_retry$warnings,
            "glm.fit: algorithm did not converge"
          ) &&
          !isTRUE(first_retry$fit$converged)
      )
      if (!exact_first_retry_nonconvergence) {
        stop(
          "logreg100第一层重试出现未批准异常；警告=",
          paste(unique(first_retry$warnings), collapse = " | "),
          "；converged=",
          isTRUE(first_retry$fit$converged)
        )
      }

      # 只有两层均为同一精确非收敛时，才允许一次20000次最终救援。
      final_retry <- fit_with_warning_capture_v2(20000L)
      final_beta <- stats::coef(final_retry$fit)
      final_covariance <- tryCatch(
        summary.glm(final_retry$fit)$cov.unscaled,
        error = function(condition) NULL
      )
      coefficient_all_finite <- all(is.finite(final_beta))
      covariance_all_finite <- (
        !is.null(final_covariance) &&
          length(final_covariance) > 0L &&
          all(is.finite(final_covariance))
      )
      covariance_chol_ok <- (
        covariance_all_finite &&
          !inherits(
            try(
              chol(mice:::sym(final_covariance)),
              silent = TRUE
            ),
            "try-error"
          )
      )
      final_probability <- if (coefficient_all_finite) {
        stats::plogis(drop(x[wy, , drop = FALSE] %*% final_beta))
      } else {
        rep(NA_real_, sum(wy))
      }
      probability_all_finite <- (
        length(final_probability) > 0L &&
          all(is.finite(final_probability))
      )
      rescue_verified <- (
        length(final_retry$warnings) == 0L &&
          isTRUE(final_retry$fit$converged) &&
          coefficient_all_finite &&
          covariance_all_finite &&
          covariance_chol_ok &&
          probability_all_finite
      )
      .logreg100_escalation_audit_v2$rows[[
        length(.logreg100_escalation_audit_v2$rows) + 1L
      ]] <- data.frame(
        retry_strategy = "nephrotoxic_logreg100_escalation_verified_v2",
        target_field = "nephrotoxic_drug_exposure",
        locked_method = "logreg",
        call_index = current_call_index,
        initial_maxit = 100L,
        initial_iterations = initial_iterations,
        initial_warning = paste(initial$warnings, collapse = " | "),
        initial_converged = isTRUE(initial$fit$converged),
        first_retry_maxit = 2000L,
        first_retry_iterations = as.integer(first_retry$fit$iter),
        first_retry_warning = paste(
          first_retry$warnings,
          collapse = " | "
        ),
        first_retry_converged = isTRUE(first_retry$fit$converged),
        final_retry_maxit = 20000L,
        final_retry_iterations = as.integer(final_retry$fit$iter),
        final_retry_warning = paste(final_retry$warnings, collapse = " | "),
        final_retry_converged = isTRUE(final_retry$fit$converged),
        observed_n = sum(ry),
        missing_n = sum(wy),
        observed_event_n = sum(as.numeric(y[ry]) == 2L),
        coefficient_all_finite = coefficient_all_finite,
        covariance_all_finite = covariance_all_finite,
        covariance_chol_ok = covariance_chol_ok,
        probability_all_finite = probability_all_finite,
        probability_min = if (probability_all_finite) {
          min(final_probability)
        } else {
          NA_real_
        },
        probability_max = if (probability_all_finite) {
          max(final_probability)
        } else {
          NA_real_
        },
        rescue_verified = rescue_verified,
        stringsAsFactors = FALSE
      )
      if (!rescue_verified) {
        stop(
          "logreg100最终层受控救援失败；警告=",
          paste(unique(final_retry$warnings), collapse = " | "),
          "；converged=",
          isTRUE(final_retry$fit$converged)
        )
      }

      # 扁平既有审计记录最终成功层；详细的2000失败保存在V2专用审计。
      .logreg100_retry_audit_environment$rows[[
        length(.logreg100_retry_audit_environment$rows) + 1L
      ]] <- data.frame(
        target_field = "nephrotoxic_drug_exposure",
        locked_method = "logreg",
        call_index = current_call_index,
        initial_maxit = 100L,
        initial_iterations = initial_iterations,
        retry_maxit = 20000L,
        retry_iterations = as.integer(final_retry$fit$iter),
        observed_n = sum(ry),
        missing_n = sum(wy),
        observed_event_n = sum(as.numeric(y[ry]) == 2L),
        stringsAsFactors = FALSE
      )
      fit <- final_retry$fit
    }
  } else if (
    length(initial$warnings) > 0L ||
      !isTRUE(fit$converged) ||
      any(!is.finite(stats::coef(fit)))
  ) {
    stop(
      "logreg100内部glm.fit出现未批准异常；警告=",
      paste(unique(initial$warnings), collapse = " | "),
      "；converged=",
      isTRUE(fit$converged)
    )
  }

  # 官方mice.impute.logreg的后验系数随机抽取与Bernoulli插补保持不变。
  fit_summary <- summary.glm(fit)
  beta <- stats::coef(fit)
  random_variation <- t(chol(mice:::sym(fit_summary$cov.unscaled)))
  beta_star <- beta + random_variation %*% stats::rnorm(
    ncol(random_variation)
  )
  probability <- stats::plogis(
    drop(x[wy, , drop = FALSE] %*% beta_star)
  )
  if (any(!is.finite(probability))) {
    stop("logreg100最终插补概率出现非有限值。")
  }
  imputed <- stats::runif(length(probability)) <= probability
  imputed[imputed] <- 1
  if (is.factor(y)) {
    imputed <- factor(imputed, c(0, 1), levels(y))
  }
  imputed
}

# 定义race_white条件重试审计环境。
# CBPS-ATE 单次条件式重试审计环境。
.cbps_retry_audit_v2 <- new.env(parent = emptyenv())
.cbps_retry_audit_v2$call_n <- 0L
.cbps_retry_audit_v2$rows <- list()

reset_cbps_retry_audit_v2 <- function() {
  .cbps_retry_audit_v2$call_n <- 0L
  .cbps_retry_audit_v2$rows <- list()
  invisible(NULL)
}

empty_cbps_retry_audit_v2 <- function() {
  data.frame(
    retry_strategy = character(0),
    trial = character(0),
    bootstrap_index = integer(0),
    specification = character(0),
    imputation_index = integer(0),
    call_index = integer(0),
    locked_method = character(0),
    estimand = character(0),
    formula = character(0),
    initial_maxit = integer(0),
    initial_warning = character(0),
    initial_convergence_code = integer(0),
    initial_objective = numeric(0),
    initial_function_evaluations = integer(0),
    initial_gradient_evaluations = integer(0),
    retry_maxit = integer(0),
    retry_warning = character(0),
    retry_convergence_code = integer(0),
    retry_objective = numeric(0),
    retry_function_evaluations = integer(0),
    retry_gradient_evaluations = integer(0),
    ps_min = numeric(0),
    ps_max = numeric(0),
    weight_min = numeric(0),
    weight_max = numeric(0),
    ps_lt_1e8 = integer(0),
    weight_gt10 = integer(0),
    weight_gt20 = integer(0),
    ess_continuation = numeric(0),
    ess_deescalation = numeric(0),
    max_absolute_covariate_smd = numeric(0),
    max_smd_term = character(0),
    unbalanced_covariate_n = integer(0),
    retry_verified = logical(0),
    stringsAsFactors = FALSE
  )
}

get_cbps_retry_audit_v2 <- function() {
  if (length(.cbps_retry_audit_v2$rows) == 0L) {
    return(empty_cbps_retry_audit_v2())
  }
  do.call(rbind, .cbps_retry_audit_v2$rows)
}

find_call_frame_value_v2 <- function(name, default = NULL) {
  frames <- rev(sys.frames())
  for (frame in frames) {
    if (exists(name, envir = frame, inherits = FALSE)) {
      return(get(name, envir = frame, inherits = FALSE))
    }
  }
  default
}

extract_optim_metric_v2 <- function(fit, name) {
  if (identical(name, "objective")) {
    value <- fit$obj$value
  } else {
    counts <- fit$obj$counts
    value <- if (is.null(counts)) NULL else counts[[name]]
  }
  if (is.null(value) || length(value) != 1L || !is.finite(value)) {
    return(if (identical(name, "objective")) NA_real_ else NA_integer_)
  }
  if (identical(name, "objective")) as.numeric(value) else as.integer(value)
}

run_cbps_weightit_v2 <- function(data, formula, maxit = NULL) {
  warnings <- character(0)
  arguments <- list(
    formula = formula,
    data = data,
    method = "cbps",
    estimand = "ATE",
    stabilize = FALSE,
    over = FALSE,
    include.obj = TRUE
  )
  if (!is.null(maxit)) arguments$maxit <- as.integer(maxit)
  fit <- withCallingHandlers(
    do.call(WeightIt::weightit, arguments),
    warning = function(condition) {
      warnings <<- c(warnings, conditionMessage(condition))
      invokeRestart("muffleWarning")
    }
  )
  list(
    fit = fit,
    warnings = warnings,
    convergence_code = extract_cbps_convergence_boot(fit),
    objective = extract_optim_metric_v2(fit, "objective"),
    function_evaluations = extract_optim_metric_v2(fit, "function"),
    gradient_evaluations = extract_optim_metric_v2(fit, "gradient")
  )
}

fit_bootstrap_weights_locked_v2 <- fit_bootstrap_weights
fit_bootstrap_weights <- function(data, formula, model) {
  if (!identical(model, "cbps_ate")) {
    return(fit_bootstrap_weights_locked_v2(data, formula, model))
  }

  .cbps_retry_audit_v2$call_n <- .cbps_retry_audit_v2$call_n + 1L
  call_index <- .cbps_retry_audit_v2$call_n
  treatment <- as.integer(as.character(data$deescalation))
  if (anyNA(treatment) || !setequal(unique(treatment), c(0L, 1L))) {
    stop("CBPS 数据缺少任一策略组。")
  }

  exact_warning <- paste0(
    "The optimization failed to converge; try again with a higher value of ",
    "`maxit`."
  )
  initial <- run_cbps_weightit_v2(data, formula)
  initial_converged <- identical(initial$convergence_code, 0L)
  retry_triggered <- (
    length(initial$warnings) == 1L &&
      identical(initial$warnings[[1]], exact_warning) &&
      !is.na(initial$convergence_code) &&
      !initial_converged
  )

  if (length(initial$warnings) == 0L && initial_converged) {
    selected <- initial
    retry <- NULL
  } else if (retry_triggered) {
    context <- list(
      trial = find_call_frame_value_v2("trial"),
      bootstrap_index = find_call_frame_value_v2("bootstrap_index"),
      specification = find_call_frame_value_v2("specification"),
      imputation_index = find_call_frame_value_v2("imputation_index")
    )
    if (
      is.null(context$trial) ||
        is.null(context$bootstrap_index) ||
        is.null(context$specification) ||
        is.null(context$imputation_index)
    ) {
      stop("CBPS 重试无法确定 trial/bootstrap/specification/imputation 上下文。")
    }
    retry <- run_cbps_weightit_v2(data, formula, maxit = 20000L)
    if (
      length(retry$warnings) != 0L ||
        !identical(retry$convergence_code, 0L)
    ) {
      stop(
        "CBPS maxit=20000 单次重试失败；convergence=",
        retry$convergence_code,
        "；warnings=",
        paste(retry$warnings, collapse = " | ")
      )
    }
    selected <- retry
  } else {
    stop(
      "CBPS 出现未批准的警告或收敛状态；convergence=",
      initial$convergence_code,
      "；warnings=",
      paste(initial$warnings, collapse = " | ")
    )
  }

  fit <- selected$fit
  propensity_score <- as.numeric(fit$ps)
  raw_weight <- as.numeric(fit$weights)
  marginal_treated <- mean(treatment == 1L)
  weight <- ifelse(
    treatment == 1L,
    raw_weight * marginal_treated,
    raw_weight * (1 - marginal_treated)
  )
  if (
    length(propensity_score) != nrow(data) ||
      any(!is.finite(propensity_score)) ||
      any(propensity_score <= 0 | propensity_score >= 1) ||
      length(weight) != nrow(data) ||
      any(!is.finite(weight)) ||
      any(weight <= 0)
  ) {
    stop("CBPS 收敛后 PS 或权重 QA 失败。")
  }

  balance <- cobalt::bal.tab(
    formula,
    data = data,
    weights = weight,
    estimand = "ATE",
    s.d.denom = "pooled",
    binary = "std",
    continuous = "std",
    un = TRUE
  )$Balance
  balance <- balance[balance$Type != "Distance", , drop = FALSE]
  absolute_smd <- abs(balance$Diff.Adj)
  if (
    nrow(balance) == 0L ||
      length(absolute_smd) == 0L ||
      any(!is.finite(absolute_smd))
  ) {
    stop("CBPS 协变量平衡 QA 无法计算。")
  }
  maximum_smd_index <- which.max(absolute_smd)
  ess <- vapply(
    c(0L, 1L),
    function(strategy) {
      effective_sample_size_boot(weight[treatment == strategy])
    },
    numeric(1)
  )

  if (retry_triggered) {
    unbalanced_n <- sum(absolute_smd >= 0.1)
    if (unbalanced_n > 0L) {
      stop(
        "CBPS 单次重试后存在绝对 SMD >= 0.1 的协变量；n=",
        unbalanced_n
      )
    }
    .cbps_retry_audit_v2$rows[[
      length(.cbps_retry_audit_v2$rows) + 1L
    ]] <- data.frame(
      retry_strategy = "cbps_ate_maxit_retry_verified_v1",
      trial = as.character(context$trial),
      bootstrap_index = as.integer(context$bootstrap_index),
      specification = as.character(context$specification),
      imputation_index = as.integer(context$imputation_index),
      call_index = as.integer(call_index),
      locked_method = "cbps",
      estimand = "ATE",
      formula = paste(deparse(formula), collapse = " "),
      initial_maxit = 5000L,
      initial_warning = initial$warnings[[1]],
      initial_convergence_code = as.integer(initial$convergence_code),
      initial_objective = initial$objective,
      initial_function_evaluations = initial$function_evaluations,
      initial_gradient_evaluations = initial$gradient_evaluations,
      retry_maxit = 20000L,
      retry_warning = "",
      retry_convergence_code = as.integer(retry$convergence_code),
      retry_objective = retry$objective,
      retry_function_evaluations = retry$function_evaluations,
      retry_gradient_evaluations = retry$gradient_evaluations,
      ps_min = min(propensity_score),
      ps_max = max(propensity_score),
      weight_min = min(weight),
      weight_max = max(weight),
      ps_lt_1e8 = sum(propensity_score < 1e-8),
      weight_gt10 = sum(weight > 10),
      weight_gt20 = sum(weight > 20),
      ess_continuation = unname(ess[1]),
      ess_deescalation = unname(ess[2]),
      max_absolute_covariate_smd = absolute_smd[maximum_smd_index],
      max_smd_term = rownames(balance)[maximum_smd_index],
      unbalanced_covariate_n = as.integer(unbalanced_n),
      retry_verified = TRUE,
      stringsAsFactors = FALSE
    )
  }

  list(
    weight = weight,
    propensity_score = propensity_score,
    max_absolute_smd = absolute_smd[maximum_smd_index],
    max_absolute_smd_term = rownames(balance)[maximum_smd_index],
    ess_continuation = unname(ess[1]),
    ess_deescalation = unname(ess[2]),
    max_weight = max(weight),
    warning_n = 0L
  )
}

.race_white_logreg_audit_v2 <- new.env(parent = emptyenv())
.race_white_logreg_audit_v2$call_n <- 0L
.race_white_logreg_audit_v2$rows <- list()

reset_race_white_logreg_audit_v2 <- function() {
  .race_white_logreg_audit_v2$call_n <- 0L
  .race_white_logreg_audit_v2$rows <- list()
  invisible(NULL)
}

empty_race_white_retry_audit_v2 <- function() {
  data.frame(
    target_field = character(0),
    locked_method = character(0),
    call_index = integer(0),
    initial_maxit = integer(0),
    initial_iterations = integer(0),
    initial_warning = character(0),
    initial_converged = logical(0),
    retry_maxit = integer(0),
    retry_iterations = integer(0),
    retry_warning = character(0),
    retry_converged = logical(0),
    observed_n = integer(0),
    missing_n = integer(0),
    observed_event_n = integer(0),
    stringsAsFactors = FALSE
  )
}

empty_existing_logreg_retry_audit_v2 <- function() {
  data.frame(
    target_field = character(0),
    locked_method = character(0),
    call_index = integer(0),
    initial_maxit = integer(0),
    initial_iterations = integer(0),
    initial_warning = character(0),
    initial_converged = logical(0),
    retry_maxit = integer(0),
    retry_iterations = integer(0),
    retry_warning = character(0),
    retry_converged = logical(0),
    observed_n = integer(0),
    missing_n = integer(0),
    observed_event_n = integer(0),
    stringsAsFactors = FALSE
  )
}

get_race_white_logreg_audit_v2 <- function() {
  if (length(.race_white_logreg_audit_v2$rows) == 0L) {
    return(empty_race_white_retry_audit_v2())
  }
  do.call(rbind, .race_white_logreg_audit_v2$rows)
}

# 官方mice.impute.logreg的目标特异数值收敛适配器。
mice.impute.racewhitelogreg <- function(y, ry, x, wy = NULL, ...) {
  .race_white_logreg_audit_v2$call_n <-
    .race_white_logreg_audit_v2$call_n + 1L
  current_call_index <- .race_white_logreg_audit_v2$call_n

  if (is.null(wy)) {
    wy <- !ry
  }
  augmented <- mice:::augment(y, ry, x, wy)
  x <- augmented$x
  y <- augmented$y
  ry <- augmented$ry
  wy <- augmented$wy
  weight <- augmented$w
  x <- cbind(1, as.matrix(x))

  initial_warnings <- character(0)
  fit <- withCallingHandlers(
    stats::glm.fit(
      x = x[ry, , drop = FALSE],
      y = y[ry],
      family = stats::quasibinomial(link = "logit"),
      weights = weight[ry],
      control = stats::glm.control(maxit = 25)
    ),
    warning = function(condition) {
      initial_warnings <<- c(
        initial_warnings,
        conditionMessage(condition)
      )
      invokeRestart("muffleWarning")
    }
  )

  exact_nonconvergence <- (
    length(initial_warnings) == 1L &&
      identical(
        initial_warnings,
        "glm.fit: algorithm did not converge"
      ) &&
      !isTRUE(fit$converged)
  )
  if (exact_nonconvergence) {
    retry_warnings <- character(0)
    retry_fit <- withCallingHandlers(
      stats::glm.fit(
        x = x[ry, , drop = FALSE],
        y = y[ry],
        family = stats::quasibinomial(link = "logit"),
        weights = weight[ry],
        control = stats::glm.control(maxit = 2000)
      ),
      warning = function(condition) {
        retry_warnings <<- c(
          retry_warnings,
          conditionMessage(condition)
        )
        invokeRestart("muffleWarning")
      }
    )
    retry_ok <- (
      length(retry_warnings) == 0L &&
        isTRUE(retry_fit$converged) &&
        all(is.finite(stats::coef(retry_fit)))
    )
    .race_white_logreg_audit_v2$rows[[
      length(.race_white_logreg_audit_v2$rows) + 1L
    ]] <- data.frame(
      target_field = "race_white",
      locked_method = "logreg",
      call_index = current_call_index,
      initial_maxit = 25L,
      initial_iterations = as.integer(fit$iter),
      initial_warning = paste(initial_warnings, collapse = " | "),
      initial_converged = isTRUE(fit$converged),
      retry_maxit = 2000L,
      retry_iterations = as.integer(retry_fit$iter),
      retry_warning = paste(retry_warnings, collapse = " | "),
      retry_converged = isTRUE(retry_fit$converged),
      observed_n = sum(ry),
      missing_n = sum(wy),
      observed_event_n = sum(as.numeric(y[ry]) == 2L),
      stringsAsFactors = FALSE
    )
    if (!retry_ok) {
      stop(
        "race_white内部glm.fit条件式重试失败；警告=",
        paste(unique(retry_warnings), collapse = " | "),
        "；converged=",
        isTRUE(retry_fit$converged)
      )
    }
    fit <- retry_fit
  } else if (
    length(initial_warnings) > 0L ||
      !isTRUE(fit$converged) ||
      any(!is.finite(stats::coef(fit)))
  ) {
    stop(
      "race_white内部glm.fit出现未批准异常；警告=",
      paste(unique(initial_warnings), collapse = " | "),
      "；converged=",
      isTRUE(fit$converged)
    )
  }

  fit_summary <- summary.glm(fit)
  beta <- stats::coef(fit)
  random_variation <- t(chol(mice:::sym(fit_summary$cov.unscaled)))
  beta_star <- beta + random_variation %*% stats::rnorm(
    ncol(random_variation)
  )
  probability <- 1 / (
    1 + exp(-(x[wy, , drop = FALSE] %*% beta_star))
  )
  if (any(!is.finite(probability))) {
    stop("race_white条件式logreg产生非有限插补概率。")
  }
  imputed <- stats::runif(nrow(probability)) <= probability
  imputed[imputed] <- 1
  if (is.factor(y)) {
    imputed <- factor(imputed, c(0, 1), levels(y))
  }
  imputed
}

# 保存当前NEW_MS03 MICE审计，供单重复结果绑定。
.new_ms03_audit_v2 <- new.env(parent = emptyenv())
.new_ms03_audit_v2$current <- NULL

empty_mice_logged_events_v2 <- function() {
  data.frame(
    it = integer(0),
    im = integer(0),
    dep = character(0),
    meth = character(0),
    out = character(0),
    stringsAsFactors = FALSE
  )
}

summarize_logged_events_v2 <- function(events) {
  if (nrow(events) == 0L) {
    return(data.frame(
      dep = character(0),
      meth = character(0),
      out = character(0),
      event_n = integer(0),
      stringsAsFactors = FALSE
    ))
  }
  stats::aggregate(
    list(event_n = rep.int(1L, nrow(events))),
    by = events[c("dep", "meth", "out")],
    FUN = sum
  )
}

audit_logged_events_v2 <- function(events, work_data, method) {
  if (nrow(events) == 0L) {
    return(logical(0))
  }
  allowed_ordered_contrast <- (
    events$dep == "urine_output_rate_24h_ml_kg_h" &
      events$meth == "pmm" &
      grepl(
        paste0(
          "^other_broad_spectrum_antibiotics\\.(L|Q)",
          "(, other_broad_spectrum_antibiotics\\.(L|Q))*$"
        ),
        events$out
      )
  )
  allowed_zero_count_death <- vapply(
    seq_len(nrow(events)),
    function(event_index) {
      dependent_field <- events$dep[event_index]
      if (
        !dependent_field %in% names(work_data) ||
          !dependent_field %in% names(method) ||
          events$out[event_index] !=
            "aki_competing_statedeath_before_aki" ||
          events$meth[event_index] !=
            unname(method[dependent_field])
      ) {
        return(FALSE)
      }
      observed_target <- !is.na(work_data[[dependent_field]])
      death_before_aki <- (
        as.character(work_data$aki_competing_state) == "death_before_aki"
      )
      sum(observed_target & death_before_aki, na.rm = TRUE) == 0L
    },
    logical(1)
  )
  allowed_ordered_contrast | allowed_zero_count_death
}

empty_race_white_lindep_summary_v3 <- function() {
  data.frame(
    trial = character(0),
    bootstrap_index = integer(0),
    strategy = character(0),
    it = integer(0),
    im = integer(0),
    dep = character(0),
    meth = character(0),
    out = character(0),
    locked_method = character(0),
    runtime_method = character(0),
    race_white_class = character(0),
    race_white_levels = character(0),
    ordered_antibiotic_class = character(0),
    ordered_antibiotic_levels = character(0),
    predictor_included = logical(0),
    exact_out = logical(0),
    event_verified = logical(0),
    stringsAsFactors = FALSE
  )
}

empty_mice_event_route_summary_v4 <- function() {
  data.frame(
    dep = character(0),
    meth = character(0),
    out = character(0),
    route = character(0),
    event_n = integer(0),
    existing_allowed_n = integer(0),
    race_lindep_allowed_n = integer(0),
    corrected_allowed_n = integer(0),
    stringsAsFactors = FALSE
  )
}

summarize_mice_event_routes_v4 <- function(events, route) {
  if (nrow(events) == 0L) return(empty_mice_event_route_summary_v4())
  if (length(route) != nrow(events)) {
    stop("MICE事件数与路由向量长度不一致。")
  }
  stats::aggregate(
    list(
      event_n = rep.int(1L, nrow(events)),
      existing_allowed_n = as.integer(route %in% c("existing_v2", "both")),
      race_lindep_allowed_n = as.integer(
        route %in% c("race_white_lindep_v3", "both")
      ),
      corrected_allowed_n = as.integer(route != "disallowed")
    ),
    by = c(events[c("dep", "meth", "out")], list(route = route)),
    FUN = sum
  )
}

# 在既有两类白名单之外，仅按逐事件、逐结构条件审核race_white有序对比列。
audit_logged_events_controlled_v3 <- function(
  events,
  work_data,
  method,
  predictor_matrix,
  trial,
  bootstrap_index,
  race_white_locked_method
) {
  existing_allowed <- audit_logged_events_v2(events, work_data, method)
  if (nrow(events) == 0L) {
    return(list(
      allowed = logical(0),
      existing_allowed = logical(0),
      race_lindep_allowed = logical(0),
      route = character(0),
      route_summary = empty_mice_event_route_summary_v4(),
      race_white_lindep_summary = empty_race_white_lindep_summary_v3(),
      race_white_lindep_event_n = 0L,
      all_race_white_lindep_events_verified = TRUE,
      race_white_locked_method = race_white_locked_method,
      race_white_runtime_method = unname(method["race_white"]),
      race_white_factor_levels = levels(work_data$race_white),
      ordered_antibiotic_factor_levels =
        levels(work_data$other_broad_spectrum_antibiotics),
      ordered_antibiotic_predictor_included = FALSE
    ))
  }

  required_fields_present <- all(c(
    "race_white",
    "other_broad_spectrum_antibiotics"
  ) %in% names(work_data))
  matrix_fields_present <-
    "race_white" %in% rownames(predictor_matrix) &&
    "other_broad_spectrum_antibiotics" %in% colnames(predictor_matrix)
  predictor_included <- isTRUE(matrix_fields_present) &&
    identical(
      as.numeric(predictor_matrix[
        "race_white",
        "other_broad_spectrum_antibiotics"
      ]),
      1
    )

  race_levels <- if (required_fields_present) {
    levels(work_data$race_white)
  } else {
    character(0)
  }
  antibiotic_levels <- if (required_fields_present) {
    levels(work_data$other_broad_spectrum_antibiotics)
  } else {
    character(0)
  }
  race_class <- if (required_fields_present) {
    class(work_data$race_white)
  } else {
    character(0)
  }
  antibiotic_class <- if (required_fields_present) {
    class(work_data$other_broad_spectrum_antibiotics)
  } else {
    character(0)
  }

  exact_out <- grepl(
    paste0(
      "^other_broad_spectrum_antibiotics\\.(L|Q)",
      "(, other_broad_spectrum_antibiotics\\.(L|Q))*$"
    ),
    events$out
  )
  is_race_lindep_event <- events$dep == "race_white" &
    events$meth == "racewhitelogreg" &
    exact_out
  shared_structure_verified <-
    trial %in% c("anti_mrsa", "anti_psa") &&
    bootstrap_index %in% formal_bootstrap_indices[[trial]] &&
    identical(race_white_locked_method, "logreg") &&
    identical(unname(method["race_white"]), "racewhitelogreg") &&
    identical(race_class, "factor") &&
    identical(race_levels, c("0", "1")) &&
    identical(antibiotic_class, c("ordered", "factor")) &&
    identical(antibiotic_levels, c("0", "1", "2")) &&
    isTRUE(predictor_included)
  race_lindep_allowed <-
    is_race_lindep_event & shared_structure_verified

  race_indices <- which(is_race_lindep_event)
  race_summary <- if (length(race_indices) == 0L) {
    empty_race_white_lindep_summary_v3()
  } else {
    data.frame(
      trial = rep(trial, length(race_indices)),
      bootstrap_index = rep(as.integer(bootstrap_index), length(race_indices)),
      strategy = rep(race_white_lindep_strategy_v3, length(race_indices)),
      it = as.integer(events$it[race_indices]),
      im = as.integer(events$im[race_indices]),
      dep = events$dep[race_indices],
      meth = events$meth[race_indices],
      out = events$out[race_indices],
      locked_method = rep(race_white_locked_method, length(race_indices)),
      runtime_method = rep(
        unname(method["race_white"]),
        length(race_indices)
      ),
      race_white_class = rep(
        paste(race_class, collapse = ","),
        length(race_indices)
      ),
      race_white_levels = rep(
        paste(race_levels, collapse = ","),
        length(race_indices)
      ),
      ordered_antibiotic_class = rep(
        paste(antibiotic_class, collapse = ","),
        length(race_indices)
      ),
      ordered_antibiotic_levels = rep(
        paste(antibiotic_levels, collapse = ","),
        length(race_indices)
      ),
      predictor_included = rep(
        isTRUE(predictor_included),
        length(race_indices)
      ),
      exact_out = exact_out[race_indices],
      event_verified = race_lindep_allowed[race_indices],
      stringsAsFactors = FALSE
    )
  }

  route <- ifelse(
    existing_allowed & race_lindep_allowed,
    "both",
    ifelse(
      existing_allowed,
      "existing_v2",
      ifelse(race_lindep_allowed, "race_white_lindep_v3", "disallowed")
    )
  )

  list(
    allowed = existing_allowed | race_lindep_allowed,
    existing_allowed = existing_allowed,
    race_lindep_allowed = race_lindep_allowed,
    route = route,
    route_summary = summarize_mice_event_routes_v4(events, route),
    race_white_lindep_summary = race_summary,
    race_white_lindep_event_n = as.integer(nrow(race_summary)),
    all_race_white_lindep_events_verified = (
      nrow(race_summary) == 0L || all(race_summary$event_verified)
    ),
    race_white_locked_method = race_white_locked_method,
    race_white_runtime_method = unname(method["race_white"]),
    race_white_factor_levels = race_levels,
    ordered_antibiotic_factor_levels = antibiotic_levels,
    ordered_antibiotic_predictor_included = isTRUE(predictor_included)
  )
}

run_new_ms03_mice_controlled_v2 <- function(
  work_data,
  method,
  predictor_matrix,
  visit_sequence,
  seed,
  trial,
  bootstrap_index
) {
  race_white_locked_method <- unname(method["race_white"])
  if (!identical(race_white_locked_method, "logreg")) {
    stop("NEW_MS03 race_white锁定方法不是logreg。")
  }
  method["race_white"] <- "racewhitelogreg"

  if (identical(trial, "anti_psa")) {
    if (!identical(
      unname(method["source_control_status"]),
      "logreg"
    )) {
      stop("NEW_MS03 anti-PSA source_control_status锁定方法不是logreg。")
    }
    method["source_control_status"] <- "sourcecontrollogreg"
    if (!identical(
      unname(method["nephrotoxic_drug_exposure"]),
      "logreg100"
    )) {
      stop("NEW_MS03 anti-PSA nephrotoxic_drug_exposure活动方法不是logreg100。")
    }
  }

  if (
    !exists("reset_logreg100_retry_audit", mode = "function") ||
      !exists("get_logreg100_retry_audit", mode = "function")
  ) {
    stop("缺少既有logreg条件式重试审计函数。")
  }
  reset_logreg100_retry_audit()
  reset_logreg100_escalation_audit_v2()
  reset_race_white_logreg_audit_v2()
  .new_ms03_audit_v2$current <- NULL

  captured_warnings <- character(0)
  # 与下面 mice(seed = seed) 使用同一个脚本派生锁定种子；显式设置仅用于
  # 可审计地声明自定义 logreg 方法中的 rnorm/runif 由该 MICE 随机流控制，
  # mice() 随后会再次设置同一数值，因此不改变既有随机序列。
  set.seed(as.integer(seed))
  imputation <- withCallingHandlers(
    mice::mice(
      work_data,
      m = mice_m,
      maxit = mice_maxit,
      method = method,
      predictorMatrix = predictor_matrix,
      visitSequence = visit_sequence,
      donors = mice_donors,
      seed = seed,
      printFlag = FALSE
    ),
    warning = function(condition) {
      warning_message <- conditionMessage(condition)
      if (grepl(
        "^Number of logged events: [0-9]+$",
        warning_message
      )) {
        captured_warnings <<- c(captured_warnings, warning_message)
        invokeRestart("muffleWarning")
      } else {
        stop(
          "NEW_MS03 MICE出现未分类警告：",
          warning_message,
          call. = FALSE
        )
      }
    }
  )

  events <- imputation$loggedEvents
  if (is.null(events)) {
    events <- empty_mice_logged_events_v2()
  }
  event_audit <- audit_logged_events_controlled_v3(
    events = events,
    work_data = work_data,
    method = method,
    predictor_matrix = predictor_matrix,
    trial = trial,
    bootstrap_index = bootstrap_index,
    race_white_locked_method = race_white_locked_method
  )
  allowed_event <- event_audit$allowed
  expected_logged_event_warning <- if (nrow(events) == 0L) {
    character(0)
  } else {
    paste0("Number of logged events: ", nrow(events))
  }
  if (
    imputation$m != mice_m ||
      !identical(captured_warnings, expected_logged_event_warning) ||
      (nrow(events) > 0L && !all(allowed_event)) ||
      !isTRUE(event_audit$all_race_white_lindep_events_verified)
  ) {
    stop(
      "NEW_MS03 MICE未产生完整MIDS或出现未批准loggedEvents；事件数=",
      nrow(events)
    )
  }

  race_retry <- get_race_white_logreg_audit_v2()
  existing_retry <- get_logreg100_retry_audit()
  logreg100_escalation <- get_logreg100_escalation_audit_v2()
  .new_ms03_audit_v2$current <- list(
    status = "success",
    trial = trial,
    bootstrap_index = as.integer(bootstrap_index),
    specification = "NEW_MS03",
    mice_seed = as.integer(seed),
    mice_m = as.integer(mice_m),
    mice_maxit = as.integer(mice_maxit),
    mice_donors = as.integer(mice_donors),
    race_white_total_call_n =
      as.integer(.race_white_logreg_audit_v2$call_n),
    race_white_retry_summary = race_retry,
    existing_logreg_retry_summary = existing_retry,
    logreg100_escalation_strategy =
      "nephrotoxic_logreg100_escalation_verified_v2",
    logreg100_escalation_n = as.integer(nrow(logreg100_escalation)),
    logreg100_escalation_summary = logreg100_escalation,
    all_logreg100_escalations_verified = (
      nrow(logreg100_escalation) == 0L ||
        all(logreg100_escalation$rescue_verified)
    ),
    mice_warning_messages = captured_warnings,
    mice_logged_events = events,
    mice_logged_event_summary = summarize_logged_events_v2(events),
    mice_logged_event_allowed = allowed_event,
    mice_logged_event_existing_allowed = event_audit$existing_allowed,
    mice_logged_event_race_lindep_allowed =
      event_audit$race_lindep_allowed,
    mice_logged_event_route = event_audit$route,
    mice_logged_event_route_summary = event_audit$route_summary,
    all_logged_events_allowed = (
      nrow(events) == 0L || all(allowed_event)
    ),
    race_white_lindep_strategy = race_white_lindep_strategy_v3,
    race_white_lindep_event_n = event_audit$race_white_lindep_event_n,
    race_white_lindep_summary = event_audit$race_white_lindep_summary,
    all_race_white_lindep_events_verified =
      event_audit$all_race_white_lindep_events_verified,
    race_white_locked_method = event_audit$race_white_locked_method,
    race_white_runtime_method = event_audit$race_white_runtime_method,
    race_white_factor_levels = event_audit$race_white_factor_levels,
    ordered_antibiotic_factor_levels =
      event_audit$ordered_antibiotic_factor_levels,
    ordered_antibiotic_predictor_included =
      event_audit$ordered_antibiotic_predictor_included,
    timestamp = format(Sys.time(), tz = "Asia/Shanghai", usetz = TRUE)
  )
  imputation
}

# 从锁定函数文本中只替换NEW_MS03的MICE调用和笼统事件停止块。
estimate_start_v2 <- which(
  trimws(locked_lines) == "estimate_specification <- function("
)
estimate_end_marker_v2 <- which(
  trimws(locked_lines) ==
    "# 定义需要由主Bootstrap检查点或同一重抽样键重新计算的预设规格。"
)
if (
  length(estimate_start_v2) != 1L ||
    length(estimate_end_marker_v2) != 1L ||
    estimate_start_v2 >= estimate_end_marker_v2
) {
  stop("锁定estimate_specification函数边界不符合预期。")
}
estimate_lines_v2 <- locked_lines[
  estimate_start_v2:(estimate_end_marker_v2 - 1L)
]
mice_call_start_v2 <- which(
  trimws(estimate_lines_v2) == "analysis_imputation <- mice::mice("
)
logged_stop_start_v2 <- which(
  trimws(estimate_lines_v2) ==
    "if (!is.null(analysis_imputation$loggedEvents)) {"
)
if (
  length(mice_call_start_v2) != 1L ||
    length(logged_stop_start_v2) != 1L ||
    logged_stop_start_v2 <= mice_call_start_v2
) {
  stop("锁定NEW_MS03 MICE代码块不符合预期。")
}
logged_stop_end_v2 <- logged_stop_start_v2 + 2L
if (!identical(trimws(estimate_lines_v2[logged_stop_end_v2]), "}")) {
  stop("锁定NEW_MS03 loggedEvents停止块终点不符合预期。")
}
replacement_mice_block_v2 <- c(
  "    analysis_imputation <- run_new_ms03_mice_controlled_v2(",
  "      work_data = work_data,",
  "      method = method,",
  "      predictor_matrix = predictor_matrix,",
  "      visit_sequence = visit_sequence,",
  "      seed = sensitivity_seed,",
  "      trial = trial,",
  paste0(
    "      bootstrap_index = if (!is.null(checkpoint$bootstrap_index)) ",
    "as.integer(checkpoint$bootstrap_index) else 0L"
  ),
  "    )"
)
estimate_lines_v2 <- c(
  estimate_lines_v2[seq_len(mice_call_start_v2 - 1L)],
  replacement_mice_block_v2,
  estimate_lines_v2[
    (logged_stop_end_v2 + 1L):length(estimate_lines_v2)
  ]
)
eval(parse(text = estimate_lines_v2), envir = .GlobalEnv)

# 创建空的零RR边界审计表，保证有无边界时字段结构一致。
empty_zero_rr_boundary_summary_v2 <- function() {
  data.frame(
    trial = character(0),
    bootstrap_index = integer(0),
    specification = character(0),
    outcome = character(0),
    risk_0 = numeric(0),
    risk_1 = numeric(0),
    rd = numeric(0),
    rr = numeric(0),
    denominator_positive = logical(0),
    risk_in_range = logical(0),
    rd_identity_error = numeric(0),
    rr_identity_error = numeric(0),
    boundary_verified = logical(0),
    stringsAsFactors = FALSE
  )
}

zero_rr_absolute_tolerance_v2 <- 1e-12
zero_rr_relative_tolerance_v2 <- 1e-10
zero_rr_boundary_strategy_v2 <- "zero_risk_ratio_boundary_verified_v1"
repeat_result_validation_version_v2 <-
  "strict_risk_rd_rr_identity_zero_boundary_v1"

# 判断数值是否在预设绝对/相对容差内一致。
is_close_numeric_v2 <- function(observed, expected) {
  abs(observed - expected) <=
    zero_rr_absolute_tolerance_v2 +
      zero_rr_relative_tolerance_v2 * abs(expected)
}

# 从通过结构检查的单重复结果提取零RR边界行及其验证证据。
build_zero_rr_boundary_summary_v2 <- function(
  result,
  trial,
  bootstrap_index
) {
  outcome_definitions <- list(
    aki_7d = c(
      risk_0 = "risk_aki_0",
      risk_1 = "risk_aki_1",
      rd = "rd_aki",
      rr = "rr_aki"
    ),
    death_30d = c(
      risk_0 = "risk_death30_0",
      risk_1 = "risk_death30_1",
      rd = "rd_death30",
      rr = "rr_death30"
    )
  )
  rows <- list()
  for (outcome_name in names(outcome_definitions)) {
    fields <- outcome_definitions[[outcome_name]]
    zero_index <- which(result[[fields[["rr"]]]] == 0)
    if (length(zero_index) == 0L) {
      next
    }
    for (row_index in zero_index) {
      risk_0 <- result[[fields[["risk_0"]]]][row_index]
      risk_1 <- result[[fields[["risk_1"]]]][row_index]
      rd <- result[[fields[["rd"]]]][row_index]
      rr <- result[[fields[["rr"]]]][row_index]
      expected_rd <- risk_1 - risk_0
      expected_rr <- risk_1 / risk_0
      denominator_positive <- is.finite(risk_0) && risk_0 > 0
      risk_in_range <- (
        is.finite(risk_0) &&
          is.finite(risk_1) &&
          risk_0 >= 0 && risk_0 <= 1 &&
          risk_1 >= 0 && risk_1 <= 1
      )
      rd_identity_error <- abs(rd - expected_rd)
      rr_identity_error <- abs(rr - expected_rr)
      boundary_verified <- (
        denominator_positive &&
          risk_in_range &&
          identical(risk_1, 0) &&
          identical(rr, 0) &&
          is_close_numeric_v2(rd, expected_rd) &&
          is_close_numeric_v2(rr, expected_rr)
      )
      rows[[length(rows) + 1L]] <- data.frame(
        trial = trial,
        bootstrap_index = as.integer(bootstrap_index),
        specification = as.character(result$specification[row_index]),
        outcome = outcome_name,
        risk_0 = risk_0,
        risk_1 = risk_1,
        rd = rd,
        rr = rr,
        denominator_positive = denominator_positive,
        risk_in_range = risk_in_range,
        rd_identity_error = rd_identity_error,
        rr_identity_error = rr_identity_error,
        boundary_verified = boundary_verified,
        stringsAsFactors = FALSE
      )
    }
  }
  if (length(rows) == 0L) {
    return(empty_zero_rr_boundary_summary_v2())
  }
  do.call(rbind, rows)
}

# 覆盖受控V1验证边界；不修改10b文件及任何风险或效应估计值。
validate_repeat_result <- function(result, trial, bootstrap_index) {
  required_columns <- c(
    "trial",
    "bootstrap_index",
    "specification",
    effect_fields_controlled,
    diagnostic_fields_controlled
  )
  stopifnot(
    is.data.frame(result),
    nrow(result) == length(expected_specifications),
    all(required_columns %in% names(result)),
    all(result$trial == trial),
    all(as.integer(result$bootstrap_index) == bootstrap_index),
    setequal(result$specification, expected_specifications),
    !anyDuplicated(result$specification)
  )
  numeric_fields <- c(effect_fields_controlled, diagnostic_fields_controlled)
  if (!all(vapply(result[numeric_fields], is.numeric, logical(1)))) {
    stop("单重复敏感性结果的数值字段类型不正确。")
  }
  if (!all(vapply(
    result[numeric_fields],
    function(x) all(is.finite(x)),
    logical(1)
  ))) {
    stop("单重复敏感性结果包含非有限值。")
  }

  risk_fields <- c(
    "risk_aki_0",
    "risk_aki_1",
    "risk_death30_0",
    "risk_death30_1"
  )
  if (any(result[risk_fields] < 0 | result[risk_fields] > 1)) {
    stop("单重复敏感性结果包含超出[0,1]的风险。")
  }
  if (any(result$rd_aki < -1 | result$rd_aki > 1 |
      result$rd_death30 < -1 | result$rd_death30 > 1)) {
    stop("单重复敏感性结果包含超出[-1,1]的RD。")
  }
  if (any(result$risk_aki_0 <= 0 | result$risk_death30_0 <= 0)) {
    stop("单重复敏感性结果的RR分母风险不是严格正值。")
  }
  if (any(result$rr_aki < 0 | result$rr_death30 < 0)) {
    stop("单重复敏感性结果包含负RR。")
  }

  expected_rd_aki <- result$risk_aki_1 - result$risk_aki_0
  expected_rd_death30 <-
    result$risk_death30_1 - result$risk_death30_0
  expected_rr_aki <- result$risk_aki_1 / result$risk_aki_0
  expected_rr_death30 <-
    result$risk_death30_1 / result$risk_death30_0
  if (!all(is_close_numeric_v2(result$rd_aki, expected_rd_aki)) ||
      !all(is_close_numeric_v2(
        result$rd_death30,
        expected_rd_death30
      ))) {
    stop("单重复敏感性结果的RD与两策略风险不一致。")
  }
  if (!all(is_close_numeric_v2(result$rr_aki, expected_rr_aki)) ||
      !all(is_close_numeric_v2(
        result$rr_death30,
        expected_rr_death30
      ))) {
    stop("单重复敏感性结果的RR与两策略风险不一致。")
  }

  zero_rr_summary <- build_zero_rr_boundary_summary_v2(
    result,
    trial,
    bootstrap_index
  )
  if (nrow(zero_rr_summary) > 0L &&
      !all(zero_rr_summary$boundary_verified)) {
    stop("单重复敏感性结果包含未通过审计的零RR边界。")
  }
  invisible(zero_rr_summary)
}

# 将NEW_MS03审计绑定到单重复结果属性，保证结果和审计原子同存。
calculate_repeat_result_locked_v1 <- calculate_repeat_result
calculate_repeat_result <- function(
  checkpoint,
  raw_data,
  trial,
  bootstrap_index
) {
  .new_ms03_audit_v2$current <- NULL
  reset_cbps_retry_audit_v2()
  result <- calculate_repeat_result_locked_v1(
    checkpoint,
    raw_data,
    trial,
    bootstrap_index
  )
  zero_rr_boundary_summary <- validate_repeat_result(
    result,
    trial,
    bootstrap_index
  )
  audit <- .new_ms03_audit_v2$current
  if (is.null(audit) || !identical(audit$status, "success")) {
    stop("单重复未产生完整NEW_MS03审计对象。")
  }
  cbps_retry_summary <- get_cbps_retry_audit_v2()
  audit$cbps_retry_strategy <- "cbps_ate_maxit_retry_verified_v1"
  audit$cbps_total_call_n <- as.integer(.cbps_retry_audit_v2$call_n)
  audit$cbps_retry_n <- as.integer(nrow(cbps_retry_summary))
  audit$cbps_retry_summary <- cbps_retry_summary
  audit$all_cbps_retries_verified <- (
    nrow(cbps_retry_summary) == 0L ||
      all(cbps_retry_summary$retry_verified)
  )
  audit$zero_rr_boundary_strategy <- zero_rr_boundary_strategy_v2
  audit$zero_rr_boundary_n <- as.integer(nrow(zero_rr_boundary_summary))
  audit$zero_rr_boundary_summary <- zero_rr_boundary_summary
  audit$all_zero_rr_boundaries_verified <- (
    nrow(zero_rr_boundary_summary) == 0L ||
      all(zero_rr_boundary_summary$boundary_verified)
  )
  audit$repeat_result_validation_version <-
    repeat_result_validation_version_v2
  audit$result_specification_n <- nrow(result)
  audit$result_validation_passed <- TRUE
  audit$locked_sensitivity_script_md5 <-
    unname(tools::md5sum(locked_script_path))
  audit$controlled_v1_script_md5 <-
    unname(tools::md5sum(controlled_v1_path))
  audit$controlled_v2_script_md5 <-
    unname(tools::md5sum(controlled_v2_worker_path))
  audit$mice_extension_md5 <-
    unname(tools::md5sum(mice_extension_path_v2))
  attr(result, "controlled_v2_audit") <- audit
  result
}

build_legacy_repeat_audit_v2 <- function(
  result,
  repeat_path,
  trial,
  bootstrap_index
) {
  list(
    status = "success",
    execution_version = "controlled_v1_reused",
    trial = trial,
    bootstrap_index = as.integer(bootstrap_index),
    result_path = normalizePath(
      repeat_path,
      winslash = "/",
      mustWork = TRUE
    ),
    result_md5 = unname(tools::md5sum(repeat_path)),
    result_specification_n = nrow(result),
    result_validation_passed = TRUE,
    cbps_retry_strategy = "cbps_ate_maxit_retry_verified_v1",
    cbps_total_call_n = NA_integer_,
    cbps_retry_n = 0L,
    cbps_retry_summary = empty_cbps_retry_audit_v2(),
    all_cbps_retries_verified = TRUE,
    race_white_total_call_n = NA_integer_,
    race_white_retry_summary = empty_race_white_retry_audit_v2(),
    existing_logreg_retry_summary = empty_existing_logreg_retry_audit_v2(),
    logreg100_escalation_strategy =
      "nephrotoxic_logreg100_escalation_verified_v2",
    logreg100_escalation_n = 0L,
    logreg100_escalation_summary = empty_logreg100_escalation_audit_v2(),
    all_logreg100_escalations_verified = TRUE,
    zero_rr_boundary_strategy = zero_rr_boundary_strategy_v2,
    zero_rr_boundary_n = 0L,
    zero_rr_boundary_summary = empty_zero_rr_boundary_summary_v2(),
    all_zero_rr_boundaries_verified = TRUE,
    repeat_result_validation_version =
      "legacy_strict_positive_rr_implies_no_zero_boundary",
    mice_warning_messages = character(0),
    mice_logged_events = empty_mice_logged_events_v2(),
    mice_logged_event_summary = summarize_logged_events_v2(
      empty_mice_logged_events_v2()
    ),
    mice_logged_event_allowed = logical(0),
    mice_logged_event_existing_allowed = logical(0),
    mice_logged_event_race_lindep_allowed = logical(0),
    mice_logged_event_route = character(0),
    mice_logged_event_route_summary = empty_mice_event_route_summary_v4(),
    all_logged_events_allowed = TRUE,
    race_white_lindep_strategy = race_white_lindep_strategy_v3,
    race_white_lindep_event_n = 0L,
    race_white_lindep_summary = empty_race_white_lindep_summary_v3(),
    all_race_white_lindep_events_verified = TRUE,
    race_white_locked_method = NA_character_,
    race_white_runtime_method = NA_character_,
    race_white_factor_levels = character(0),
    ordered_antibiotic_factor_levels = character(0),
    ordered_antibiotic_predictor_included = FALSE,
    provenance_note = paste0(
      "10b已在任一glm.fit警告或NEW_MS03 logged event时硬停止；",
      "本对象仅补充既有成功RDS的执行来源，不重新计算。"
    ),
    timestamp = format(Sys.time(), tz = "Asia/Shanghai", usetz = TRUE)
  )
}

validate_repeat_audit_v2 <- function(audit, trial, bootstrap_index) {
  stopifnot(
    is.list(audit),
    identical(audit$status, "success"),
    identical(audit$trial, trial),
    identical(as.integer(audit$bootstrap_index), bootstrap_index),
    isTRUE(audit$result_validation_passed),
    identical(
      as.integer(audit$result_specification_n),
      length(expected_specifications)
    ),
    isTRUE(audit$all_logged_events_allowed)
  )
  mice_logged_events <- if (is.null(audit$mice_logged_events)) {
    empty_mice_logged_events_v2()
  } else {
    audit$mice_logged_events
  }
  route_field_names <- c(
    "mice_logged_event_existing_allowed",
    "mice_logged_event_race_lindep_allowed",
    "mice_logged_event_route",
    "mice_logged_event_route_summary"
  )
  route_field_present <- vapply(
    route_field_names,
    function(field) !is.null(audit[[field]]),
    logical(1)
  )
  if (any(route_field_present) && !all(route_field_present)) {
    stop("MICE事件路由审计字段不完整。")
  }
  has_route_audit <- all(route_field_present)
  if (has_route_audit) {
    existing_route_allowed <-
      as.logical(audit$mice_logged_event_existing_allowed)
    race_route_allowed <-
      as.logical(audit$mice_logged_event_race_lindep_allowed)
    event_route <- as.character(audit$mice_logged_event_route)
    route_summary <- audit$mice_logged_event_route_summary
    # R 4.4.2 的 ifelse(logical(0), ...) 返回 logical(0)；零事件路由锁定为
    # character(0)，以便与正式审计对象执行严格类型恒等校验。
    expected_route <- if (length(event_route) == 0L) {
      character(0)
    } else {
      ifelse(
        existing_route_allowed & race_route_allowed,
        "both",
        ifelse(
          existing_route_allowed,
          "existing_v2",
          ifelse(race_route_allowed, "race_white_lindep_v3", "disallowed")
        )
      )
    }
    recalculated_route_summary <- summarize_mice_event_routes_v4(
      mice_logged_events,
      event_route
    )
    stopifnot(
      length(existing_route_allowed) == nrow(mice_logged_events),
      length(race_route_allowed) == nrow(mice_logged_events),
      length(event_route) == nrow(mice_logged_events),
      identical(
        as.logical(audit$mice_logged_event_allowed),
        existing_route_allowed | race_route_allowed
      ),
      identical(event_route, expected_route),
      !any(event_route == "disallowed"),
      identical(route_summary, recalculated_route_summary)
    )
  } else {
    existing_route_allowed <- logical(0)
    race_route_allowed <- logical(0)
    event_route <- character(0)
    route_summary <- empty_mice_event_route_summary_v4()
  }
  race_lindep <- if (is.null(audit$race_white_lindep_summary)) {
    legacy_unverified_race_lindep <- nrow(mice_logged_events) > 0L &&
      any(
        mice_logged_events$dep == "race_white" &
          mice_logged_events$meth == "racewhitelogreg" &
          grepl(
            paste0(
              "^other_broad_spectrum_antibiotics\\.(L|Q)",
              "(, other_broad_spectrum_antibiotics\\.(L|Q))*$"
            ),
            mice_logged_events$out
          )
      )
    if (legacy_unverified_race_lindep) {
      stop("旧审计中存在未经V3逐事件核验的race_white logged event。")
    }
    empty_race_white_lindep_summary_v3()
  } else {
    audit$race_white_lindep_summary
  }
  race_lindep_n <- if (is.null(audit$race_white_lindep_event_n)) {
    nrow(race_lindep)
  } else {
    as.integer(audit$race_white_lindep_event_n)
  }
  stopifnot(
    identical(race_lindep_n, as.integer(nrow(race_lindep))),
    is.null(audit$all_race_white_lindep_events_verified) ||
      isTRUE(audit$all_race_white_lindep_events_verified)
  )
  if (nrow(race_lindep) > 0L) {
    stopifnot(
      identical(
        audit$race_white_lindep_strategy,
        race_white_lindep_strategy_v3
      ),
      all(race_lindep$trial == trial),
      all(as.integer(race_lindep$bootstrap_index) == bootstrap_index),
      all(race_lindep$strategy == race_white_lindep_strategy_v3),
      all(race_lindep$dep == "race_white"),
      all(race_lindep$meth == "racewhitelogreg"),
      all(grepl(
        paste0(
          "^other_broad_spectrum_antibiotics\\.(L|Q)",
          "(, other_broad_spectrum_antibiotics\\.(L|Q))*$"
        ),
        race_lindep$out
      )),
      all(race_lindep$locked_method == "logreg"),
      all(race_lindep$runtime_method == "racewhitelogreg"),
      all(race_lindep$race_white_class == "factor"),
      all(race_lindep$race_white_levels == "0,1"),
      all(race_lindep$ordered_antibiotic_class == "ordered,factor"),
      all(race_lindep$ordered_antibiotic_levels == "0,1,2"),
      all(race_lindep$predictor_included),
      all(race_lindep$exact_out),
      all(race_lindep$event_verified),
      sum(
        mice_logged_events$dep == "race_white" &
          mice_logged_events$meth == "racewhitelogreg" &
          grepl(
            paste0(
              "^other_broad_spectrum_antibiotics\\.(L|Q)",
              "(, other_broad_spectrum_antibiotics\\.(L|Q))*$"
            ),
            mice_logged_events$out
          )
      ) == nrow(race_lindep)
    )
  }
  if (identical(trial, "anti_mrsa") && bootstrap_index == 584L) {
    stopifnot(
      nrow(race_lindep) == 46L,
      sum(race_lindep$out ==
        "other_broad_spectrum_antibiotics.L") == 20L,
      sum(race_lindep$out ==
        "other_broad_spectrum_antibiotics.Q") == 26L
    )
  }
  if (identical(trial, "anti_mrsa") && bootstrap_index == 635L) {
    expected_route_keys <- sort(c(
      paste(
        "urine_output_rate_6h_ml_kg_h",
        "pmm",
        "aki_competing_statedeath_before_aki",
        sep = "|"
      ),
      paste(
        "urine_output_rate_12h_ml_kg_h",
        "pmm",
        "aki_competing_statedeath_before_aki",
        sep = "|"
      ),
      paste(
        "urine_output_rate_24h_ml_kg_h",
        "pmm",
        "aki_competing_statedeath_before_aki",
        sep = "|"
      ),
      paste(
        "other_broad_spectrum_antibiotics",
        "polr",
        "aki_competing_statedeath_before_aki",
        sep = "|"
      ),
      paste(
        "race_white",
        "racewhitelogreg",
        "aki_competing_statedeath_before_aki",
        sep = "|"
      )
    ))
    observed_route_keys <- sort(paste(
      route_summary$dep,
      route_summary$meth,
      route_summary$out,
      sep = "|"
    ))
    stopifnot(
      has_route_audit,
      nrow(mice_logged_events) == 5000L,
      nrow(race_lindep) == 0L,
      all(existing_route_allowed),
      !any(race_route_allowed),
      all(event_route == "existing_v2"),
      nrow(route_summary) == 5L,
      all(route_summary$route == "existing_v2"),
      all(as.integer(route_summary$event_n) == 1000L),
      identical(observed_route_keys, expected_route_keys),
      identical(
        audit$mice_warning_messages,
        "Number of logged events: 5000"
      )
    )
  }
  cbps_retry <- if (is.null(audit$cbps_retry_summary)) {
    empty_cbps_retry_audit_v2()
  } else {
    audit$cbps_retry_summary
  }
  cbps_retry_n <- if (is.null(audit$cbps_retry_n)) {
    nrow(cbps_retry)
  } else {
    as.integer(audit$cbps_retry_n)
  }
  stopifnot(
    identical(cbps_retry_n, as.integer(nrow(cbps_retry))),
    is.null(audit$all_cbps_retries_verified) ||
      isTRUE(audit$all_cbps_retries_verified)
  )
  if (nrow(cbps_retry) > 0L) {
    stopifnot(
      all(cbps_retry$retry_strategy ==
        "cbps_ate_maxit_retry_verified_v1"),
      all(cbps_retry$trial == trial),
      all(cbps_retry$bootstrap_index == bootstrap_index),
      all(cbps_retry$specification %in% expected_specifications),
      all(cbps_retry$imputation_index >= 1L &
        cbps_retry$imputation_index <= mice_m),
      all(cbps_retry$locked_method == "cbps"),
      all(cbps_retry$estimand == "ATE"),
      all(nzchar(cbps_retry$formula)),
      all(cbps_retry$initial_maxit == 5000L),
      all(cbps_retry$initial_warning == paste0(
        "The optimization failed to converge; try again with a higher value of ",
        "`maxit`."
      )),
      all(cbps_retry$initial_convergence_code != 0L),
      all(cbps_retry$retry_maxit == 20000L),
      all(cbps_retry$retry_warning == ""),
      all(cbps_retry$retry_convergence_code == 0L),
      all(is.finite(cbps_retry$ps_min)),
      all(is.finite(cbps_retry$ps_max)),
      all(cbps_retry$ps_min > 0 & cbps_retry$ps_max < 1),
      all(is.finite(cbps_retry$weight_min)),
      all(is.finite(cbps_retry$weight_max)),
      all(cbps_retry$weight_min > 0),
      all(is.finite(cbps_retry$ess_continuation)),
      all(is.finite(cbps_retry$ess_deescalation)),
      all(is.finite(cbps_retry$max_absolute_covariate_smd)),
      all(cbps_retry$max_absolute_covariate_smd < 0.1),
      all(cbps_retry$unbalanced_covariate_n == 0L),
      all(cbps_retry$retry_verified)
    )
  }
  if (!is.null(audit$race_white_retry_summary) &&
      nrow(audit$race_white_retry_summary) > 0L) {
    retry <- audit$race_white_retry_summary
    stopifnot(
      all(retry$target_field == "race_white"),
      all(retry$locked_method == "logreg"),
      all(retry$initial_maxit == 25L),
      all(retry$initial_warning ==
        "glm.fit: algorithm did not converge"),
      all(!retry$initial_converged),
      all(retry$retry_maxit == 2000L),
      all(retry$retry_warning == ""),
      all(retry$retry_converged)
    )
  }
  escalation <- if (is.null(audit$logreg100_escalation_summary)) {
    empty_logreg100_escalation_audit_v2()
  } else {
    audit$logreg100_escalation_summary
  }
  escalation_n <- if (is.null(audit$logreg100_escalation_n)) {
    nrow(escalation)
  } else {
    as.integer(audit$logreg100_escalation_n)
  }
  stopifnot(
    identical(escalation_n, as.integer(nrow(escalation))),
    is.null(audit$all_logreg100_escalations_verified) ||
      isTRUE(audit$all_logreg100_escalations_verified)
  )
  if (identical(trial, "anti_psa") && bootstrap_index == 79L) {
    stopifnot(nrow(escalation) >= 1L)
  }
  if (nrow(escalation) > 0L) {
    stopifnot(
      identical(
        audit$logreg100_escalation_strategy,
        "nephrotoxic_logreg100_escalation_verified_v2"
      ),
      all(escalation$retry_strategy ==
        "nephrotoxic_logreg100_escalation_verified_v2"),
      all(escalation$target_field == "nephrotoxic_drug_exposure"),
      all(escalation$locked_method == "logreg"),
      all(escalation$initial_maxit == 100L),
      all(escalation$initial_warning ==
        "glm.fit: algorithm did not converge"),
      all(!escalation$initial_converged),
      all(escalation$first_retry_maxit == 2000L),
      all(escalation$first_retry_warning ==
        "glm.fit: algorithm did not converge"),
      all(!escalation$first_retry_converged),
      all(escalation$final_retry_maxit == 20000L),
      all(escalation$final_retry_warning == ""),
      all(escalation$final_retry_converged),
      all(escalation$coefficient_all_finite),
      all(escalation$covariance_all_finite),
      all(escalation$covariance_chol_ok),
      all(escalation$probability_all_finite),
      all(is.finite(escalation$probability_min)),
      all(is.finite(escalation$probability_max)),
      all(escalation$probability_min >= 0),
      all(escalation$probability_max <= 1),
      all(escalation$rescue_verified)
    )
    existing_retry <- audit$existing_logreg_retry_summary
    escalation_existing <- existing_retry[
      existing_retry$target_field == "nephrotoxic_drug_exposure" &
        existing_retry$retry_maxit == 20000L,
      ,
      drop = FALSE
    ]
    stopifnot(
      nrow(escalation_existing) == nrow(escalation),
      identical(
        sort(as.integer(escalation_existing$call_index)),
        sort(as.integer(escalation$call_index))
      )
    )
  }
  zero_rr_boundary <- if (is.null(audit$zero_rr_boundary_summary)) {
    empty_zero_rr_boundary_summary_v2()
  } else {
    audit$zero_rr_boundary_summary
  }
  zero_rr_boundary_n <- if (is.null(audit$zero_rr_boundary_n)) {
    nrow(zero_rr_boundary)
  } else {
    as.integer(audit$zero_rr_boundary_n)
  }
  stopifnot(
    identical(
      zero_rr_boundary_n,
      as.integer(nrow(zero_rr_boundary))
    ),
    is.null(audit$all_zero_rr_boundaries_verified) ||
      isTRUE(audit$all_zero_rr_boundaries_verified)
  )
  if (nrow(zero_rr_boundary) > 0L) {
    stopifnot(
      identical(
        audit$zero_rr_boundary_strategy,
        zero_rr_boundary_strategy_v2
      ),
      all(zero_rr_boundary$trial == trial),
      all(as.integer(zero_rr_boundary$bootstrap_index) == bootstrap_index),
      all(zero_rr_boundary$specification %in% expected_specifications),
      all(zero_rr_boundary$outcome %in% c("aki_7d", "death_30d")),
      all(is.finite(zero_rr_boundary$risk_0)),
      all(is.finite(zero_rr_boundary$risk_1)),
      all(is.finite(zero_rr_boundary$rd)),
      all(is.finite(zero_rr_boundary$rr)),
      all(zero_rr_boundary$risk_0 > 0),
      all(zero_rr_boundary$risk_1 == 0),
      all(zero_rr_boundary$rr == 0),
      all(zero_rr_boundary$denominator_positive),
      all(zero_rr_boundary$risk_in_range),
      all(zero_rr_boundary$rd_identity_error <=
        zero_rr_absolute_tolerance_v2),
      all(zero_rr_boundary$rr_identity_error <=
        zero_rr_absolute_tolerance_v2),
      all(zero_rr_boundary$boundary_verified)
    )
  }
  if (identical(trial, "anti_mrsa") && bootstrap_index == 373L) {
    stopifnot(
      nrow(zero_rr_boundary) == 1L,
      identical(zero_rr_boundary$specification, "complete_case"),
      identical(zero_rr_boundary$outcome, "death_30d")
    )
  }
  invisible(TRUE)
}

is_reviewed_b0031_failure_v2 <- function(path) {
  if (!file.exists(path)) {
    return(FALSE)
  }
  if (!identical(
    tolower(unname(tools::md5sum(path))),
    expected_md5_v2[["anti_mrsa_b0031_failure"]]
  )) {
    return(FALSE)
  }
  evidence <- readRDS(path)
  identical(evidence$status, "failed") &&
    identical(evidence$trial, "anti_mrsa") &&
    identical(as.integer(evidence$bootstrap_index), 31L) &&
    identical(
      evidence$message,
      "未分类警告：glm.fit: algorithm did not converge"
    ) &&
    identical(
      tolower(evidence$locked_script_md5),
      expected_md5_v2[["locked_sensitivity"]]
    )
}

get_argument_v2 <- function(name, default = NULL) {
  arguments <- commandArgs(trailingOnly = TRUE)
  prefix <- paste0("--", name, "=")
  hit <- arguments[startsWith(arguments, prefix)]
  if (length(hit) == 0L) {
    return(default)
  }
  if (length(hit) != 1L) {
    stop("参数重复：", name)
  }
  sub(prefix, "", hit, fixed = TRUE)
}

ensure_existing_repeat_audits_v2 <- function() {
  for (current_trial in names(formal_bootstrap_indices)) {
    for (current_b in formal_bootstrap_indices[[current_trial]]) {
      repeat_path <- file.path(
        repeat_dir,
        sprintf("%s_b%04d_sensitivity.rds", current_trial, current_b)
      )
      if (!file.exists(repeat_path)) {
        next
      }
      result <- readRDS(repeat_path)
      validate_repeat_result(result, current_trial, current_b)
      audit_path <- file.path(
        repeat_audit_dir_v2,
        sprintf("%s_b%04d_sensitivity_audit.rds", current_trial, current_b)
      )
      if (!file.exists(audit_path)) {
        embedded_audit <- attr(result, "controlled_v2_audit")
        audit <- if (is.null(embedded_audit)) {
          build_legacy_repeat_audit_v2(
            result,
            repeat_path,
            current_trial,
            current_b
          )
        } else {
          embedded_audit$result_path <- normalizePath(
            repeat_path,
            winslash = "/",
            mustWork = TRUE
          )
          embedded_audit$result_md5 <- unname(tools::md5sum(repeat_path))
          embedded_audit$execution_version <- "controlled_v2"
          embedded_audit
        }
        atomic_save_rds_controlled(audit, audit_path)
      }
      validate_repeat_audit_v2(
        readRDS(audit_path),
        current_trial,
        current_b
      )
    }
  }
  invisible(TRUE)
}

run_trial_worker_v2 <- function() {
  trial <- get_argument_v2("trial")
  if (is.null(trial) || !trial %in% names(formal_bootstrap_indices)) {
    stop("必须指定--trial=anti_mrsa或--trial=anti_psa。")
  }
  requested_start <- as.integer(get_argument_v2(
    "b_start",
    min(formal_bootstrap_indices[[trial]])
  ))
  requested_end <- as.integer(get_argument_v2(
    "b_end",
    max(formal_bootstrap_indices[[trial]])
  ))
  if (
    is.na(requested_start) ||
      is.na(requested_end) ||
      requested_start > requested_end
  ) {
    stop("b_start/b_end参数无效。")
  }
  current_indices <- formal_bootstrap_indices[[trial]]
  current_indices <- current_indices[
    current_indices >= requested_start &
      current_indices <= requested_end
  ]
  if (length(current_indices) == 0L) {
    stop("指定范围不包含正式Bootstrap编号。")
  }
  if (identical(trial, "anti_psa") && 934L %in% current_indices) {
    stop("anti-PSA b=934不得进入正式敏感性worker。")
  }

  b0031_md5_before <- unname(tools::md5sum(
    reviewed_b0031_failure_path_v2
  ))
  b0934_md5_before <- unname(tools::md5sum(preserved_failure_path))
  reviewed_cbps_failure_md5_before <- vapply(
    reviewed_cbps_failure_registry_v2$path,
    function(path) unname(tools::md5sum(path)),
    character(1)
  )
  b0079_failure_md5_before <- unname(tools::md5sum(
    reviewed_b0079_failure_path_v2
  ))
  b0079_diagnostic_md5_before <- unname(tools::md5sum(
    reviewed_b0079_diagnostic_path_v2
  ))
  b0373_failure_md5_before <- unname(tools::md5sum(
    reviewed_b0373_failure_path_v2
  ))
  b0373_diagnostic_md5_before <- unname(tools::md5sum(
    reviewed_b0373_diagnostic_path_v2
  ))
  b0584_failure_md5_before <- unname(tools::md5sum(
    reviewed_b0584_failure_path_v2
  ))
  b0584_diagnostic_md5_before <- unname(tools::md5sum(
    reviewed_b0584_diagnostic_path_v2
  ))
  b0584_preflight_md5_before <- unname(tools::md5sum(
    reviewed_b0584_preflight_path_v2
  ))
  b0635_failure_md5_before <- unname(tools::md5sum(
    reviewed_b0635_failure_path_v2
  ))
  b0635_diagnostic_md5_before <- unname(tools::md5sum(
    reviewed_b0635_diagnostic_path_v2
  ))
  b0635_diagnostic_script_md5_before <- unname(tools::md5sum(
    reviewed_b0635_diagnostic_script_path_v2
  ))
  ensure_existing_repeat_audits_v2()

  raw_data <- utils::read.csv(
    file.path(data_dir, paste0("analysis_dataset_", trial, ".csv")),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  raw_data <- augment_creatinine_only_outcome(raw_data, trial)

  for (position in seq_along(current_indices)) {
    bootstrap_index <- current_indices[position]
    repeat_path <- file.path(
      repeat_dir,
      sprintf("%s_b%04d_sensitivity.rds", trial, bootstrap_index)
    )
    audit_path <- file.path(
      repeat_audit_dir_v2,
      sprintf("%s_b%04d_sensitivity_audit.rds", trial, bootstrap_index)
    )
    legacy_failure_path <- file.path(
      failure_dir,
      sprintf("%s_b%04d_sensitivity_failure.rds", trial, bootstrap_index)
    )
    v2_failure_path <- file.path(
      failure_dir_v2,
      sprintf(
        "%s_b%04d_sensitivity_controlled_v2_failure.rds",
        trial,
        bootstrap_index
      )
    )
    if (
      file.exists(v2_failure_path) &&
        !is_reviewed_cbps_failure_v2(
          v2_failure_path,
          trial,
          bootstrap_index
        ) &&
        !(
          identical(trial, "anti_psa") &&
            bootstrap_index == 79L &&
            is_reviewed_b0079_failure_v2(v2_failure_path)
        ) &&
        !(
          identical(trial, "anti_mrsa") &&
            bootstrap_index == 373L &&
            is_reviewed_b0373_failure_v2(v2_failure_path)
        ) &&
        !(
          identical(trial, "anti_mrsa") &&
            bootstrap_index == 584L &&
            is_reviewed_b0584_failure_v2(v2_failure_path)
        ) &&
        !(
          identical(trial, "anti_mrsa") &&
            bootstrap_index == 635L &&
            is_reviewed_b0635_failure_v2(v2_failure_path)
        )
    ) {
      stop("存在尚未审核的V2敏感性失败证据：", v2_failure_path)
    }
    if (file.exists(legacy_failure_path)) {
      reviewed <- (
        identical(trial, "anti_mrsa") &&
          bootstrap_index == 31L &&
          is_reviewed_b0031_failure_v2(legacy_failure_path)
      )
      if (!reviewed) {
        stop("存在尚未审核的历史敏感性失败证据：", legacy_failure_path)
      }
    }

    if (file.exists(repeat_path)) {
      result <- readRDS(repeat_path)
      validate_repeat_result(result, trial, bootstrap_index)
      if (!file.exists(audit_path)) {
        embedded_audit <- attr(result, "controlled_v2_audit")
        if (is.null(embedded_audit)) {
          embedded_audit <- build_legacy_repeat_audit_v2(
            result,
            repeat_path,
            trial,
            bootstrap_index
          )
        }
        embedded_audit$result_path <- normalizePath(
          repeat_path,
          winslash = "/",
          mustWork = TRUE
        )
        embedded_audit$result_md5 <- unname(tools::md5sum(repeat_path))
        atomic_save_rds_controlled(embedded_audit, audit_path)
      }
      validate_repeat_audit_v2(
        readRDS(audit_path),
        trial,
        bootstrap_index
      )
    } else {
      checkpoint_path <- file.path(
        checkpoint_dir,
        sprintf("%s_b%04d_checkpoint.rds", trial, bootstrap_index)
      )
      if (!file.exists(checkpoint_path)) {
        stop("缺少正式Bootstrap checkpoint：", checkpoint_path)
      }
      checkpoint <- readRDS(checkpoint_path)
      stopifnot(
        identical(checkpoint$status, "success"),
        identical(checkpoint$trial, trial),
        identical(as.integer(checkpoint$bootstrap_index), bootstrap_index)
      )

      start_time <- Sys.time()
      calculation <- tryCatch(
        withCallingHandlers(
          calculate_repeat_result(
            checkpoint,
            raw_data,
            trial,
            bootstrap_index
          ),
          warning = function(warning_condition) {
            stop(
              "未分类警告：",
              conditionMessage(warning_condition),
              call. = FALSE
            )
          }
        ),
        error = function(error_condition) error_condition
      )
      if (inherits(calculation, "error")) {
        failure_output_path <- if (file.exists(v2_failure_path)) {
          file.path(
            failure_dir_v2,
            sprintf(
              "%s_b%04d_sensitivity_controlled_v2_retry_failure_%s.rds",
              trial,
              bootstrap_index,
              format(Sys.time(), "%Y%m%d_%H%M%S")
            )
          )
        } else {
          v2_failure_path
        }
        atomic_save_rds_controlled(
          list(
            status = "failed",
            execution_version = "controlled_v2",
            trial = trial,
            bootstrap_index = bootstrap_index,
            checkpoint_path = checkpoint_path,
            checkpoint_md5 = unname(tools::md5sum(checkpoint_path)),
            locked_script_md5 = unname(tools::md5sum(locked_script_path)),
            controlled_v2_script_md5 =
              unname(tools::md5sum(controlled_v2_worker_path)),
            message = conditionMessage(calculation),
            timestamp = format(
              Sys.time(),
              tz = "Asia/Shanghai",
              usetz = TRUE
            )
          ),
          failure_output_path
        )
        stop(
          "V2敏感性分析在",
          trial,
          " b=",
          bootstrap_index,
          "停止：",
          conditionMessage(calculation)
        )
      }
      end_time <- Sys.time()
      validate_repeat_result(calculation, trial, bootstrap_index)
      embedded_audit <- attr(calculation, "controlled_v2_audit")
      if (is.null(embedded_audit)) {
        stop("V2成功结果缺少嵌入审计。")
      }
      embedded_audit$execution_version <- "controlled_v2"
      embedded_audit$checkpoint_path <- normalizePath(
        checkpoint_path,
        winslash = "/",
        mustWork = TRUE
      )
      embedded_audit$checkpoint_md5 <-
        unname(tools::md5sum(checkpoint_path))
      embedded_audit$started_at <- format(
        start_time,
        tz = "Asia/Shanghai",
        usetz = TRUE
      )
      embedded_audit$completed_at <- format(
        end_time,
        tz = "Asia/Shanghai",
        usetz = TRUE
      )
      embedded_audit$elapsed_seconds <- as.numeric(difftime(
        end_time,
        start_time,
        units = "secs"
      ))
      attr(calculation, "controlled_v2_audit") <- embedded_audit

      atomic_save_rds_controlled(calculation, repeat_path)
      embedded_audit$result_path <- normalizePath(
        repeat_path,
        winslash = "/",
        mustWork = TRUE
      )
      embedded_audit$result_md5 <- unname(tools::md5sum(repeat_path))
      atomic_save_rds_controlled(embedded_audit, audit_path)
      validate_repeat_audit_v2(
        readRDS(audit_path),
        trial,
        bootstrap_index
      )
    }

    if (position %% 10L == 0L || position == length(current_indices)) {
      message(
        "CONTROLLED_ST08_V2_PROGRESS trial=",
        trial,
        " completed=",
        position,
        "/",
        length(current_indices),
        " last_b=",
        bootstrap_index
      )
    }
  }

  if (!identical(
    b0031_md5_before,
    unname(tools::md5sum(reviewed_b0031_failure_path_v2))
  )) {
    stop("anti-MRSA b=31原失败证据MD5发生变化。")
  }
  if (!identical(
    b0934_md5_before,
    unname(tools::md5sum(preserved_failure_path))
  )) {
    stop("anti-PSA b=934原失败证据MD5发生变化。")
  }
  reviewed_cbps_failure_md5_after <- vapply(
    reviewed_cbps_failure_registry_v2$path,
    function(path) unname(tools::md5sum(path)),
    character(1)
  )
  if (!identical(
    reviewed_cbps_failure_md5_before,
    reviewed_cbps_failure_md5_after
  )) {
    stop("已审核 CBPS V2 failure 证据 MD5 发生变化。")
  }
  if (!identical(
    b0079_failure_md5_before,
    unname(tools::md5sum(reviewed_b0079_failure_path_v2))
  )) {
    stop("anti-PSA b=79原失败证据MD5发生变化。")
  }
  if (!identical(
    b0079_diagnostic_md5_before,
    unname(tools::md5sum(reviewed_b0079_diagnostic_path_v2))
  )) {
    stop("anti-PSA b=79只读诊断证据MD5发生变化。")
  }
  if (!identical(
    b0373_failure_md5_before,
    unname(tools::md5sum(reviewed_b0373_failure_path_v2))
  )) {
    stop("anti-MRSA b=373原失败证据MD5发生变化。")
  }
  if (!identical(
    b0373_diagnostic_md5_before,
    unname(tools::md5sum(reviewed_b0373_diagnostic_path_v2))
  )) {
    stop("anti-MRSA b=373零RR诊断证据MD5发生变化。")
  }
  if (!identical(
    b0584_failure_md5_before,
    unname(tools::md5sum(reviewed_b0584_failure_path_v2))
  )) {
    stop("anti-MRSA b=584原失败证据MD5发生变化。")
  }
  if (!identical(
    b0584_diagnostic_md5_before,
    unname(tools::md5sum(reviewed_b0584_diagnostic_path_v2))
  )) {
    stop("anti-MRSA b=584只读诊断证据MD5发生变化。")
  }
  if (!identical(
    b0584_preflight_md5_before,
    unname(tools::md5sum(reviewed_b0584_preflight_path_v2))
  )) {
    stop("anti-MRSA b=584候选预检证据MD5发生变化。")
  }
  if (!identical(
    b0635_failure_md5_before,
    unname(tools::md5sum(reviewed_b0635_failure_path_v2))
  )) {
    stop("anti-MRSA b=635原始V2失败证据MD5发生变化。")
  }
  if (!identical(
    b0635_diagnostic_md5_before,
    unname(tools::md5sum(reviewed_b0635_diagnostic_path_v2))
  )) {
    stop("anti-MRSA b=635只读事件路由诊断证据MD5发生变化。")
  }
  if (!identical(
    b0635_diagnostic_script_md5_before,
    unname(tools::md5sum(reviewed_b0635_diagnostic_script_path_v2))
  )) {
    stop("anti-MRSA b=635只读诊断脚本MD5发生变化。")
  }
  if (
    file.exists(forbidden_success_path) ||
      file.exists(file.path(
        repeat_dir,
        "anti_psa_b0934_sensitivity.rds"
      ))
  ) {
    stop("anti-PSA b=934出现禁止的成功对象。")
  }
  message(
    "CONTROLLED_ST08_V2_WORKER_STATUS=success; trial=",
    trial,
    "; range=",
    min(current_indices),
    "-",
    max(current_indices)
  )
  invisible(TRUE)
}

# 被10d作为函数库source时不自动启动worker。
if (!identical(Sys.getenv("ST08_CONTROLLED_V2_LIBRARY_ONLY"), "1")) {
  run_trial_worker_v2()
}
