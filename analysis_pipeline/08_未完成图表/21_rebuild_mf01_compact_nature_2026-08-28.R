# ============================================================
# Analysis: Rebuild MF-01 as a compact six-step cohort flow diagram
# Date: 2026-08-28
# Random seed: 42 (fixed for reproducibility; no stochastic operation is used)
# R: 4.4.2
# Key packages: ggplot2; digest
#
# 代码来源：项目既有 09_cohort_flow_diagram.R 与 plot_bmc_compliance_figures.R。
# 引用/commit：RD-4.6-ANALYSIS-OUTPUT-PLAN-V1.0；当前目录不是 Git 仓库，故无 commit ID。
# 访问日期：2026-08-28。
# 适配说明：仅合并展示步骤并采用低饱和 Nature 风格配色；不改原始数据、模型或锁定估计值。
# 外部实现核对：审阅 ggconsort、ggflowchart、consort 与 flowchart 的公开GitHub实现；
# 本图因双试验并列、共享队列脚注及自定义合并排除节点，保留项目内 ggplot2 实现，未复制外部代码。
# 参考链接：https://github.com/tgerke/ggconsort；https://github.com/nrennie/ggflowchart；
# https://github.com/adayim/consort；https://github.com/bruigtp/flowchart（访问日期：2026-08-28）。
# 验证命令：D:/software/R4.4.2/R-4.4.2/bin/x64/Rscript.exe 02_代码/R/08_未完成图表/21_rebuild_mf01_compact_nature_2026-08-28.R
# ============================================================

options(stringsAsFactors = FALSE, warn = 1)
set.seed(42)

for (pkg in c("ggplot2", "digest")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Package '", pkg, "' is required.")
  }
}

project_dir <- Sys.getenv("TTE_PROJECT_ROOT", unset = ".")
project_dir <- normalizePath(project_dir, winslash = "/", mustWork = TRUE)

source_flow_path <- file.path(
  project_dir, "04_主文图表", "09_tables", "MF-01_cohort_flow_diagram_data.csv"
)
strategy_path <- file.path(
  project_dir, "04_主文图表", "09_tables", "MF-01_final_strategy_counts_QA.csv"
)
if (!file.exists(source_flow_path)) stop("Missing locked flow source: ", source_flow_path)
if (!file.exists(strategy_path)) stop("Missing locked strategy-count source: ", strategy_path)

source_flow <- utils::read.csv(
  source_flow_path, check.names = FALSE, stringsAsFactors = FALSE, fileEncoding = "UTF-8"
)
strategy_qa <- utils::read.csv(
  strategy_path, check.names = FALSE, stringsAsFactors = FALSE, fileEncoding = "UTF-8"
)

required_flow_fields <- c("trial", "step_no", "excluded_from_previous", "remaining_n")
if (!all(required_flow_fields %in% names(source_flow))) {
  stop("Locked flow source is missing required fields.")
}
if (!all(c("trial", "treatment_strategy", "icu_stay_id") %in% names(strategy_qa))) {
  stop("Strategy-count source is missing required fields.")
}

trial_spec <- data.frame(
  trial = c("anti_mrsa", "anti_psa"),
  panel = c("A. Anti-MRSA trial", "B. Antipseudomonal trial"),
  accent = c("#3C5488", "#00A087"),
  pale = c("#E8EEF7", "#E2F3EF"),
  expected_start = c(5329L, 5074L),
  expected_final = c(268L, 372L),
  stringsAsFactors = FALSE
)

fmt_n <- function(x) format(as.integer(x), big.mark = ",", scientific = FALSE, trim = TRUE)
wrap_text <- function(x, width) {
  paragraphs <- strsplit(x, "\n", fixed = TRUE)[[1L]]
  wrapped <- unlist(lapply(paragraphs, strwrap, width = width), use.names = FALSE)
  paste(wrapped, collapse = "\n")
}

get_step <- function(d, step_number) {
  out <- d[d$step_no == step_number, , drop = FALSE]
  if (nrow(out) != 1L) stop("Expected exactly one source row for step ", step_number, ".")
  out
}

compact_rows <- list()
for (i in seq_len(nrow(trial_spec))) {
  spec <- trial_spec[i, , drop = FALSE]
  d <- source_flow[source_flow$trial == spec$trial, , drop = FALSE]
  d <- d[order(d$step_no), , drop = FALSE]
  if (!identical(as.integer(d$step_no), 0:12)) {
    stop("Unexpected source-step sequence for ", spec$trial, ".")
  }

  s <- lapply(0:12, function(step_number) get_step(d, step_number))
  names(s) <- as.character(0:12)

  start_n <- s[["0"]]$remaining_n
  excluded_t0 <- sum(vapply(s[c("1", "2")], function(x) x$excluded_from_previous, numeric(1)))
  after_t0 <- s[["2"]]$remaining_n

  excluded_renal_components <- c(
    esrd_or_chronic_dialysis = s[["3"]]$excluded_from_previous,
    unresolved_dialysis_encounter = s[["4"]]$excluded_from_previous,
    missing_kdigo_assessment = s[["5"]]$excluded_from_previous,
    pre_t0_kdigo_aki = s[["6"]]$excluded_from_previous,
    pre_t0_rrt = s[["7"]]$excluded_from_previous,
    final_urine_output_review = s[["12"]]$excluded_from_previous
  )
  excluded_renal <- sum(excluded_renal_components)
  after_renal <- s[["7"]]$remaining_n - s[["12"]]$excluded_from_previous

  excluded_coverage <- s[["8"]]$excluded_from_previous
  after_coverage <- after_renal - excluded_coverage
  excluded_mdro <- s[["9"]]$excluded_from_previous
  after_mdro <- after_coverage - excluded_mdro
  excluded_strategy <- s[["10"]]$excluded_from_previous
  after_strategy <- after_mdro - excluded_strategy
  excluded_episode <- s[["11"]]$excluded_from_previous
  final_n <- after_strategy - excluded_episode

  strategy_d <- strategy_qa[strategy_qa$trial == spec$trial, , drop = FALSE]
  continuation_n <- strategy_d$icu_stay_id[strategy_d$treatment_strategy == 0]
  deescalation_n <- strategy_d$icu_stay_id[strategy_d$treatment_strategy == 1]
  if (length(continuation_n) != 1L || length(deescalation_n) != 1L) {
    stop("Strategy counts are not unique for ", spec$trial, ".")
  }

  # 审计断言：合并展示必须保持原顺序算术与最终队列不变。
  stopifnot(
    start_n == spec$expected_start,
    start_n - excluded_t0 == after_t0,
    after_t0 - excluded_renal == after_renal,
    after_renal - excluded_coverage == after_coverage,
    after_coverage - excluded_mdro == after_mdro,
    after_mdro - excluded_strategy == after_strategy,
    after_strategy - excluded_episode == final_n,
    final_n == s[["12"]]$remaining_n,
    final_n == spec$expected_final,
    continuation_n + deescalation_n == final_n,
    excluded_renal_components[["unresolved_dialysis_encounter"]] == 0,
    excluded_renal_components[["missing_kdigo_assessment"]] == 0
  )

  renal_detail <- paste0(
    "Excluded, n = ", fmt_n(excluded_renal), "\n",
    "ESKD or chronic dialysis: ", fmt_n(excluded_renal_components[["esrd_or_chronic_dialysis"]]), "\n",
    "Pre-T0 KDIGO AKI: ", fmt_n(excluded_renal_components[["pre_t0_kdigo_aki"]]), "\n",
    "Pre-T0 RRT: ", fmt_n(excluded_renal_components[["pre_t0_rrt"]]), "\n",
    "Additional urine-output AKI on final review: ",
    fmt_n(excluded_renal_components[["final_urine_output_review"]])
  )

  node_text <- c(
    paste0("Eligible first target-antibiotic administration after sepsis\nN = ", fmt_n(start_n)),
    paste0("Alive and observable in hospital at time zero (T0)\nN = ", fmt_n(after_t0)),
    paste0("Met renal eligibility criteria at T0\nN = ", fmt_n(after_renal)),
    paste0("Eligible target coverage administered during hours 48–72\nN = ", fmt_n(after_coverage)),
    paste0("No trial-specific target MDRO positivity by T0\nN = ", fmt_n(after_mdro)),
    paste0("Time-zero strategy classified as de-escalation or continuation\nN = ", fmt_n(after_strategy)),
    paste0(
      "Consistent with the index sepsis treatment episode\nFinal cohort: N = ", fmt_n(final_n), "\n",
      "De-escalation: ", fmt_n(deescalation_n), "; continuation: ", fmt_n(continuation_n)
    )
  )
  exclusion_text <- c(
    NA_character_,
    paste0(
      "Excluded, n = ", fmt_n(excluded_t0),
      "\nNot observable in hospital at T0 (including death or discharge)"
    ),
    renal_detail,
    paste0("Excluded, n = ", fmt_n(excluded_coverage), "\nNo eligible target coverage during hours 48–72"),
    paste0("Excluded, n = ", fmt_n(excluded_mdro), "\nTrial-specific target MDRO positivity by T0"),
    paste0("Excluded, n = ", fmt_n(excluded_strategy), "\nStrategy not classifiable at T0"),
    paste0("Excluded, n = ", fmt_n(excluded_episode), "\nNot consistent with the index sepsis treatment episode")
  )
  source_steps <- c("0", "1–2", "3–7 and 12", "8", "9", "10", "11")
  remaining_n <- c(start_n, after_t0, after_renal, after_coverage, after_mdro, after_strategy, final_n)
  excluded_n <- c(NA_integer_, excluded_t0, excluded_renal, excluded_coverage, excluded_mdro, excluded_strategy, excluded_episode)

  compact_rows[[i]] <- data.frame(
    trial = spec$trial,
    panel = spec$panel,
    display_step = 0:6,
    source_steps = source_steps,
    remaining_n = as.integer(remaining_n),
    excluded_n = as.integer(excluded_n),
    node_text = node_text,
    exclusion_text = exclusion_text,
    accent = spec$accent,
    pale = spec$pale,
    stringsAsFactors = FALSE
  )
}

compact <- do.call(rbind, compact_rows)
compact$y <- rep(c(8.35, 7.15, 5.45, 3.85, 2.75, 1.65, 0.45), times = nrow(trial_spec))
compact$panel <- factor(compact$panel, levels = trial_spec$panel)
compact$node_text_wrapped <- vapply(compact$node_text, wrap_text, character(1), width = 47L)
compact$exclusion_text_wrapped <- vapply(
  compact$exclusion_text,
  function(x) if (is.na(x)) NA_character_ else wrap_text(x, 44L),
  character(1)
)

audit_dir <- file.path(project_dir, "04_主文图表", "09_tables")
dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)
audit_path <- file.path(audit_dir, "MF-01_compact_cohort_flow_diagram_data.csv")
audit_output <- compact[c(
  "trial", "panel", "display_step", "source_steps", "remaining_n", "excluded_n",
  "node_text", "exclusion_text"
)]
audit_output$source_flow_sha256 <- digest::digest(file = source_flow_path, algo = "sha256")
audit_output$strategy_count_sha256 <- digest::digest(file = strategy_path, algo = "sha256")
utils::write.csv(audit_output, audit_path, row.names = FALSE, fileEncoding = "UTF-8")

main_nodes <- compact
exclusion_nodes <- compact[compact$display_step > 0, , drop = FALSE]
vertical_arrows <- do.call(rbind, lapply(split(compact, compact$trial), function(d) {
  d <- d[order(d$display_step), , drop = FALSE]
  data.frame(
    panel = d$panel[-nrow(d)],
    accent = d$accent[-nrow(d)],
    y = d$y[-nrow(d)] - 0.32,
    yend = d$y[-1L] + 0.32,
    stringsAsFactors = FALSE
  )
}))
branch_lines <- data.frame(
  panel = exclusion_nodes$panel,
  accent = exclusion_nodes$accent,
  y = exclusion_nodes$y,
  stringsAsFactors = FALSE
)

fig1 <- ggplot2::ggplot() +
  ggplot2::geom_segment(
    data = vertical_arrows,
    ggplot2::aes(x = 0, xend = 0, y = y, yend = yend, colour = accent),
    linewidth = 0.58,
    arrow = grid::arrow(length = grid::unit(0.065, "inches"), type = "closed")
  ) +
  ggplot2::geom_segment(
    data = branch_lines,
    ggplot2::aes(x = 0.55, xend = 0.94, y = y, yend = y, colour = accent),
    linewidth = 0.48,
    arrow = grid::arrow(length = grid::unit(0.055, "inches"), type = "closed")
  ) +
  ggplot2::geom_label(
    data = main_nodes,
    ggplot2::aes(x = 0, y = y, label = node_text_wrapped, fill = pale, colour = accent),
    family = "Times New Roman",
    fontface = "plain",
    size = 3.15,
    lineheight = 0.96,
    linewidth = 0.65,
    label.r = grid::unit(0.10, "lines"),
    label.padding = grid::unit(0.32, "lines")
  ) +
  ggplot2::geom_label(
    data = exclusion_nodes,
    ggplot2::aes(x = 1.58, y = y, label = exclusion_text_wrapped),
    family = "Times New Roman",
    size = 2.65,
    lineheight = 0.94,
    linewidth = 0.45,
    label.r = grid::unit(0.08, "lines"),
    label.padding = grid::unit(0.26, "lines"),
    fill = "#FAF4F1",
    colour = "#6B4D43"
  ) +
  ggplot2::facet_wrap(~panel, nrow = 1) +
  ggplot2::scale_fill_identity() +
  ggplot2::scale_colour_identity() +
  ggplot2::coord_cartesian(xlim = c(-0.68, 2.36), ylim = c(0.02, 8.82), clip = "off") +
  ggplot2::labs(
    caption = paste0(
      "The two trial cohorts were not mutually exclusive: 171 participants were included in both cohorts. ",
      "The renal-eligibility step includes the final review of pre-time-zero urine-output AKI."
    )
  ) +
  ggplot2::theme_void(base_family = "Times New Roman", base_size = 10) +
  ggplot2::theme(
    strip.text = ggplot2::element_text(face = "bold", size = 12.2, colour = "#252525"),
    plot.caption = ggplot2::element_text(size = 9.2, hjust = 0, colour = "#4A4A4A", margin = ggplot2::margin(t = 8)),
    panel.spacing = grid::unit(0.42, "inches"),
    plot.margin = ggplot2::margin(10, 20, 12, 20)
  )

staging_dir <- file.path(project_dir, "manuscript", "plot_new")
package_dir <- file.path(project_dir, "09_英文投稿统一版", "01_main_figures")
core_dir <- file.path(project_dir, "04_主文图表", "08_figures")
data_dir <- file.path(project_dir, "data")
for (path in c(staging_dir, package_dir, core_dir, data_dir)) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

staging_paths <- c(
  png = file.path(staging_dir, "Figure1_cohort_flow_bmc.png"),
  pdf = file.path(staging_dir, "Figure1_cohort_flow_bmc.pdf"),
  tiff = file.path(staging_dir, "Figure1_cohort_flow_bmc.tiff")
)
ggplot2::ggsave(staging_paths[["png"]], fig1, width = 14, height = 10.2, units = "in", dpi = 320, bg = "white")
ggplot2::ggsave(
  staging_paths[["pdf"]], fig1, width = 14, height = 10.2, units = "in",
  device = grDevices::cairo_pdf, bg = "white"
)
ggplot2::ggsave(
  staging_paths[["tiff"]], fig1, width = 14, height = 10.2, units = "in", dpi = 600,
  device = "tiff", compression = "lzw", bg = "white"
)
if (any(!file.exists(staging_paths)) || any(file.info(staging_paths)$size <= 0)) {
  stop("One or more Figure 1 staging outputs are missing or empty.")
}

copy_targets <- list(
  package = file.path(package_dir, paste0("MF-01_cohort_flow_diagram.", c("png", "pdf", "tiff"))),
  core = file.path(core_dir, paste0("MF-01_cohort_flow_diagram.", c("png", "pdf", "tiff"))),
  data = file.path(data_dir, paste0("Figure 1_Cohort Flow Diagram.", c("png", "pdf", "tiff")))
)
names(staging_paths) <- c("png", "pdf", "tiff")
for (target_group in copy_targets) {
  target_group <- target_group[match(names(staging_paths), tools::file_ext(target_group))]
  copied <- file.copy(unname(staging_paths), unname(target_group), overwrite = TRUE, copy.date = TRUE)
  if (!all(copied)) stop("Failed to synchronise one or more Figure 1 outputs.")
}

hash_rows <- do.call(rbind, lapply(
  c(unname(staging_paths), unlist(copy_targets, use.names = FALSE)),
  function(path) data.frame(
    path = normalizePath(path, winslash = "/", mustWork = TRUE),
    size_bytes = file.info(path)$size,
    sha256 = digest::digest(file = path, algo = "sha256"),
    stringsAsFactors = FALSE
  )
))
hash_audit_path <- file.path(audit_dir, "MF-01_compact_output_sha256.csv")
utils::write.csv(hash_rows, hash_audit_path, row.names = FALSE, fileEncoding = "UTF-8")

message("MF-01 compact flow diagram generated and synchronised.")
message("Compact audit data: ", audit_path)
message("Output hash audit: ", hash_audit_path)
