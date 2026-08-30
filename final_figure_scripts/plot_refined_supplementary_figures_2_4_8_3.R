# =============================================================================
# 分析：Supplementary Figures 2, 4, 8, and 3 的纯展示层重绘
# 日期：2026-08-28
# 随机性：无；本脚本不抽样、不插补、不拟合模型、不运行 bootstrap
# R：4.4.2
#
# 数据来源：
# 1) ST-05_missing_data_long_audit.csv（锁定缺失比例汇总）
# 2) SF-04_SF-05_mice_chain_statistics_*.csv（锁定 MICE 链统计）
# 3) ST-08_sensitivity_percentile_ci.csv 与归档出院处理敏感性结果
# 4) SF-03_continuous_correlation_data.csv（锁定 Spearman 相关系数）
# 版本/引用：RD-4.6-ANALYSIS-OUTPUT-PLAN-V1.0；controlled reporting outputs
# 访问日期：2026-08-28
# 适配说明：仅调整布局、标签、颜色、留白和直接数值标注；所有数值原样读取。
# 验证命令：
# D:/software/R4.4.2/R-4.4.2/bin/x64/Rscript.exe manuscript/plot_refined_supplementary_figures_2_4_8_3.R SF02
# 可选参数：SF02、SF04、SF08、SF03 或 ALL；ALL 仍按 SF02 -> SF04 -> SF08 -> SF03 执行。
# =============================================================================

options(stringsAsFactors = FALSE, warn = 1)
suppressPackageStartupMessages({
  library(ggplot2)
  library(gridExtra)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("Cannot determine the current script path.")
script_path <- normalizePath(sub("^--file=", "", script_arg[[1L]]), winslash = "/", mustWork = TRUE)
project_root <- normalizePath(Sys.getenv("TTE_PROJECT_ROOT", unset = "."), winslash = "/", mustWork = TRUE)
analysis_root <- normalizePath(Sys.getenv("TTE_RUNTIME_ROOT", unset = "."), winslash = "/", mustWork = TRUE)
locked_dir <- file.path(analysis_root, "09_输出结果", "09_tables", "controlled_reporting")
out_dir <- file.path(project_root, "manuscript", "plot_new")
audit_dir <- file.path(project_root, "05_补充图", "09_tables")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)

read_locked <- function(path) {
  if (!file.exists(path)) stop("Missing locked source: ", path)
  out <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE, fileEncoding = "UTF-8")
  if (nrow(out) == 0L) stop("Locked source has no rows: ", path)
  out
}

save_figure <- function(plot_object, stem, width, height) {
  ggsave(file.path(out_dir, paste0(stem, ".png")), plot_object,
         width = width, height = height, units = "in", dpi = 320, bg = "white")
  ggsave(file.path(out_dir, paste0(stem, ".pdf")), plot_object,
         width = width, height = height, units = "in", device = grDevices::cairo_pdf, bg = "white")
  ggsave(file.path(out_dir, paste0(stem, ".tiff")), plot_object,
         width = width, height = height, units = "in", dpi = 600,
         device = "tiff", compression = "lzw", bg = "white")
}

dictionary <- read_locked(file.path(
  analysis_root, "02_配置与变量字典", "publication_display_labels_v1.0.csv"
))
display_by_field <- setNames(dictionary$display_label, dictionary$clean_field)
trial_labels <- c(
  anti_mrsa = "A. Anti-MRSA trial",
  anti_psa = "B. Antipseudomonal trial"
)

selected <- toupper(commandArgs(trailingOnly = TRUE))
if (length(selected) == 0L) selected <- "ALL"
valid <- c("SF02", "SF04", "SF08", "SF03", "ALL")
if (length(selected) != 1L || !selected %in% valid) {
  stop("Select one of: SF02, SF04, SF08, SF03, ALL.")
}
should_run <- function(id) identical(selected, "ALL") || identical(selected, id)

# -----------------------------------------------------------------------------
# SF02：仅显示至少一个策略分层存在缺失的变量，保留所有锁定百分比。
# -----------------------------------------------------------------------------
if (should_run("SF02")) {
  missing_all <- read_locked(file.path(locked_dir, "ST-05_missing_data_long_audit.csv"))
  missing_all$missing_pct <- as.numeric(missing_all$missing_pct)
  missing_strategy <- missing_all[missing_all$stratum %in% c("continued", "deescalated"), , drop = FALSE]

  max_by_field <- aggregate(missing_pct ~ field, missing_strategy, max, na.rm = TRUE)
  nonzero_fields <- max_by_field$field[max_by_field$missing_pct > 0]
  display <- missing_strategy[missing_strategy$field %in% nonzero_fields, , drop = FALSE]
  stopifnot(
    length(unique(missing_all$field)) == 38L,
    length(nonzero_fields) == 10L,
    nrow(display) == 40L
  )

  display$variable_label <- unname(display_by_field[display$field])
  if (anyNA(display$variable_label)) stop("SF02 has an unmapped variable label.")
  field_max <- setNames(max_by_field$missing_pct, max_by_field$field)
  ordered_fields <- names(sort(field_max[nonzero_fields], decreasing = TRUE))
  ordered_labels <- unname(display_by_field[ordered_fields])
  display$variable_label <- factor(display$variable_label, levels = rev(ordered_labels))
  display$trial_label <- factor(display$trial, levels = names(trial_labels), labels = unname(trial_labels))
  display$strategy_label <- factor(
    display$stratum,
    levels = c("continued", "deescalated"),
    labels = c("Continuation", "De-escalation")
  )
  display$text_colour <- ifelse(display$missing_pct >= 18, "white", "#252525")
  omitted_n <- length(unique(missing_all$field)) - length(nonzero_fields)

  sf02 <- ggplot(display, aes(strategy_label, variable_label, fill = missing_pct)) +
    geom_tile(colour = "white", linewidth = 0.8, width = 0.96, height = 0.94) +
    geom_text(
      aes(label = sprintf("%.1f%%", missing_pct), colour = text_colour),
      size = 3.4, family = "serif", fontface = "bold"
    ) +
    facet_wrap(~ trial_label, nrow = 1) +
    scale_colour_identity() +
    scale_fill_gradient(
      low = "#F7FBFF", high = "#2171B5", limits = c(0, 32),
      breaks = c(0, 8, 16, 24, 32), name = "Missing (%)"
    ) +
    labs(
      title = "Missing data by trial and treatment strategy",
      subtitle = "Variables with any nonzero missingness in the displayed strategy strata",
      x = NULL, y = NULL,
      caption = sprintf(
        "%d additional prespecified analysis variables had 0%% missingness in every displayed stratum and are not plotted.",
        omitted_n
      )
    ) +
    theme_minimal(base_size = 10.5, base_family = "serif") +
    theme(
      plot.title = element_text(face = "bold", size = 13.2, margin = margin(b = 3)),
      plot.subtitle = element_text(size = 9.5, colour = "grey30", margin = margin(b = 8)),
      strip.text = element_text(face = "bold", size = 11.2, hjust = 0),
      axis.text.x = element_text(size = 10, colour = "#252525", margin = margin(t = 5)),
      axis.text.y = element_text(size = 9.2, colour = "#252525"),
      panel.grid = element_blank(),
      panel.spacing.x = grid::unit(0.42, "inches"),
      legend.position = "right",
      legend.title = element_text(face = "bold"),
      plot.caption = element_text(size = 8.7, hjust = 0, colour = "grey30", margin = margin(t = 8)),
      plot.margin = margin(10, 12, 9, 10)
    )

  write.csv(display, file.path(audit_dir, "SF-02_missing_data_refined_display_data.csv"), row.names = FALSE)
  save_figure(sf02, "SF02_missing_data_refined", 11.3, 6.2)
  message("SF02 complete: 40 locked strategy-specific cells; 10 displayed variables.")
}

# -----------------------------------------------------------------------------
# SF04：保留全部 20×50 链轨迹，并增加跨链中位数以帮助判断稳定性。
# -----------------------------------------------------------------------------
if (should_run("SF04")) {
  chain <- rbind(
    read_locked(file.path(locked_dir, "SF-04_SF-05_mice_chain_statistics_anti_mrsa.csv")),
    read_locked(file.path(locked_dir, "SF-04_SF-05_mice_chain_statistics_anti_psa.csv"))
  )
  chain$iteration <- as.integer(chain$iteration)
  chain$chain <- as.integer(chain$chain)
  chain$chain_mean <- as.numeric(chain$chain_mean)
  stopifnot(
    nrow(chain) == 13000L,
    length(unique(chain$chain)) == 20L,
    max(chain$iteration, na.rm = TRUE) == 50L
  )
  chain$variable_label <- unname(display_by_field[chain$variable])
  if (anyNA(chain$variable_label)) stop("SF04 has an unmapped variable label.")
  chain$chain_group <- interaction(chain$trial, chain$variable, chain$chain, drop = TRUE)

  median_data <- aggregate(
    chain_mean ~ trial + variable + variable_label + iteration,
    chain,
    median,
    na.rm = TRUE
  )

  make_chain_panel <- function(trial_key, title, show_x_title) {
    d <- chain[chain$trial == trial_key, , drop = FALSE]
    med <- median_data[median_data$trial == trial_key, , drop = FALSE]
    variable_order <- unique(d$variable_label)
    d$variable_label <- factor(d$variable_label, levels = variable_order)
    med$variable_label <- factor(med$variable_label, levels = variable_order)

    ggplot(d, aes(iteration, chain_mean, group = chain_group)) +
      geom_line(colour = "#7E9FB2", alpha = 0.24, linewidth = 0.30, na.rm = TRUE) +
      geom_line(
        data = med, aes(iteration, chain_mean, group = variable_label),
        inherit.aes = FALSE, colour = "#17365D", linewidth = 0.68, na.rm = TRUE
      ) +
      facet_wrap(~ variable_label, ncol = 4, scales = "free_y") +
      scale_x_continuous(breaks = seq(0, 50, by = 10), limits = c(1, 50)) +
      labs(title = title, x = if (show_x_title) "Iteration" else NULL, y = "Chain mean") +
      theme_minimal(base_size = 9.2, base_family = "serif") +
      theme(
        plot.title = element_text(face = "bold", size = 12.2, hjust = 0, margin = margin(b = 5)),
        strip.text = element_text(face = "bold", size = 8.6, lineheight = 0.95),
        axis.text = element_text(size = 7.8, colour = "#3A3A3A"),
        axis.title = element_text(size = 9.5),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = "#E8E8E8", linewidth = 0.28),
        panel.spacing = grid::unit(0.52, "lines"),
        plot.margin = margin(5, 8, 5, 8)
      )
  }

  panel_mrsa <- make_chain_panel("anti_mrsa", "A. Anti-MRSA trial", FALSE)
  panel_psa <- make_chain_panel("anti_psa", "B. Antipseudomonal trial", TRUE)
  sf04 <- arrangeGrob(
    panel_mrsa,
    panel_psa,
    ncol = 1,
    heights = c(1.08, 1.00),
    top = grid::textGrob(
      "MICE chain means across 20 imputations",
      x = 0.01, hjust = 0, gp = grid::gpar(fontfamily = "serif", fontsize = 14, fontface = "bold")
    ),
    bottom = grid::textGrob(
      "Thin lines are individual imputation chains; the dark line is the across-chain median. Discrete traces reflect binary/categorical variables with few missing observations.",
      x = 0.01, hjust = 0, gp = grid::gpar(fontfamily = "serif", fontsize = 9, col = "#4A4A4A")
    )
  )

  write.csv(chain, file.path(audit_dir, "SF-04_mice_chain_mean_refined_display_data.csv"), row.names = FALSE)
  write.csv(median_data, file.path(audit_dir, "SF-04_mice_chain_mean_refined_median_trace.csv"), row.names = FALSE)
  save_figure(sf04, "SF04_mice_chain_mean_refined", 13.6, 10.5)
  message("SF04 complete: 13,000 locked chain rows; 20 chains and 50 iterations retained.")
}

# -----------------------------------------------------------------------------
# SF08：用视觉分组区分八项 CBPS-ATE 与一项 ATO，并增加数值列。
# -----------------------------------------------------------------------------
if (should_run("SF08")) {
  sensitivity_source <- read_locked(file.path(
    analysis_root, "09_输出结果", "06_敏感性分析", "controlled_reporting",
    "ST-08_sensitivity_percentile_ci.csv"
  ))
  core <- sensitivity_source[
    sensitivity_source$estimand_field == "rd_aki" &
      !sensitivity_source$specification %in% c(
        "main_cbps_ate", "overlap_weight_ato_brlogit", "death_calendar_boundary"
      ),
    , drop = FALSE
  ]
  stopifnot(nrow(core) == 16L)
  d <- data.frame(
    trial = core$trial,
    specification = core$specification,
    bootstrap_B = as.integer(core$bootstrap_B),
    ci_lower = as.numeric(core$ci_lower),
    point_estimate = as.numeric(core$point_estimate),
    ci_upper = as.numeric(core$ci_upper),
    estimand_group = "CBPS-ATE",
    stringsAsFactors = FALSE
  )

  discharge_table <- read_locked(file.path(
    project_root, "manuscript", "archive_redundant_main_tables_2026-08-27",
    "Former_MT-04_discharge_handling_sensitivity.csv"
  ))
  discharge <- discharge_table[
    discharge_table$Analysis == "Sensitivity" &
      discharge_table$`Handling of live discharge` == "Censoring",
    , drop = FALSE
  ]
  stopifnot(nrow(discharge) == 2L)
  discharge_trial <- c("Anti-MRSA trial" = "anti_mrsa", "Anti-PSA trial" = "anti_psa")
  discharge <- discharge[match(names(discharge_trial), discharge$Trial), , drop = FALSE]
  d <- rbind(
    d,
    data.frame(
      trial = unname(discharge_trial),
      specification = "discharge_as_censoring",
      bootstrap_B = c(1000L, 999L),
      ci_lower = as.numeric(discharge$`RD 95% CI lower`),
      point_estimate = as.numeric(discharge$`7-day AKI risk difference`),
      ci_upper = as.numeric(discharge$`RD 95% CI upper`),
      estimand_group = "ATO",
      stringsAsFactors = FALSE
    )
  )
  # Read the corresponding CBPS-ATE and ATO base analyses from the same locked
  # controlled-reporting source. This adds visual anchors without re-estimation.
  base <- sensitivity_source[
    sensitivity_source$estimand_field == "rd_aki" &
      sensitivity_source$specification %in% c("main_cbps_ate", "overlap_weight_ato_brlogit"),
    , drop = FALSE
  ]
  stopifnot(
    nrow(base) == 4L,
    setequal(base$trial, c("anti_mrsa", "anti_psa")),
    setequal(base$specification, c("main_cbps_ate", "overlap_weight_ato_brlogit"))
  )
  d <- rbind(
    data.frame(
      trial = base$trial,
      specification = base$specification,
      bootstrap_B = as.integer(base$bootstrap_B),
      ci_lower = as.numeric(base$ci_lower),
      point_estimate = as.numeric(base$point_estimate),
      ci_upper = as.numeric(base$ci_upper),
      estimand_group = ifelse(base$specification == "main_cbps_ate", "CBPS-ATE", "ATO"),
      stringsAsFactors = FALSE
    ),
    d
  )
  stopifnot(nrow(d) == 22L)

  specification_order <- c(
    "main_cbps_ate", "aki_cross_t0_uo_boundary", "weight_truncation_1_99", "weight_truncation_5_95",
    "creatinine_only_phenotype", "NEW_MS03", "ABX_TIME_proxy",
    "MICRO_STATE_proxy", "complete_case", "overlap_weight_ato_brlogit", "discharge_as_censoring"
  )
  specification_labels <- c(
    main_cbps_ate = "Complementary CBPS-ATE basic analysis",
    aki_cross_t0_uo_boundary = "Urine-output window crossing time zero",
    weight_truncation_1_99 = "Weights truncated at 1st/99th percentiles",
    weight_truncation_5_95 = "Weights truncated at 5th/95th percentiles",
    creatinine_only_phenotype = "Creatinine-only AKI phenotype",
    NEW_MS03 = "Alternative DAG adjustment set (with race)",
    ABX_TIME_proxy = "Adjustment with antibiotic-timing proxy",
    MICRO_STATE_proxy = "Microbiology-status proxy perturbation",
    complete_case = "Complete-case analysis",
    overlap_weight_ato_brlogit = "Primary ATO competing-risk analysis",
    discharge_as_censoring = "Live discharge treated as censoring"
  )
  row_y <- setNames(c(10.70, 9.70, 8.70, 7.70, 6.70, 5.70, 4.70, 3.70, 2.70, 1.25, 0.25), specification_order)
  d$y <- unname(row_y[d$specification])
  d$specification_label <- unname(specification_labels[d$specification])
  if (anyNA(d$y) || anyNA(d$specification_label)) stop("SF08 has an unmapped row.")
  d$trial_label <- factor(d$trial, levels = names(trial_labels), labels = unname(trial_labels))
  d$estimand_group <- factor(d$estimand_group, levels = c("CBPS-ATE", "ATO"))
  d$row_type <- factor(
    ifelse(d$specification %in% c("main_cbps_ate", "overlap_weight_ato_brlogit"), "Corresponding base analysis", "Sensitivity variant"),
    levels = c("Corresponding base analysis", "Sensitivity variant")
  )
  d$estimate_pp <- 100 * d$point_estimate
  d$lower_pp <- 100 * d$ci_lower
  d$upper_pp <- 100 * d$ci_upper
  d$estimate_ci <- sprintf("%+.1f (%+.1f to %+.1f)", d$estimate_pp, d$lower_pp, d$upper_pp)

  backgrounds <- expand.grid(
    trial_label = levels(d$trial_label),
    group = c("CBPS-ATE", "ATO"),
    stringsAsFactors = FALSE
  )
  backgrounds$ymin <- ifelse(backgrounds$group == "CBPS-ATE", 2.15, -0.20)
  backgrounds$ymax <- ifelse(backgrounds$group == "CBPS-ATE", 11.15, 2.15)
  backgrounds$fill <- ifelse(backgrounds$group == "CBPS-ATE", "#F4F7FA", "#FFF1E6")

  header_data <- expand.grid(
    trial_label = levels(d$trial_label),
    stringsAsFactors = FALSE
  )
  header_data$x <- 43.5
  header_data$y <- 11.48
  header_data$label <- "RD (95% CI), pp"
  group_headers <- rbind(
    data.frame(trial_label = levels(d$trial_label), x = -31.5, y = 11.48,
               label = "CBPS-ATE analyses | full eligible population", stringsAsFactors = FALSE),
    data.frame(trial_label = levels(d$trial_label), x = -31.5, y = 1.92,
               label = "ATO analyses | overlap population", stringsAsFactors = FALSE)
  )

  y_breaks <- unname(row_y)
  y_labels <- unname(specification_labels[specification_order])
  sf08 <- ggplot(d, aes(y = y)) +
    geom_rect(
      data = backgrounds,
      aes(xmin = -33, xmax = 56, ymin = ymin, ymax = ymax, fill = fill),
      inherit.aes = FALSE, colour = NA
    ) +
    scale_fill_identity() +
    geom_hline(yintercept = 2.15, colour = "#B8B8B8", linewidth = 0.45) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "#777777", linewidth = 0.55) +
    geom_vline(xintercept = 31.5, colour = "#D0D0D0", linewidth = 0.40) +
    geom_segment(aes(x = lower_pp, xend = upper_pp, yend = y, colour = estimand_group), linewidth = 0.78) +
    geom_segment(aes(x = lower_pp, xend = lower_pp, y = y - 0.08, yend = y + 0.08, colour = estimand_group), linewidth = 0.65) +
    geom_segment(aes(x = upper_pp, xend = upper_pp, y = y - 0.08, yend = y + 0.08, colour = estimand_group), linewidth = 0.65) +
    geom_point(aes(x = estimate_pp, colour = estimand_group, shape = row_type, size = row_type)) +
    geom_text(aes(x = 43.5, label = estimate_ci), family = "serif", size = 3.15, colour = "#252525") +
    geom_text(
      data = group_headers, aes(x = x, y = y, label = label),
      inherit.aes = FALSE, hjust = 0, family = "serif", fontface = "bold", size = 3.05,
      colour = "#3A3A3A"
    ) +
    geom_text(
      data = header_data, aes(x = x, y = y, label = label),
      inherit.aes = FALSE, family = "serif", fontface = "bold", size = 3.20
    ) +
    facet_wrap(~ trial_label, nrow = 1) +
    scale_colour_manual(values = c("CBPS-ATE" = "#3C5488", "ATO" = "#D55E00"), guide = "none") +
    scale_shape_manual(values = c("Corresponding base analysis" = 16, "Sensitivity variant" = 18), guide = "none") +
    scale_size_manual(values = c("Corresponding base analysis" = 3.65, "Sensitivity variant" = 3.15), guide = "none") +
    scale_y_continuous(
      breaks = y_breaks, labels = y_labels, limits = c(-0.20, 11.78),
      expand = expansion(mult = c(0, 0))
    ) +
    scale_x_continuous(
      breaks = seq(-30, 30, by = 10), limits = c(-33, 56),
      labels = function(x) sprintf("%+d", x), expand = expansion(mult = c(0, 0))
    ) +
    labs(
      x = "7-day AKI risk difference, percentage points (de-escalation minus continuation)",
      y = NULL,
      caption = "Within each estimand, the corresponding base analysis is shown before its sensitivity variants. CBPS-ATE uses the full eligible population; ATO uses the overlap population. Points and bars are locked estimates and 95% bootstrap percentile intervals."
    ) +
    theme_classic(base_size = 10.4, base_family = "serif") +
    theme(
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = 11.8, hjust = 0.5, margin = margin(b = 6)),
      axis.text.y = element_text(size = 8.7, colour = "#303030", lineheight = 0.98),
      axis.text.x = element_text(size = 9),
      axis.title.x = element_text(size = 10.3, margin = margin(t = 8)),
      panel.grid.major.x = element_line(colour = "#E4E4E4", linewidth = 0.30),
      panel.grid.major.y = element_blank(),
      panel.spacing.x = grid::unit(0.45, "inches"),
      plot.caption = element_text(size = 8.5, hjust = 0, colour = "grey30", margin = margin(t = 8)),
      plot.margin = margin(10, 12, 8, 10)
    )

  write.csv(d, file.path(audit_dir, "SF-08_sensitivity_refined_display_data.csv"), row.names = FALSE)
  save_figure(sf08, "SF08_sensitivity_forest_refined", 15.2, 9.6)
  message("SF08 complete: 18 locked CBPS-ATE rows and 4 locked ATO rows retained, including both corresponding base analyses.")
}

# -----------------------------------------------------------------------------
# SF03：下三角相关矩阵，去除对称重复单元格，保留全部锁定系数。
# -----------------------------------------------------------------------------
if (should_run("SF03")) {
  corr <- read_locked(file.path(locked_dir, "SF-03_continuous_correlation_data.csv"))
  corr$spearman_rho <- as.numeric(corr$spearman_rho)
  variable_order <- c(
    "age_years", "liver_disease_severity", "baseline_creatinine_mg_dl",
    "target_antibiotic_administration_count_72h", "recent_nonrenal_sofa",
    "nonrenal_sofa_change", "creatinine_change_mg_dl", "urine_output_rate_6h_ml_kg_h"
  )
  short_labels <- c(
    age_years = "Age, years",
    liver_disease_severity = "Liver disease severity",
    baseline_creatinine_mg_dl = "Baseline creatinine",
    target_antibiotic_administration_count_72h = "Target-antibiotic administrations, 72 h",
    recent_nonrenal_sofa = "Recent nonrenal SOFA",
    nonrenal_sofa_change = "Change in nonrenal SOFA",
    creatinine_change_mg_dl = "Pre-index creatinine change",
    urine_output_rate_6h_ml_kg_h = "6-h urine-output rate"
  )
  stopifnot(
    setequal(unique(corr$variable_1), variable_order),
    setequal(unique(corr$variable_2), variable_order),
    nrow(corr) == 128L
  )
  corr$index_1 <- match(corr$variable_1, variable_order)
  corr$index_2 <- match(corr$variable_2, variable_order)
  display <- corr[corr$index_2 >= corr$index_1, , drop = FALSE]
  stopifnot(nrow(display) == 72L)
  display$label_1 <- factor(unname(short_labels[display$variable_1]), levels = unname(short_labels[variable_order]))
  display$label_2 <- factor(unname(short_labels[display$variable_2]), levels = rev(unname(short_labels[variable_order])))
  display$trial_label <- factor(display$trial, levels = names(trial_labels), labels = unname(trial_labels))
  display$text_colour <- ifelse(abs(display$spearman_rho) >= 0.55, "white", "#252525")
  display$rho_label <- ifelse(
    abs(display$spearman_rho) < 0.005,
    "0.00",
    sprintf("%.2f", display$spearman_rho)
  )

  sf03 <- ggplot(display, aes(label_1, label_2, fill = spearman_rho)) +
    geom_tile(colour = "white", linewidth = 0.75) +
    geom_text(
      aes(label = rho_label, colour = text_colour),
      family = "serif", size = 3.05
    ) +
    facet_wrap(~ trial_label, nrow = 1) +
    scale_colour_identity() +
    scale_fill_gradient2(
      low = "#2166AC", mid = "white", high = "#B2182B",
      midpoint = 0, limits = c(-1, 1), breaks = c(-1, -0.5, 0, 0.5, 1),
      name = "Spearman rho"
    ) +
    coord_fixed() +
    labs(
      title = "Spearman correlations among continuous variables",
      subtitle = "Lower triangle and diagonal; coefficients are calculated from the locked analysis data",
      x = NULL, y = NULL,
      caption = "SOFA, Sequential Organ Failure Assessment."
    ) +
    theme_minimal(base_size = 10.1, base_family = "serif") +
    theme(
      plot.title = element_text(face = "bold", size = 13.2, margin = margin(b = 3)),
      plot.subtitle = element_text(size = 9.4, colour = "grey30", margin = margin(b = 7)),
      strip.text = element_text(face = "bold", size = 11.2, hjust = 0.5, margin = margin(b = 5)),
      axis.text.x = element_text(angle = 42, hjust = 1, size = 8.4, colour = "#303030"),
      axis.text.y = element_text(size = 8.7, colour = "#303030"),
      panel.grid = element_blank(),
      panel.spacing.x = grid::unit(0.36, "inches"),
      legend.position = "right",
      legend.title = element_text(face = "bold"),
      plot.caption = element_text(size = 8.5, hjust = 0, colour = "grey30", margin = margin(t = 7)),
      plot.margin = margin(10, 12, 8, 10)
    )

  write.csv(display, file.path(audit_dir, "SF-03_correlation_refined_display_data.csv"), row.names = FALSE)
  save_figure(sf03, "SF03_continuous_correlation_heatmap_refined", 12.4, 7.0)
  message("SF03 complete: 72 lower-triangle/diagonal cells retained from 128 locked matrix rows.")
}
