# =============================================================================
# 脚本用途：以锁定的主分析 ATO（overlap weighting）权重对象重绘补充图 5–7。
# 方法来源：锁定的主分析权重对象（bias-reduced logistic PS + overlap weighting）；
#   平衡诊断遵循 cobalt::bal.tab 的标准化均数差计算。
# 引用 / commit：RD-4.6-ANALYSIS-OUTPUT-PLAN-V1.0；项目锁定权重对象。
# 访问日期：2026-08-20。
# 适配说明：仅重建展示层诊断图；不重新估计倾向评分、插补、结局模型或效应量。
# 验证命令：Rscript --vanilla plot_primary_ato_diagnostics_2026-08-20.R
# =============================================================================

options(stringsAsFactors = FALSE, warn = 1)
set.seed(42)

suppressPackageStartupMessages({
  library(cobalt)
  library(ggplot2)
})

# 以当前脚本位置定位项目，避免依赖工作目录。
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("Cannot determine the current script path.")
script_path <- normalizePath(sub("^--file=", "", script_arg[[1L]]), winslash = "/", mustWork = TRUE)
project_root <- normalizePath(Sys.getenv("TTE_PROJECT_ROOT", unset = "."), winslash = "/", mustWork = TRUE)
data_dir <- file.path(project_root, "data")
runtime_root <- normalizePath(Sys.getenv("TTE_RUNTIME_ROOT", unset = "."), winslash = "/", mustWork = TRUE)
weight_dir <- file.path(runtime_root, "09_输出结果", "03_权重")

if (!dir.exists(weight_dir)) stop("Locked weight directory was not found: ", weight_dir)

save_figure <- function(plot_object, stem, width, height) {
  # 输出提交所需三种格式；图形仅由锁定的权重对象派生。
  ggsave(file.path(data_dir, paste0(stem, ".png")), plot_object,
         width = width, height = height, units = "in", dpi = 320, bg = "white")
  ggsave(file.path(data_dir, paste0(stem, ".tiff")), plot_object,
         width = width, height = height, units = "in", dpi = 600,
         device = "tiff", compression = "lzw", bg = "white")
  ggsave(file.path(data_dir, paste0(stem, ".pdf")), plot_object,
         width = width, height = height, units = "in",
         device = grDevices::cairo_pdf, bg = "white")
}

trial_labels <- c(anti_mrsa = "A. Anti-MRSA trial", anti_psa = "B. Antipseudomonal trial")
strategy_labels <- c("0" = "Continuation", "1" = "De-escalation")
colour_strategy <- c("Continuation" = "#0072B2", "De-escalation" = "#D55E00")
colour_balance <- c("Unweighted" = "#D55E00", "ATO weighted" = "#0072B2")

theme_submission <- theme_minimal(base_size = 10.5, base_family = "Times New Roman") +
  theme(
    plot.title = element_text(size = 13, face = "bold"),
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 9),
    strip.text = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 9),
    legend.title = element_text(size = 9),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank()
  )

# 使用既有受控 Love plot 文件仅作为英文显示标签和排序字典；不读取或复用其 CBPS 数值。
display_dictionary <- read.csv(
  file.path(runtime_root, "09_输出结果", "09_tables", "controlled_reporting", "SF-06_love_plot_data_controlled.csv"),
  check.names = FALSE
)
display_dictionary <- unique(display_dictionary[, c("term", "display_order", "display_label")])
display_dictionary <- display_dictionary[order(display_dictionary$display_order), ]
if (!any(display_dictionary$display_label == "Age, y")) {
  stop("Expected display-layer age label was not found.")
}
display_dictionary$display_label[display_dictionary$display_label == "Age, y"] <- "Age, years"

balance_rows <- list()
ps_rows <- list()
weight_rows <- list()

for (trial in names(trial_labels)) {
  for (imputation in seq_len(20L)) {
    # 读取锁定的插补-权重对象；不对对象、模型或数据作任何写入。
    weight_file <- file.path(weight_dir, sprintf("weights_%s_imp%02d.rds", trial, imputation))
    if (!file.exists(weight_file)) stop("Missing locked weight object: ", weight_file)
    wt <- readRDS(weight_file)
    dat <- wt$completed_data
    ato_weight <- as.numeric(wt$ato_brlogit$weight)
    ato_ps <- as.numeric(wt$ato_brlogit$propensity_score)
    if (length(ato_weight) != nrow(dat) || length(ato_ps) != nrow(dat)) {
      stop("ATO weight/PS length mismatch in ", basename(weight_file))
    }
    if (any(!is.finite(ato_weight)) || any(!is.finite(ato_ps))) {
      stop("Non-finite ATO weight or propensity score in ", basename(weight_file))
    }

    # 按主分析倾向评分公式计算加权前后 SMD；结果用于 Love plot 的展示汇总。
    balance <- cobalt::bal.tab(
      x = wt$ps_formula,
      data = dat,
      weights = ato_weight,
      un = TRUE,
      s.d.denom = "pooled"
    )$Balance
    balance_rows[[length(balance_rows) + 1L]] <- rbind(
      data.frame(trial = trial, imputation = imputation, term = rownames(balance),
                 stage = "Unweighted", smd = abs(balance$Diff.Un), row.names = NULL),
      data.frame(trial = trial, imputation = imputation, term = rownames(balance),
                 stage = "ATO weighted", smd = abs(balance$Diff.Adj), row.names = NULL)
    )

    # 诊断用患者级倾向评分跨 20 次插补平均，避免将同一患者重复计入密度曲线。
    if (is.null(wt$key_map$icu_stay_id)) stop("ICU-stay identifier unavailable in ", basename(weight_file))
    ps_rows[[length(ps_rows) + 1L]] <- data.frame(
      trial = trial,
      imputation = imputation,
      patient = as.character(wt$key_map$icu_stay_id),
      strategy = as.character(dat$deescalation),
      ps = ato_ps,
      stringsAsFactors = FALSE
    )

    # 权重密度使用全部 20 个锁定插补对象，反映主分析权重分布。
    weight_rows[[length(weight_rows) + 1L]] <- data.frame(
      trial = trial,
      imputation = imputation,
      strategy = as.character(dat$deescalation),
      weight = ato_weight,
      stringsAsFactors = FALSE
    )
  }
}

# -----------------------------------------------------------------------------
# Supplementary Figure 5: ATO Love plot (unweighted versus overlap-weighted).
# -----------------------------------------------------------------------------
balance_all <- do.call(rbind, balance_rows)
balance_all <- merge(balance_all, display_dictionary, by = "term", all.x = TRUE, sort = FALSE)
if (anyNA(balance_all$display_label) || anyNA(balance_all$display_order)) {
  stop("An ATO balance term is missing from the display dictionary.")
}

love_key <- interaction(
  balance_all$trial, balance_all$stage, balance_all$term,
  balance_all$display_order, balance_all$display_label, drop = TRUE, lex.order = TRUE
)
love_split <- split(balance_all$smd, love_key)
love_values <- t(vapply(
  love_split,
  function(x) c(minimum = min(x), median = median(x), maximum = max(x)),
  numeric(3)
))
love_first <- balance_all[match(names(love_split), love_key),
                          c("trial", "stage", "term", "display_order", "display_label")]
love_summary <- cbind(love_first, as.data.frame(love_values, row.names = NULL))
love_summary$stage <- factor(love_summary$stage, levels = c("Unweighted", "ATO weighted"))
love_summary$trial_label <- factor(love_summary$trial, levels = names(trial_labels), labels = unname(trial_labels))
label_levels <- unique(love_summary$display_label[order(love_summary$display_order)])
love_summary$display_label <- factor(love_summary$display_label, levels = rev(label_levels))

p_love <- ggplot(love_summary, aes(x = median, y = display_label, colour = stage, shape = stage)) +
  geom_vline(xintercept = 0.10, linetype = "dashed", colour = "#666666", linewidth = 0.45) +
  geom_segment(aes(x = minimum, xend = maximum, yend = display_label), linewidth = 0.65, alpha = 0.78) +
  geom_point(size = 2.45) +
  facet_wrap(~trial_label, nrow = 1, scales = "free_y") +
  scale_colour_manual(values = colour_balance) +
  scale_shape_manual(values = c("Unweighted" = 1, "ATO weighted" = 16)) +
  scale_x_continuous(limits = c(0, max(0.36, max(love_summary$maximum) * 1.05)),
                     breaks = seq(0, 0.4, by = 0.1), expand = expansion(mult = c(0, 0.02))) +
  labs(
    x = "Absolute standardised mean difference", y = NULL, colour = NULL, shape = NULL,
    caption = "Points are medians and horizontal lines are ranges across 20 imputed datasets. Dashed line: SMD = 0.10."
  ) +
  theme_submission +
  theme(
    legend.position = "top",
    axis.text.y = element_text(size = 8.2, lineheight = 1.02),
    panel.spacing = grid::unit(0.7, "lines"),
    plot.caption = element_text(size = 8.2, colour = "grey35"),
    plot.margin = margin(8, 10, 8, 10)
  )

save_figure(p_love, "Supplementary Figure 5_ATO Love Plot", width = 13.5, height = 7.5)

# -----------------------------------------------------------------------------
# Supplementary Figure 6: primary ATO propensity-score overlap.
# -----------------------------------------------------------------------------
ps_all <- do.call(rbind, ps_rows)
ps_summary <- stats::aggregate(ps ~ trial + patient + strategy, data = ps_all, FUN = mean)
if (anyDuplicated(ps_summary[c("trial", "patient")])) stop("Duplicate patient rows after PS pooling.")
ps_summary$strategy_label <- factor(ps_summary$strategy, levels = names(strategy_labels), labels = unname(strategy_labels))
ps_summary$trial_label <- factor(ps_summary$trial, levels = names(trial_labels), labels = unname(trial_labels))

p_ps <- ggplot(ps_summary, aes(x = ps, colour = strategy_label, fill = strategy_label, linetype = strategy_label)) +
  geom_density(alpha = 0.20, linewidth = 1.0, adjust = 1.05) +
  facet_wrap(~trial_label, nrow = 1) +
  scale_colour_manual(values = colour_strategy) +
  scale_fill_manual(values = colour_strategy) +
  scale_linetype_manual(values = c("Continuation" = "solid", "De-escalation" = "22")) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2), expand = c(0, 0)) +
  labs(x = "Propensity score for de-escalation", y = "Density", colour = NULL, fill = NULL, linetype = NULL) +
  theme_submission +
  theme(
    legend.position = "bottom",
    panel.grid.major.y = element_line(colour = "#E8E8E8", linewidth = 0.25),
    panel.spacing.x = grid::unit(0.28, "in")
  )

save_figure(p_ps, "Supplementary Figure 6_ATO Propensity Score Overlap", width = 11.5, height = 4.8)

# -----------------------------------------------------------------------------
# Supplementary Figure 7: ATO weights only; no CBPS-ATE diagnostic panel.
# -----------------------------------------------------------------------------
weight_all <- do.call(rbind, weight_rows)
if (any(weight_all$weight <= 0)) stop("Nonpositive overlap weight prevents log10 display.")
weight_all$strategy_label <- factor(weight_all$strategy, levels = names(strategy_labels), labels = unname(strategy_labels))
weight_all$trial_label <- factor(
  weight_all$trial,
  levels = names(trial_labels),
  labels = unname(trial_labels)
)
weight_all$log10_weight <- log10(weight_all$weight)

p_weight <- ggplot(weight_all, aes(x = log10_weight, colour = strategy_label, fill = strategy_label, linetype = strategy_label)) +
  geom_density(alpha = 0.20, linewidth = 1.0, adjust = 1.05) +
  facet_wrap(~trial_label, nrow = 1) +
  scale_colour_manual(values = colour_strategy) +
  scale_fill_manual(values = colour_strategy) +
  scale_linetype_manual(values = c("Continuation" = "solid", "De-escalation" = "22")) +
  labs(x = expression(log[10] * " overlap weight"), y = "Density", colour = NULL, fill = NULL, linetype = NULL) +
  theme_submission +
  theme(
    legend.position = "bottom",
    panel.grid.major.y = element_line(colour = "#E8E8E8", linewidth = 0.25),
    panel.spacing.x = grid::unit(0.24, "in"),
    strip.text = element_text(size = 10.5, face = "bold"),
    plot.margin = margin(12, 18, 12, 18)
  )

save_figure(p_weight, "Supplementary Figure 7_ATO Weight Distribution", width = 6.7, height = 3.15)

message("Created: Supplementary Figures 5–7 for the primary ATO analysis only.")
