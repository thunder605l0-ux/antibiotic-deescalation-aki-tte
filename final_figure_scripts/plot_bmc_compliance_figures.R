# ============================================================
# 用途：为 BMC 投稿生成 Figure 1 及 Supplementary Figures 2、3、4、7 的展示性修订图。
# 数据来源：RD-4.6-ANALYSIS-OUTPUT-PLAN-V1.0 锁定 CSV 与经审计的权重 RDS。
# 引用/版本：publication display dictionary v1.0；原始绘图逻辑 20b_generate_publication_figures_en.R。
# 访问日期：2026-08-10。
# 适配说明：仅调整文字、单位、面板布局和版式；不拟合模型、不变更数据或结果数值。
# 验证命令：D:/software/R4.4.2/R-4.4.2/bin/x64/Rscript.exe manuscript/plot_bmc_compliance_figures.R
# ============================================================

options(stringsAsFactors = FALSE, warn = 1)
suppressMessages({
  library(ggplot2)
  library(gridExtra)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("Cannot determine the current script path.")
script_path <- normalizePath(sub("^--file=", "", script_arg[[1L]]), winslash = "/", mustWork = TRUE)
public_code_root <- dirname(dirname(script_path))
project_dir <- normalizePath(Sys.getenv("TTE_PROJECT_ROOT", unset = "."), winslash = "/", mustWork = TRUE)
runtime_dir <- normalizePath(Sys.getenv("TTE_RUNTIME_ROOT", unset = "."), winslash = "/", mustWork = TRUE)
out_dir <- file.path(project_dir, "manuscript", "plot_new")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

trial_labels <- c(anti_mrsa = "Anti-MRSA trial", anti_psa = "Antipseudomonal trial")
strategy_colors <- c("Continuation" = "#0072B2", "De-escalation" = "#D55E00")
strategy_linetypes <- c("Continuation" = "solid", "De-escalation" = "22")

save_three <- function(plot_object, stem, width, height) {
  paths <- c(
    png = file.path(out_dir, paste0(stem, ".png")),
    pdf = file.path(out_dir, paste0(stem, ".pdf")),
    tiff = file.path(out_dir, paste0(stem, ".tiff"))
  )
  ggsave(paths[["png"]], plot_object, width = width, height = height, dpi = 320, bg = "white")
  ggsave(paths[["pdf"]], plot_object, width = width, height = height,
         device = grDevices::cairo_pdf, bg = "white")
  ggsave(paths[["tiff"]], plot_object, width = width, height = height, dpi = 300,
         device = "tiff", compression = "lzw", bg = "white")
  if (any(file.info(paths)$size <= 0)) stop("A candidate figure is empty: ", stem)
  paths
}

base_theme <- theme_minimal(base_size = 9.5, base_family = "sans") +
  theme(
    axis.title = element_text(size = 10.5),
    axis.text = element_text(size = 8.5),
    strip.text = element_text(size = 10, face = "bold"),
    panel.grid.minor = element_blank(),
    plot.margin = margin(8, 10, 8, 10)
  )

# Figure 1: compact six-step cohort flow diagram.
# 独立脚本从原13步锁定CSV派生六步展示数据，完成算术断言并同步正式输出。
source(
  file.path(
    public_code_root, "analysis_pipeline", "08_未完成图表",
    "21_rebuild_mf01_compact_nature_2026-08-28.R"
  ),
  encoding = "UTF-8"
)

# Supplementary Figure 2: missing-data heat map, with explicit percent unit.
missing_long <- read.csv(file.path(runtime_dir, "09_输出结果", "09_tables", "controlled_reporting",
                                   "ST-05_missing_data_long_audit.csv"),
                         check.names = FALSE, stringsAsFactors = FALSE)
missing_long <- missing_long[missing_long$stratum != "overall", , drop = FALSE]
missing_long$trial_label <- factor(missing_long$trial, levels = names(trial_labels), labels = unname(trial_labels))
missing_long$stratum_label <- factor(missing_long$stratum,
                                     levels = c("continued", "deescalated"),
                                     labels = c("Continuation", "De-escalation"))
variable_order <- unique(missing_long$variable[order(missing_long$display_order)])
missing_long$variable <- factor(missing_long$variable, levels = rev(variable_order))
sf2 <- ggplot(missing_long, aes(stratum_label, variable, fill = missing_pct)) +
  geom_tile(colour = "white", linewidth = 0.25) +
  geom_text(aes(label = ifelse(missing_pct == 0, "", sprintf("%.1f", missing_pct))), family = "sans", size = 2.1) +
  facet_wrap(~trial_label, nrow = 1) +
  scale_fill_gradient(low = "white", high = "#0072B2", name = "Missing (%)") +
  labs(x = NULL, y = NULL) + base_theme +
  theme(axis.text.x = element_text(angle = 30, hjust = 1), panel.grid = element_blank())
save_three(sf2, "SF02_missing_data_bmc", 12, 12)

# Supplementary Figure 3: correlation heat map, with the statistic named explicitly.
corr <- read.csv(file.path(runtime_dir, "09_输出结果", "09_tables", "controlled_reporting",
                           "SF-03_continuous_correlation_data.csv"),
                 check.names = FALSE, stringsAsFactors = FALSE)
dictionary <- read.csv(file.path(project_dir, "03_数据", "03_变量字典与配置",
                                 "publication_display_labels_v1.0.csv"),
                       check.names = FALSE, stringsAsFactors = FALSE)
display_by_field <- setNames(dictionary$display_label, dictionary$clean_field)
corr$label1 <- unname(display_by_field[corr$variable_1])
corr$label2 <- unname(display_by_field[corr$variable_2])
if (anyNA(corr$label1) || anyNA(corr$label2)) stop("Supplementary Figure 3 contains an unmapped variable.")
corr$trial_label <- factor(corr$trial, levels = names(trial_labels), labels = unname(trial_labels))
cor_order <- unique(corr$label1[corr$trial == "anti_mrsa"])
corr$label1 <- factor(corr$label1, levels = cor_order)
corr$label2 <- factor(corr$label2, levels = rev(cor_order))
sf3 <- ggplot(corr, aes(label1, label2, fill = spearman_rho)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.2f", spearman_rho)), family = "sans", size = 2.35) +
  facet_wrap(~trial_label, nrow = 1) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
                       limits = c(-1, 1), name = "Spearman correlation (rho)") +
  labs(x = NULL, y = NULL) + base_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank())
save_three(sf3, "SF03_correlation_heatmap_bmc", 12, 7)

# Supplementary Figure 4: trial-specific factor order prevents anti-PSA variables becoming NA.
chain <- rbind(
  read.csv(file.path(project_dir, "05_补充图", "09_tables", "controlled_reporting",
                     "SF-04_SF-05_mice_chain_statistics_anti_mrsa.csv"), check.names = FALSE),
  read.csv(file.path(project_dir, "05_补充图", "09_tables", "controlled_reporting",
                     "SF-04_SF-05_mice_chain_statistics_anti_psa.csv"), check.names = FALSE)
)
chain$variable_label <- unname(display_by_field[chain$variable])
if (anyNA(chain$variable_label)) stop("Supplementary Figure 4 contains an unmapped variable.")
build_chain_block <- function(trial_id, block_label, show_x_title) {
  d <- chain[chain$trial == trial_id, , drop = FALSE]
  variable_order <- unique(d$variable_label)
  d$variable_label <- factor(d$variable_label, levels = variable_order)
  d$chain_group <- interaction(d$variable_label, d$chain, drop = TRUE)
  ggplot(d, aes(iteration, chain_mean, group = chain_group)) +
    geom_line(colour = "#4D4D4D", alpha = 0.38, linewidth = 0.32, na.rm = TRUE) +
    facet_wrap(~variable_label, ncol = 4, scales = "free_y") +
    labs(title = block_label, x = if (show_x_title) "Iteration" else NULL, y = "Chain mean") +
    base_theme +
    theme(plot.title = element_text(size = 11, face = "bold", hjust = 0),
          axis.text.x = element_text(size = 6.7), axis.text.y = element_text(size = 6.7),
          strip.text = element_text(size = 7.2, face = "bold"),
          panel.grid.major = element_line(colour = "#ECECEC", linewidth = 0.2),
          panel.spacing = grid::unit(0.22, "inches"), plot.margin = margin(6, 8, 6, 8))
}
mrsa_block <- build_chain_block("anti_mrsa", "A. Anti-MRSA trial", FALSE)
psa_block <- build_chain_block("anti_psa", "B. Anti-PSA trial", TRUE)
sf4 <- gridExtra::arrangeGrob(mrsa_block, psa_block, ncol = 1)
save_three(sf4, "SF04_mice_chain_means_bmc", 13.5, 10.8)

# Supplementary Figure 7: untruncated weight distributions; no recalculation of analysis weights.
weight_dir <- file.path(runtime_dir, "09_输出结果", "03_权重")
weight_rows <- list()
for (trial in names(trial_labels)) {
  for (imp in seq_len(20L)) {
    path <- file.path(weight_dir, sprintf("weights_%s_imp%02d.rds", trial, imp))
    if (!file.exists(path)) stop("Missing locked weight object: ", path)
    saved <- readRDS(path)
    treatment <- as.integer(as.character(saved$completed_data$deescalation))
    key <- paste(saved$key_map$patient_id, saved$key_map$hospital_admission_id, saved$key_map$icu_stay_id, sep = "-")
    for (model in c("cbps_ate", "ato_brlogit")) {
      fit <- saved[[model]]
      weight_rows[[length(weight_rows) + 1L]] <- data.frame(
        trial = trial, imputation = imp, model = model, patient_key = key, strategy = treatment,
        weight = as.numeric(fit$weight), stringsAsFactors = FALSE
      )
    }
  }
}
weight_long <- do.call(rbind, weight_rows)
if (any(!is.finite(weight_long$weight)) || any(weight_long$weight <= 0)) stop("Invalid locked analysis weight.")
weight_long$trial_label <- factor(weight_long$trial, levels = names(trial_labels), labels = unname(trial_labels))
weight_long$model_label <- factor(weight_long$model, levels = c("cbps_ate", "ato_brlogit"),
                                  labels = c("CBPS-ATE", "Overlap weighting (ATO)"))
weight_long$strategy_label <- factor(ifelse(weight_long$strategy == 1L, "De-escalation", "Continuation"),
                                     levels = names(strategy_colors))
sf7 <- ggplot(weight_long, aes(log10(weight), colour = strategy_label, fill = strategy_label,
                               linetype = strategy_label)) +
  geom_density(alpha = 0.12, linewidth = 0.65) +
  facet_grid(trial_label ~ model_label, scales = "free_y") +
  scale_colour_manual(values = strategy_colors) +
  scale_fill_manual(values = strategy_colors) +
  scale_linetype_manual(values = strategy_linetypes) +
  scale_x_continuous(breaks = log10(c(0.01, 0.03, 0.1, 0.3, 1, 3, 10)),
                     labels = c("0.01", "0.03", "0.10", "0.30", "1", "3", "10")) +
  labs(x = "Weight (log10 scale)", y = "Density", colour = NULL, fill = NULL, linetype = NULL) +
  base_theme +
  theme(legend.position = "bottom", panel.grid.major.y = element_blank())
save_three(sf7, "SF07_weight_distribution_bmc", 11.5, 8)

cat("BMC candidate figures generated in: ", out_dir, "\n", sep = "")
