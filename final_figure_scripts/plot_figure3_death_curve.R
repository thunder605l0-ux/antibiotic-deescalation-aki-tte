# ============================================================
# Figure 3 新方案：加权 30 天累积死亡曲线（主分析 = Overlap weighting, ATO）
# 数据：final.data（death_30d, death_time_from_index_days, deescalation）+ 锁定 ATO 权重
# 方法：加权 Kaplan-Meier（survfit, stype=1），累积死亡 = 1 - 生存；20 插补点估计取均值
# Day-30 标注使用锁定 ATO 结果表值（风险 + 95%CI + RD/RR）
# 注意：本脚本仅使用锁定数据与权重，不修改任何数据
# 代码来源：项目锁定权重对象、survival::survfit 与既有 Figure 3 绘图脚本。
# 引用/commit：RD-4.6-ANALYSIS-OUTPUT-PLAN-V1.0；当前项目非 Git 仓库，无 commit ID。
# 访问日期：2026-08-28。
# 适配说明：仅将纵轴收紧至 0%-30%、增加注释框并统一 Nature 风格；不运行 bootstrap。
# 验证命令：D:/software/R4.4.2/R-4.4.2/bin/x64/Rscript.exe manuscript/plot_figure3_death_curve.R
# ============================================================

suppressMessages({
  library(survival)
  library(ggplot2)
})

# ---------- 输入 ----------
project_root <- normalizePath(Sys.getenv("TTE_PROJECT_ROOT", unset = "."), winslash = "/", mustWork = TRUE)
runtime_root <- normalizePath(Sys.getenv("TTE_RUNTIME_ROOT", unset = "."), winslash = "/", mustWork = TRUE)
data_dir  <- file.path(project_root, "data")
wt_dir    <- file.path(runtime_root, "09_输出结果", "03_权重")
out_dir   <- file.path(project_root, "manuscript", "plot_new")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

trials <- list(
  anti_mrsa = list(csv = "final.data_anti-MRSA.csv", label = "A. Anti-MRSA trial",
                   cont_risk = 0.137, cont_lo = 0.086, cont_hi = 0.195,
                   deesc_risk = 0.116, deesc_lo = 0.049, deesc_hi = 0.183,
                   rd = -0.021, rd_lo = -0.103, rd_hi = 0.060,
                   rr = 0.84, rr_lo = 0.34, rr_hi = 1.62),
  anti_psa  = list(csv = "final.data_anti-PSA.csv", label = "B. Antipseudomonal trial",
                   cont_risk = 0.116, cont_lo = 0.076, cont_hi = 0.170,
                   deesc_risk = 0.181, deesc_lo = 0.083, deesc_hi = 0.272,
                   rd = 0.065, rd_lo = -0.039, rd_hi = 0.161,
                   rr = 1.56, rr_lo = 0.71, rr_hi = 2.82)
)

# ---------- 计算某试验 20 插补加权累积死亡曲线 ----------
calc_death_curve <- function(trial_key, csv_name, n_expected) {
  df <- read.csv(file.path(data_dir, csv_name), check.names = FALSE)
  stopifnot(nrow(df) == n_expected)
  # 以 icu_stay_id 为连接键：final.data 的死亡时间关联到权重对象
  dt_map <- setNames(df$death_time_from_index_days, as.character(df$icu_stay_id))
  d30_map <- setNames(as.numeric(df$death_30d), as.character(df$icu_stay_id))
  deesc_map <- setNames(as.integer(df$deescalation), as.character(df$icu_stay_id))
  row_list <- list()
  for (imp in 1:20) {
    wt <- readRDS(file.path(wt_dir, sprintf("weights_%s_imp%02d.rds", trial_key, imp)))
    cd <- wt$completed_data
    km <- wt$key_map
    stopifnot(nrow(cd) == n_expected && nrow(km) == n_expected)
    keys <- as.character(km$icu_stay_id)
    # 关联死亡信息（按 ID，不假设行序）；unname 去除命名向量属性
    time <- unname(dt_map[keys])
    ev   <- unname(d30_map[keys])
    strat <- unname(deesc_map[keys])
    stopifnot(!anyNA(ev) && !anyNA(strat))
    # 验证 completed_data 与 final.data 的 death_30d 一致（按 ID；death_30d 为 factor，先转字符）
    stopifnot(identical(as.numeric(as.character(cd$death_30d)), ev))
    # 时间：死亡者 = 死亡时间；未死亡者删失于 30 天
    time <- ifelse(ev == 1, time, 30)
    time[is.na(time)] <- 30
    w  <- wt$ato_brlogit$weight
    for (s in 0:1) {
      idx <- which(strat == s)
      fit <- survfit(Surv(time[idx], ev[idx]) ~ 1, weights = w[idx], stype = 1)
      sfit <- summary(fit, times = seq(0, 30, by = 1), extend = TRUE)
      cum_death <- 1 - sfit$surv
      row_list[[length(row_list) + 1L]] <- data.frame(
        trial = trial_key, imputation = imp, strategy = s,
        time_days = seq(0, 30, by = 1), death_risk = cum_death)
    }
  }
  out <- do.call(rbind, row_list)
  agg <- aggregate(death_risk ~ trial + strategy + time_days, out, mean)
  agg
}

# ---------- 生成两试验曲线 ----------
curves <- rbind(
  calc_death_curve("anti_mrsa", "final.data_anti-MRSA.csv", 268L),
  calc_death_curve("anti_psa",  "final.data_anti-PSA.csv",  372L)
)

curves$strategy_label <- factor(curves$strategy,
                                levels = c(0L, 1L),
                                labels = c("Continuation", "De-escalation"))
curves$trial_label <- factor(curves$trial,
                             levels = c("anti_mrsa", "anti_psa"),
                        labels = c("A. Anti-MRSA trial", "B. Antipseudomonal trial"))

# ---------- Day-30 标注 ----------
ann <- do.call(rbind, lapply(names(trials), function(tk) {
  t <- trials[[tk]]
  data.frame(trial = tk,
             annotation = sprintf(
               "Day-30 risks (95%% CI)\nContinuation: %.1f%% (%.1f%% to %.1f%%)\nDe-escalation: %.1f%% (%.1f%% to %.1f%%)\nRD: %+.1f pp (%.1f to %.1f)\nRR: %.2f (%.2f to %.2f)",
               100*t$cont_risk, 100*t$cont_lo, 100*t$cont_hi,
               100*t$deesc_risk, 100*t$deesc_lo, 100*t$deesc_hi,
               100*t$rd, 100*t$rd_lo, 100*t$rd_hi,
               t$rr, t$rr_lo, t$rr_hi),
             stringsAsFactors = FALSE)
}))
ann$trial_label <- factor(ann$trial, levels = c("anti_mrsa", "anti_psa"),
                        labels = c("A. Anti-MRSA trial", "B. Antipseudomonal trial"))
ann$x <- 0.45
ann$y <- 0.292

# Save the exact display-layer curve coordinates and locked annotation text.
# These files provide a direct audit trail from plotted points to the figure.
audit_dir <- file.path(project_root, "04_主文图表", "09_tables")
dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(curves,
          file.path(audit_dir, "MF-03_weighted_mortality_curve_display_data.csv"),
          row.names = FALSE)
write.csv(ann[, c("trial", "annotation")],
          file.path(audit_dir, "MF-03_weighted_mortality_curve_locked_annotations.csv"),
          row.names = FALSE)

# ---------- 绘图 ----------
strategy_colors <- c("Continuation" = "#0072B2", "De-escalation" = "#D55E00")
strategy_linetype <- c("Continuation" = "solid", "De-escalation" = "22")

p <- ggplot(curves, aes(time_days, death_risk, colour = strategy_label,
                        linetype = strategy_label)) +
  geom_step(linewidth = 1.0, direction = "hv") +
  geom_label(data = ann, aes(x = x, y = y, label = annotation),
             inherit.aes = FALSE, hjust = 0, vjust = 1,
             size = 2.85, lineheight = 1.03, colour = "#252525",
             family = "serif", fill = "white",
             linewidth = 0.28, label.r = grid::unit(0.08, "lines"),
             label.padding = grid::unit(0.28, "lines")) +
  facet_wrap(~ trial_label, ncol = 1) +
  scale_colour_manual(values = strategy_colors, drop = FALSE) +
  scale_linetype_manual(values = strategy_linetype, drop = FALSE) +
  scale_x_continuous(breaks = seq(0, 30, by = 5), limits = c(0, 30)) +
  scale_y_continuous(labels = function(x) paste0(round(100 * x), "%"),
                     breaks = seq(0, 0.3, by = 0.1),
                     limits = c(0, 0.3), expand = expansion(mult = c(0, 0))) +
  labs(title = "Weighted cumulative incidence of 30-day all-cause mortality",
       subtitle = "Overlap weighting (ATO); Day-30 values from the locked result tables",
       x = "Days since time zero", y = "Cumulative mortality",
       colour = NULL, linetype = NULL) +
  theme_classic(base_size = 11.2, base_family = "serif") +
  theme(
    plot.title = element_text(face = "bold", size = 13.2, margin = margin(b = 3)),
    plot.subtitle = element_text(size = 9.4, colour = "grey30", margin = margin(b = 7)),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 11.5, hjust = 0, margin = margin(t = 4, b = 4)),
    legend.position = "bottom",
    legend.key.width = grid::unit(0.55, "inches"),
    panel.grid.major.y = element_line(colour = "#E3E3E3", linewidth = 0.32),
    panel.spacing.y = grid::unit(0.28, "inches"),
    plot.margin = margin(10, 12, 8, 10)
  )

ggsave(file.path(out_dir, "Figure3_30day_mortality_curve.png"), p,
       width = 8.5, height = 7.4, dpi = 320, bg = "white")
ggsave(file.path(out_dir, "Figure3_30day_mortality_curve.pdf"), p,
       width = 8.5, height = 7.4, device = grDevices::cairo_pdf)
ggsave(file.path(out_dir, "Figure3_30day_mortality_curve.tiff"), p,
       width = 8.5, height = 7.4, dpi = 600, device = "tiff", compression = "lzw", bg = "white")

# ---------- 输出 Day-30 点估计对照（自检） ----------
cat("Day-30 点估计（20 插补均值）:\n")
print(subset(curves, time_days == 30))
cat("\n锁定 Table 2/3 Day-30 值:\n")
for (tk in names(trials)) {
  t <- trials[[tk]]
  cat(sprintf("%s: 继续 %.3f, 降阶梯 %.3f\n", tk, t$cont_risk, t$deesc_risk))
}
cat("\n生成完成:", file.path(out_dir, "Figure3_30day_mortality_curve.png"), "\n")
