# ============================================================
# Figure 2 新方案：加权 7 天 AKI 累积发生曲线（主分析 = Overlap weighting, ATO）
# 数据：final.data（aki_7d, aki_followup_days, aki_7d_event_code, deescalation）
#       + 锁定 ATO (ato_brlogit) 权重
# 方法：加权 Aalen-Johansen（survfit 多状态），AKI=1、死亡=2、出院=3 为竞争事件
#       20 插补点估计取均值；不绘制逐时间点置信带
# Day-7 标注使用正文/锁定 ATO 结果表值（风险 + 95%CI + RD/RR）
# 注意：仅使用锁定数据与权重，不修改任何数据/模型
# 代码来源：项目锁定权重对象、survival::survfit 与既有 Figure 2 绘图脚本。
# 引用/commit：RD-4.6-ANALYSIS-OUTPUT-PLAN-V1.0；当前项目非 Git 仓库，无 commit ID。
# 访问日期：2026-08-28。
# 适配说明：仅将纵轴收紧至 0%-40%、增加注释框并统一 Nature 风格；不运行 bootstrap。
# 验证命令：D:/software/R4.4.2/R-4.4.2/bin/x64/Rscript.exe manuscript/plot_figure2_aki_curve_ato.R
# ============================================================

suppressMessages({
  library(survival)
  library(ggplot2)
})

project_root <- normalizePath(Sys.getenv("TTE_PROJECT_ROOT", unset = "."), winslash = "/", mustWork = TRUE)
runtime_root <- normalizePath(Sys.getenv("TTE_RUNTIME_ROOT", unset = "."), winslash = "/", mustWork = TRUE)
data_dir  <- file.path(project_root, "data")
wt_dir    <- file.path(runtime_root, "09_输出结果", "03_权重")
out_dir   <- file.path(project_root, "manuscript", "plot_new")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(20260809)

# 锁定 ATO 结果（正文/ST-07）：Day-7 风险与 CI、RD、RR
trials <- list(
  anti_mrsa = list(csv = "final.data_anti-MRSA.csv", label = "A. Anti-MRSA trial",
                   cont_risk = 0.308, cont_lo = 0.229, cont_hi = 0.388,
                   deesc_risk = 0.298, deesc_lo = 0.207, deesc_hi = 0.403,
                   rd = -0.010, rd_lo = -0.129, rd_hi = 0.114,
                   rr = 0.97, rr_lo = 0.63, rr_hi = 1.42),
  anti_psa  = list(csv = "final.data_anti-PSA.csv", label = "B. Antipseudomonal trial",
                   cont_risk = 0.240, cont_lo = 0.180, cont_hi = 0.311,
                   deesc_risk = 0.170, deesc_lo = 0.077, deesc_hi = 0.260,
                   rd = -0.070, rd_lo = -0.181, rd_hi = 0.026,
                   rr = 0.71, rr_lo = 0.30, rr_hi = 1.12)
)

# ---------- 计算某试验：20 插补 ATO 加权 AJ 曲线点估计 ----------
calc_aki_curve_point <- function(trial_key, csv_name, n_expected) {
  df <- read.csv(file.path(data_dir, csv_name), check.names = FALSE)
  stopifnot(nrow(df) == n_expected)
  aki_map  <- setNames(as.numeric(df$aki_7d), as.character(df$icu_stay_id))
  day_map  <- setNames(as.numeric(df$aki_followup_days), as.character(df$icu_stay_id))
  ec_map   <- setNames(as.integer(df$aki_7d_event_code), as.character(df$icu_stay_id))
  deesc_map<- setNames(as.integer(df$deescalation), as.character(df$icu_stay_id))
  # 事件：0=随访期满, 1=AKI, 2=死亡, 3=出院
  ev_func <- function(ec) ifelse(ec == 1, 1, ifelse(ec == 2, 2, ifelse(ec == 3, 3, 0)))
  row_list <- list()
  for (imp in 1:20) {
    wt <- readRDS(file.path(wt_dir, sprintf("weights_%s_imp%02d.rds", trial_key, imp)))
    km <- wt$key_map
    keys <- as.character(km$icu_stay_id)
    aki  <- unname(aki_map[keys]); day <- unname(day_map[keys])
    ec   <- unname(ec_map[keys]); strat <- unname(deesc_map[keys])
    stopifnot(!anyNA(aki) && !anyNA(day) && !anyNA(ec) && !anyNA(strat))
    ev <- ev_func(ec)
    w  <- wt$ato_brlogit$weight
    for (s in 0:1) {
      idx <- which(strat == s)
      if (length(idx) < 5) next
      fit <- survfit(Surv(day[idx], factor(ev[idx], levels = 0:3)) ~ 1,
                     weights = w[idx])
      sm <- summary(fit, times = seq(0, 7, by = 0.25), extend = TRUE)
      # pstate 列：状态顺序 0,1,2,3；AKI 为状态 1
      aki_risk <- sm$pstate[, match("1", colnames(sm$pstate))]
      row_list[[length(row_list) + 1L]] <- data.frame(
        trial = trial_key, imputation = imp, strategy = s,
        time_days = seq(0, 7, by = 0.25), aki_risk = aki_risk)
    }
  }
  out <- do.call(rbind, row_list)
  agg <- aggregate(aki_risk ~ trial + strategy + time_days, out, mean)
  agg
}

# ---------- 生成两试验曲线 ----------
pt <- rbind(
  calc_aki_curve_point("anti_mrsa", "final.data_anti-MRSA.csv", 268L),
  calc_aki_curve_point("anti_psa",  "final.data_anti-PSA.csv",  372L)
)
curves <- pt
curves$strategy_label <- factor(curves$strategy, levels = c(0L, 1L),
                                labels = c("Continuation", "De-escalation"))
curves$trial_label <- factor(curves$trial, levels = c("anti_mrsa", "anti_psa"),
                        labels = c("A. Anti-MRSA trial", "B. Antipseudomonal trial"))

# ---------- Day-7 标注（锁定 ATO 值） ----------
ann <- do.call(rbind, lapply(names(trials), function(tk) {
  t <- trials[[tk]]
  data.frame(trial = tk,
             annotation = sprintf(
               "Day-7 risks (95%% CI)\nContinuation: %.1f%% (%.1f%% to %.1f%%)\nDe-escalation: %.1f%% (%.1f%% to %.1f%%)\nRD: %+.1f pp (%.1f to %.1f)\nRR: %.2f (%.2f to %.2f)",
               100*t$cont_risk, 100*t$cont_lo, 100*t$cont_hi,
               100*t$deesc_risk, 100*t$deesc_lo, 100*t$deesc_hi,
               100*t$rd, 100*t$rd_lo, 100*t$rd_hi,
               t$rr, t$rr_lo, t$rr_hi),
             stringsAsFactors = FALSE)
}))
ann$trial_label <- factor(ann$trial, levels = c("anti_mrsa", "anti_psa"),
                        labels = c("A. Anti-MRSA trial", "B. Antipseudomonal trial"))
ann$x <- 0.20
ann$y <- 0.392

# Save the exact display-layer curve coordinates and locked annotation text.
# These files provide a direct audit trail from plotted points to the figure.
audit_dir <- file.path(project_root, "04_主文图表", "09_tables")
dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(curves,
          file.path(audit_dir, "MF-02_weighted_aki_curve_display_data.csv"),
          row.names = FALSE)
write.csv(ann[, c("trial", "annotation")],
          file.path(audit_dir, "MF-02_weighted_aki_curve_locked_annotations.csv"),
          row.names = FALSE)

# ---------- 绘图 ----------
strategy_colors <- c("Continuation" = "#0072B2", "De-escalation" = "#D55E00")
strategy_linetype <- c("Continuation" = "solid", "De-escalation" = "22")

p <- ggplot(curves, aes(time_days, aki_risk, colour = strategy_label,
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
  scale_x_continuous(breaks = seq(0, 7, by = 1), limits = c(0, 7)) +
  scale_y_continuous(labels = function(x) paste0(round(100 * x), "%"),
                     breaks = seq(0, 0.4, by = 0.1),
                     limits = c(0, 0.4), expand = expansion(mult = c(0, 0))) +
  labs(title = "Weighted cumulative incidence of 7-day in-hospital AKI",
       subtitle = "ATO point estimates averaged across 20 imputations; Day-7 intervals from locked nested bootstrap",
       x = "Days since time zero", y = "Cumulative AKI incidence",
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

ggsave(file.path(out_dir, "Figure2_7day_aki_curve_ATO.png"), p,
       width = 8.5, height = 7.4, dpi = 320, bg = "white")
ggsave(file.path(out_dir, "Figure2_7day_aki_curve_ATO.pdf"), p,
       width = 8.5, height = 7.4, device = grDevices::cairo_pdf)
ggsave(file.path(out_dir, "Figure2_7day_aki_curve_ATO.tiff"), p,
       width = 8.5, height = 7.4, dpi = 600, device = "tiff", compression = "lzw", bg = "white")

# ---------- 输出 Day-7 点估计对照（自检） ----------
cat("Day-7 ATO 点估计（20 插补均值）:\n")
print(subset(curves, time_days == 7)[, c("trial","strategy","aki_risk")])
cat("\n锁定正文 Day-7 值:\n")
for (tk in names(trials)) {
  t <- trials[[tk]]
  cat(sprintf("%s: 继续 %.3f, 降阶梯 %.3f\n", tk, t$cont_risk, t$deesc_risk))
}
cat("\n生成完成:", file.path(out_dir, "Figure2_7day_aki_curve_ATO.png"), "\n")
