# ============================================================
# Figure 4：初始抗假单胞菌药物探索性分层的表格化森林图
# 点估计与95%区间均为既有锁定结果；本脚本不调用上方可选展示层bootstrap函数。
# 代码来源：项目既有 plot_main_findings.R 与锁定药物分层结果。
# 引用/commit：RD-4.6-ANALYSIS-OUTPUT-PLAN-V1.0；当前项目非 Git 仓库，无 commit ID。
# 访问日期：2026-08-28。
# 适配说明：仅压缩画布、统一 Times New Roman、增加精确CI和样本量列；不修改估计值。
# 验证命令：D:/software/R4.4.2/R-4.4.2/bin/x64/Rscript.exe manuscript/plot_main_findings.R
# ============================================================
suppressMessages({library(survival); library(ggplot2)})
project_root <- normalizePath(Sys.getenv("TTE_PROJECT_ROOT", unset = "."), winslash = "/", mustWork = TRUE)
runtime_root <- normalizePath(Sys.getenv("TTE_RUNTIME_ROOT", unset = "."), winslash = "/", mustWork = TRUE)
data_dir <- file.path(project_root, "data")
wt_dir   <- file.path(runtime_root, "09_输出结果", "03_权重")
out_dir  <- file.path(project_root, "manuscript", "plot_new")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(20260809)

# 出院删失敏感性：出院=删失、死亡=竞争，7天AKI RD
calc_disch_sens <- function(trial, csv, n, B=200) {
  df <- read.csv(file.path(data_dir, csv), check.names=FALSE)
  stopifnot(nrow(df)==n)
  # 用第一份插补的权重（展示层近似；CI 为患者级重抽样）
  wt <- readRDS(file.path(wt_dir, sprintf("weights_%s_imp01.rds", trial)))
  keys <- as.character(wt$key_map$icu_stay_id)
  aki <- setNames(as.numeric(df$aki_7d), as.character(df$icu_stay_id))[keys]
  day <- setNames(as.numeric(df$aki_followup_days), as.character(df$icu_stay_id))[keys]
  ec  <- setNames(as.character(df$aki_7d_event_code), as.character(df$icu_stay_id))[keys]
  st  <- setNames(as.integer(df$deescalation), as.character(df$icu_stay_id))[keys]
  w   <- wt$ato_brlogit$weight
  ev <- ifelse(aki==1, 1, ifelse(ec=="2", 2, 0))
  rd_func <- function(idx_s) {
    r <- sapply(0:1, function(s) {
      i <- idx_s[st[idx_s] == s]
      if (length(i) < 5) return(NA_real_)
      f <- survfit(Surv(day[i], factor(ev[i], levels=0:2)) ~ 1, weights = w[i])
      s2 <- summary(f, times = 7, extend = TRUE)$pstate
      s2[1, "1"]
    })
    r[2] - r[1]
  }
  full <- rd_func(seq_len(n))
  b <- replicate(B, rd_func(sample.int(n, n, replace = TRUE)))
  c(rd = full, lo = quantile(b, .025, na.rm = TRUE), hi = quantile(b, .975, na.rm = TRUE))
}
# 主分析 ATO RD（ST-08）
main_rd <- c(anti_mrsa=-0.010, anti_psa=-0.070)
# Display-only figure refresh: do not rerun the optional bootstrap calculations above.
# The Figure 4 estimates and intervals below are locked values carried forward unchanged.

# Figure 4：药物分层（探索性，抗PSA 含/不含哌拉西林-他唑巴坦）— JAMA 风格森林图
# 正确值：完整 ATO + 竞争风险 AJ（20插补均值 + bootstrap CI）
d5 <- data.frame(
  Stratum = c("Without piperacillin-tazobactam", "With piperacillin-tazobactam"),
  RD = c(-0.071, -0.068),        # 降阶梯 - 继续
  lo = c(-0.195, -0.261),
  hi = c(0.066, 0.138),
  Cont_n = c(212, 103), Deesc_n = c(35, 22),
  stringsAsFactors=FALSE)
d5$y <- c(2, 1)
d5$estimate_ci <- sprintf("%+.1f (%+.1f to %+.1f)", 100*d5$RD, 100*d5$lo, 100*d5$hi)
d5$n_text <- sprintf("%d/%d", d5$Cont_n, d5$Deesc_n)

# 锁定值断言：防止展示层调整意外更改数字。
stopifnot(
  identical(round(d5$RD, 3), c(-0.071, -0.068)),
  identical(round(d5$lo, 3), c(-0.195, -0.261)),
  identical(round(d5$hi, 3), c(0.066, 0.138)),
  identical(d5$Cont_n, c(212, 103)),
  identical(d5$Deesc_n, c(35, 22))
)

header <- data.frame(
  x = c(22.0, 37.0), y = c(2.54, 2.54),
  label = c("RD (95% CI), pp", "n (continuation/de-escalation)"),
  stringsAsFactors = FALSE
)

p5 <- ggplot(d5, aes(y=y)) +
  geom_vline(xintercept=0, linetype="dashed", colour="#777777", linewidth=0.55) +
  geom_vline(xintercept=c(16.5, 31.5), colour="#E1E1E1", linewidth=0.35) +
  geom_segment(aes(x=lo*100, xend=hi*100, yend=y), linewidth=0.9, colour="#3C5488") +
  geom_segment(aes(x=lo*100, xend=lo*100, y=y-0.09, yend=y+0.09), linewidth=0.75, colour="#3C5488") +
  geom_segment(aes(x=hi*100, xend=hi*100, y=y-0.09, yend=y+0.09), linewidth=0.75, colour="#3C5488") +
  geom_point(aes(x=RD*100), shape=21, size=4.0, stroke=0.65, fill="#3C5488", colour="white") +
  geom_text(aes(x=22.0, label=estimate_ci), hjust=0.5, size=3.7,
            family="serif", colour="#252525") +
  geom_text(aes(x=37.0, label=n_text), hjust=0.5, size=3.7,
            family="serif", colour="#444444") +
  geom_text(data=header, aes(x=x, y=y, label=label), inherit.aes=FALSE,
            fontface="bold", size=3.25, family="serif", colour="#252525") +
  scale_y_continuous(
    breaks=c(2,1), labels=c("Without piperacillin-tazobactam", "With piperacillin-tazobactam"),
    limits=c(0.62,2.68), expand=expansion(mult=c(0,0))
  ) +
  scale_x_continuous(breaks=seq(-30,10,10), limits=c(-33,44), expand=expansion(mult=c(0,0))) +
  labs(x="7-day AKI risk difference, percentage points (de-escalation minus continuation)", y=NULL,
       title="Post hoc exploratory stratification by initial antipseudomonal agent",
       subtitle="Overlap-weighted 7-day AKI risk differences with 95% bootstrap percentile intervals") +
  theme_classic(base_size=11.2, base_family="serif") +
  theme(
    panel.grid.major.x=element_line(colour="#E8E8E8", linewidth=0.32),
    panel.grid.major.y=element_blank(),
    plot.title=element_text(face="bold", size=13.2, margin=margin(b=3)),
    plot.subtitle=element_text(size=9.4, colour="grey30", margin=margin(b=8)),
    axis.text.y=element_text(size=10.8, colour="#252525"),
    axis.text.x=element_text(size=9.5),
    axis.title.x=element_text(size=10.5, margin=margin(t=8)),
    plot.margin=margin(10, 18, 10, 10)
  )

audit_dir <- file.path(project_root, "04_主文图表", "09_tables")
dir.create(audit_dir, recursive=TRUE, showWarnings=FALSE)
write.csv(d5[c("Stratum","RD","lo","hi","Cont_n","Deesc_n","estimate_ci","n_text")],
          file.path(audit_dir,"MF-04_drug_stratification_display_data.csv"), row.names=FALSE)

ggsave(file.path(out_dir,"Figure4_drug_stratification_RD.png"), p5, width=12.0, height=3.8, dpi=320, bg="white")
ggsave(file.path(out_dir,"Figure4_drug_stratification_RD.pdf"), p5, width=12.0, height=3.8, device=grDevices::cairo_pdf)
ggsave(file.path(out_dir,"Figure4_drug_stratification_RD.tiff"), p5, width=12.0, height=3.8, dpi=600, bg="white", device=grDevices::tiff, compression="lzw")
cat("Figure 4（药物分层）已生成\n")
