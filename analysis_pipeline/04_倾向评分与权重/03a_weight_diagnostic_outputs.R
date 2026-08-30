# ==============================================================================
# Analysis: 倾向评分与权重诊断图表
# Method: cobalt标准化差异、PS重叠、权重分布和Kish ESS
# Purpose: 生成MT-01/SF-06/SF-07/SF-08/ST-07所需诊断产物
# Source type: GitHub package documentation + local_knowledge_base
# Source name or URL: https://github.com/ngreifer/cobalt ; https://github.com/ngreifer/WeightIt
# Citation or commit/version: cobalt 4.6.3; WeightIt commit 595d2812c2f4d56ff7c14aaf6fcde76c59aa8c8f
# Access date: 2026-07-29
# Adaptation notes: 20份插补分别诊断；Love图使用各协变量跨插补最差绝对SMD。
# Verification command: Rscript 04_倾向评分与权重/03a_weight_diagnostic_outputs.R
# ==============================================================================

source(
  file.path(Sys.getenv("TTE_RUNTIME_ROOT", unset = "."), "02_配置与变量字典", "analysis_config.R"),
  encoding = "UTF-8"
)
setwd(package_root)
require_packages(c("cobalt", "ggplot2", "flextable", "officer"))

weight_dir <- file.path(output_root, "03_权重")
figure_dir <- file.path(output_root, "08_figures")
table_dir <- file.path(output_root, "09_tables")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

balance_rows <- list()
weight_rows <- list()
ps_rows <- list()
diagnostic_rows <- list()

for (trial in trial_registry$trial) {
  for (imputation_index in seq_len(mice_m)) {
    object_path <- file.path(
      weight_dir,
      paste0(
        "weights_",
        trial,
        "_imp",
        sprintf("%02d", imputation_index),
        ".rds"
      )
    )
    if (!file.exists(object_path)) {
      stop("缺少权重对象：", object_path)
    }
    saved <- readRDS(object_path)
    stopifnot(
      identical(saved$trial, trial),
      saved$imputation == imputation_index
    )
    treatment <- as.integer(as.character(saved$completed_data$deescalation))

    # 主要Love图只使用锁定CBPS-ATE权重。
    balance <- cobalt::bal.tab(
      ps_formula_main,
      data = saved$completed_data,
      weights = saved$cbps_ate$weight,
      estimand = "ATE",
      s.d.denom = "pooled",
      binary = "std",
      continuous = "std",
      un = TRUE
    )$Balance
    balance <- balance[balance$Type != "Distance", , drop = FALSE]
    balance_rows[[length(balance_rows) + 1L]] <- data.frame(
      trial = trial,
      imputation = imputation_index,
      term = rownames(balance),
      absolute_smd_unweighted = abs(balance$Diff.Un),
      absolute_smd_weighted = abs(balance$Diff.Adj),
      stringsAsFactors = FALSE
    )

    for (model_name in c("cbps_ate", "ato_brlogit")) {
      fit <- if (model_name == "cbps_ate") {
        saved$cbps_ate
      } else {
        saved$ato_brlogit
      }
      group_ps_min <- tapply(
        fit$propensity_score,
        treatment,
        min
      )
      group_ps_max <- tapply(
        fit$propensity_score,
        treatment,
        max
      )
      overlap_lower <- max(group_ps_min)
      overlap_upper <- min(group_ps_max)
      continuation_n <- sum(treatment == 0L)
      deescalation_n <- sum(treatment == 1L)
      weight_rows[[length(weight_rows) + 1L]] <- data.frame(
        trial = trial,
        imputation = imputation_index,
        model = model_name,
        row_id = seq_along(fit$weight),
        strategy = treatment,
        propensity_score = fit$propensity_score,
        weight = fit$weight,
        stringsAsFactors = FALSE
      )
      diagnostic_rows[[length(diagnostic_rows) + 1L]] <- data.frame(
        trial = trial,
        imputation = imputation_index,
        model = model_name,
        ps_min = min(fit$propensity_score),
        ps_max = max(fit$propensity_score),
        ps_below_0_01_n = sum(fit$propensity_score < 0.01),
        ps_above_0_99_n = sum(fit$propensity_score > 0.99),
        deescalation_ps_below_0_01_n = sum(
          treatment == 1L & fit$propensity_score < 0.01
        ),
        continuation_ps_above_0_99_n = sum(
          treatment == 0L & fit$propensity_score > 0.99
        ),
        ps_outside_0_01_0_99_fraction = mean(
          fit$propensity_score < 0.01 |
            fit$propensity_score > 0.99
        ),
        empirical_overlap_lower = overlap_lower,
        empirical_overlap_upper = overlap_upper,
        empirical_overlap_width = overlap_upper - overlap_lower,
        empirical_overlap_absent = overlap_upper <= overlap_lower,
        weight_min = min(fit$weight),
        weight_p01 = stats::quantile(fit$weight, 0.01, names = FALSE),
        weight_p05 = stats::quantile(fit$weight, 0.05, names = FALSE),
        weight_p50 = stats::quantile(fit$weight, 0.50, names = FALSE),
        weight_p95 = stats::quantile(fit$weight, 0.95, names = FALSE),
        weight_p99 = stats::quantile(fit$weight, 0.99, names = FALSE),
        weight_max = max(fit$weight),
        weight_above_10_n = sum(fit$weight > 10),
        weight_above_10_fraction = mean(fit$weight > 10),
        continuation_n = continuation_n,
        deescalation_n = deescalation_n,
        ess_continuation = fit$ess_continuation,
        ess_deescalation = fit$ess_deescalation,
        ess_fraction_continuation =
          fit$ess_continuation / continuation_n,
        ess_fraction_deescalation =
          fit$ess_deescalation / deescalation_n,
        max_absolute_smd = fit$max_absolute_smd,
        stringsAsFactors = FALSE
      )
    }
  }
}

balance_long <- do.call(rbind, balance_rows)
weight_long <- do.call(rbind, weight_rows)
diagnostic_long <- do.call(rbind, diagnostic_rows)

atomic_write_csv(
  balance_long,
  file.path(weight_dir, "ST-07_balance_by_imputation.csv")
)
atomic_write_csv(
  diagnostic_long,
  file.path(weight_dir, "ST-07_weight_diagnostics_by_imputation.csv")
)

# 各操作变量在20份插补中的最差绝对SMD，避免用单一插补掩盖残余不平衡。
love_summary <- stats::aggregate(
  cbind(absolute_smd_unweighted, absolute_smd_weighted) ~ trial + term,
  data = balance_long,
  FUN = max
)
love_plot_data <- rbind(
  data.frame(
    trial = love_summary$trial,
    term = love_summary$term,
    stage = "Unweighted",
    absolute_smd = love_summary$absolute_smd_unweighted
  ),
  data.frame(
    trial = love_summary$trial,
    term = love_summary$term,
    stage = "Weighted",
    absolute_smd = love_summary$absolute_smd_weighted
  )
)
love_plot_data$term <- factor(
  love_plot_data$term,
  levels = rev(unique(love_summary$term))
)
love_plot <- ggplot2::ggplot(
  love_plot_data,
  ggplot2::aes(x = absolute_smd, y = term, colour = stage, shape = stage)
) +
  ggplot2::geom_vline(xintercept = 0.10, linetype = 2, colour = "#B42318") +
  ggplot2::geom_point(size = 1.8) +
  ggplot2::facet_wrap(~trial, ncol = 1, scales = "free_y") +
  ggplot2::labs(
    title = "Covariate balance before and after CBPS-ATE weighting",
    subtitle = "Worst absolute standardized mean difference across 20 imputations",
    x = "Absolute standardized mean difference",
    y = NULL
  ) +
  ggplot2::theme_minimal(base_size = 9) +
  ggplot2::theme(panel.grid.major.y = ggplot2::element_blank())
ggplot2::ggsave(
  file.path(figure_dir, "SF-06_love_plot_cbps_ate.png"),
  love_plot,
  width = 9,
  height = 12,
  dpi = 320
)
ggplot2::ggsave(
  file.path(figure_dir, "SF-06_love_plot_cbps_ate.pdf"),
  love_plot,
  width = 9,
  height = 12
)

# 对同一患者跨20份插补取平均PS，用于展示总体重叠而不重复计算样本量。
cbps_weights <- weight_long[weight_long$model == "cbps_ate", , drop = FALSE]
ps_summary <- stats::aggregate(
  propensity_score ~ trial + row_id + strategy,
  data = cbps_weights,
  FUN = mean
)
ps_summary$strategy_label <- factor(
  ps_summary$strategy,
  levels = c(0, 1),
  labels = c("Continuation", "De-escalation")
)
ps_plot <- ggplot2::ggplot(
  ps_summary,
  ggplot2::aes(
    x = propensity_score,
    colour = strategy_label,
    fill = strategy_label
  )
) +
  ggplot2::geom_density(alpha = 0.20, linewidth = 0.7) +
  ggplot2::facet_wrap(~trial, ncol = 1) +
  ggplot2::labs(
    title = "Propensity-score overlap",
    subtitle = "Patient-level mean propensity score across 20 imputations",
    x = "Propensity score for de-escalation",
    y = "Density",
    colour = "Strategy",
    fill = "Strategy"
  ) +
  ggplot2::theme_minimal(base_size = 10)
ggplot2::ggsave(
  file.path(figure_dir, "SF-07_propensity_score_overlap.png"),
  ps_plot,
  width = 8,
  height = 7,
  dpi = 320
)
ggplot2::ggsave(
  file.path(figure_dir, "SF-07_propensity_score_overlap.pdf"),
  ps_plot,
  width = 8,
  height = 7
)

weight_long$model_label <- factor(
  weight_long$model,
  levels = c("cbps_ate", "ato_brlogit"),
  labels = c("CBPS-ATE", "ATO (br.logit)")
)
weight_plot <- ggplot2::ggplot(
  weight_long,
  ggplot2::aes(x = weight, fill = factor(strategy))
) +
  ggplot2::geom_histogram(
    bins = 50,
    alpha = 0.55,
    position = "identity"
  ) +
  ggplot2::facet_grid(trial ~ model_label, scales = "free") +
  ggplot2::labs(
    title = "Distribution of analysis weights",
    subtitle = "All 20 imputations; no silent truncation",
    x = "Weight",
    y = "Count",
    fill = "Strategy"
  ) +
  ggplot2::theme_minimal(base_size = 10)
ggplot2::ggsave(
  file.path(figure_dir, "SF-08_weight_distribution.png"),
  weight_plot,
  width = 10,
  height = 7,
  dpi = 320
)
ggplot2::ggsave(
  file.path(figure_dir, "SF-08_weight_distribution.pdf"),
  weight_plot,
  width = 10,
  height = 7
)

# ST-07按试验和模型汇总20份插补中的范围与最差诊断。
st07_rows <- list()
for (trial in trial_registry$trial) {
  for (model_name in c("cbps_ate", "ato_brlogit")) {
    current <- diagnostic_long[
      diagnostic_long$trial == trial &
        diagnostic_long$model == model_name,
      ,
      drop = FALSE
    ]
    st07_rows[[length(st07_rows) + 1L]] <- data.frame(
      trial = trial,
      model = model_name,
      imputations = nrow(current),
      ps_minimum = min(current$ps_min),
      ps_maximum = max(current$ps_max),
      maximum_ps_below_0_01_n = max(current$ps_below_0_01_n),
      maximum_ps_above_0_99_n = max(current$ps_above_0_99_n),
      maximum_deescalation_ps_below_0_01_n = max(
        current$deescalation_ps_below_0_01_n
      ),
      maximum_continuation_ps_above_0_99_n = max(
        current$continuation_ps_above_0_99_n
      ),
      maximum_ps_outside_0_01_0_99_fraction = max(
        current$ps_outside_0_01_0_99_fraction
      ),
      minimum_empirical_overlap_width = min(
        current$empirical_overlap_width
      ),
      empirical_overlap_absent_in_any_imputation = any(
        current$empirical_overlap_absent
      ),
      minimum_weight = min(current$weight_min),
      minimum_weight_p01 = min(current$weight_p01),
      minimum_weight_p05 = min(current$weight_p05),
      minimum_weight_median = min(current$weight_p50),
      maximum_weight_median = max(current$weight_p50),
      maximum_weight_p95 = max(current$weight_p95),
      maximum_weight_p99 = max(current$weight_p99),
      maximum_weight = max(current$weight_max),
      maximum_weight_above_10_n = max(current$weight_above_10_n),
      maximum_weight_above_10_fraction = max(
        current$weight_above_10_fraction
      ),
      minimum_ess_continuation = min(current$ess_continuation),
      minimum_ess_deescalation = min(current$ess_deescalation),
      minimum_ess_fraction_continuation = min(
        current$ess_fraction_continuation
      ),
      minimum_ess_fraction_deescalation = min(
        current$ess_fraction_deescalation
      ),
      maximum_absolute_smd = max(current$max_absolute_smd),
      stringsAsFactors = FALSE
    )
  }
}
st07 <- do.call(rbind, st07_rows)

# 对原样本主CBPS-ATE实施Analysis Lock规定的平衡与ESS硬门槛。
main_gate <- st07[st07$model == "cbps_ate", , drop = FALSE]
if (
  any(main_gate$maximum_absolute_smd >= 0.10) ||
    any(main_gate$minimum_ess_fraction_continuation <= 0.25) ||
    any(main_gate$minimum_ess_fraction_deescalation <= 0.25)
) {
  stop("原样本CBPS-ATE未通过预设SMD或策略组ESS门槛。")
}

atomic_write_csv(st07, file.path(table_dir, "ST-07_weight_diagnostics.csv"))
st07_flex <- flextable::flextable(st07)
st07_flex <- flextable::theme_booktabs(st07_flex)
st07_flex <- flextable::autofit(st07_flex)
flextable::save_as_docx(
  "Supplementary Table 7. Propensity-score and weight diagnostics" =
    st07_flex,
  path = file.path(table_dir, "ST-07_weight_diagnostics.docx")
)

# MT-01已在Analysis Lock阶段完成；这里只做编号映射和版式输出，不重新计算。
locked_mt01_path <- file.path(
  config_dir,
  "completed_MT-01_weighted_baseline.csv"
)
if (!file.exists(locked_mt01_path)) {
  stop("缺少已锁定MT-01输入：", locked_mt01_path)
}
mt01 <- utils::read.csv(
  locked_mt01_path,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
atomic_write_csv(
  mt01,
  file.path(table_dir, "MT-01_weighted_baseline_characteristics.csv")
)
mt01_flex <- flextable::flextable(mt01)
mt01_flex <- flextable::theme_booktabs(mt01_flex)
mt01_flex <- flextable::autofit(mt01_flex)
flextable::save_as_docx(
  "Main Table 1. Weighted baseline characteristics" = mt01_flex,
  path = file.path(
    table_dir,
    "MT-01_weighted_baseline_characteristics.docx"
  )
)

message("SF-06、SF-07、SF-08和ST-07已生成。")
