# Analysis: Generate English publication-ready figures from locked plot data and weight objects
# Date: 2026-08-08
# Random seed: 20260808 for deterministic DAG label placement only
# R: 4.4.2
#
# Method: deterministic plotting of locked aggregate CSVs and locked propensity-score/weight RDS objects
# Source: controlled MF-01/MF-02/SF-01 to SF-08 plot data and audited weight objects
# Citation/version: RD-4.6-ANALYSIS-OUTPUT-PLAN-V1.0; publication display dictionary v1.0
# Access date: 2026-08-08
# Adaptation: Times New Roman, consistent strategy encodings, reader-facing English labels, PNG/TIFF/PDF output.
# Verification command: D:/software/R4.4.2/R-4.4.2/bin/x64/Rscript.exe 08_未完成图表/20b_generate_publication_figures_en.R

options(stringsAsFactors = FALSE, warn = 1)

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg[[1L]]), winslash = "/", mustWork = TRUE)
project_root <- dirname(dirname(script_path))
workspace_root <- dirname(project_root)
source(file.path(project_root, "02_配置与变量字典", "publication_style_helpers.R"), encoding = "UTF-8")

for (pkg in c("ggplot2", "ggrepel", "gridExtra", "digest")) {
  if (!requireNamespace(pkg, quietly = TRUE)) stop("Package '", pkg, "' is required.")
}

output_root <- file.path(project_root, "09_输出结果", "11_publication_ready_en")
main_dir <- file.path(output_root, "01_main_figures")
supp_dir <- file.path(output_root, "03_supplementary_figures")
audit_dir <- file.path(output_root, "_audit")
dir.create(main_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(supp_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)

# Keep grid/ggrepel metric calculations off R's implicit Rplots.pdf device.
metric_device_path <- tempfile("publication_font_metrics_", fileext = ".pdf")
grDevices::cairo_pdf(metric_device_path, width = 7, height = 7, family = publication_font_family)
on.exit({
  if (grDevices::dev.cur() > 1L) grDevices::dev.off()
  unlink(metric_device_path, force = TRUE)
}, add = TRUE)

dictionary <- load_publication_dictionary(project_root)
display_by_field <- setNames(dictionary$display_label, dictionary$clean_field)

read_locked_csv <- function(relative_path) {
  path <- file.path(project_root, relative_path)
  if (!file.exists(path)) stop("Missing locked figure source: ", path)
  utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE, fileEncoding = "UTF-8")
}

trial_labels <- c(anti_mrsa = "Anti-MRSA trial", anti_psa = "Anti-PSA trial")
strategy_labels <- c(`0` = "Continuation", `1` = "De-escalation")

with_font_metric_guard <- function(expr) {
  withCallingHandlers(
    expr,
    warning = function(warning_condition) {
      known_null_device_warning <- grepl(
        "font family 'Times New Roman' not found in PostScript font database",
        conditionMessage(warning_condition), fixed = TRUE
      )
      if (known_null_device_warning) invokeRestart("muffleWarning")
    }
  )
}

save_figure <- function(plot_object, id, stem, width, height) {
  target_dir <- if (startsWith(id, "MF")) main_dir else supp_dir
  paths <- file.path(target_dir, paste0(stem, c(".png", ".tiff", ".pdf")))
  # ggrepel queries the null PostScript device while arranging labels on Windows.
  # The saved PDF font is independently verified below with pdftools::pdf_fonts().
  with_font_metric_guard(ggplot2::ggsave(paths[[1L]], plot_object, width = width, height = height, units = "in", dpi = 320, bg = "white"))
  with_font_metric_guard(ggplot2::ggsave(paths[[2L]], plot_object, width = width, height = height, units = "in", dpi = 600, device = "tiff", compression = "lzw", bg = "white"))
  with_font_metric_guard(ggplot2::ggsave(paths[[3L]], plot_object, width = width, height = height, units = "in", device = grDevices::cairo_pdf, bg = "white"))
  paths
}

label_audit_rows <- list()
register_labels <- function(id, label_type, values) {
  values <- unique(as.character(values))
  values <- values[!is.na(values) & nzchar(values)]
  bad <- grepl("race_white|White race|factor\\(|anti_mrsa|anti_psa|cbps_ate|ato_brlogit|NEW_MS03|ABX_TIME_proxy|MICRO_STATE_proxy", values)
  bad <- bad | grepl("[A-Za-z0-9]+_[A-Za-z0-9_]+", values)
  label_audit_rows[[length(label_audit_rows) + 1L]] <<- data.frame(
    artifact_id = id, label_type = label_type, label = values, pass = !bad, stringsAsFactors = FALSE
  )
  if (any(bad)) stop(id, " contains a forbidden reader-facing label: ", values[which(bad)[1L]])
}

generated <- list()

# MF-01: cohort flow diagram.
flow <- read_locked_csv("09_输出结果/09_tables/MF-01_cohort_flow_diagram_data.csv")
flow$trial_label <- factor(flow$trial, levels = names(trial_labels), labels = paste0(c("A. ", "B. "), unname(trial_labels)))
flow$node_text <- gsub("Target-antibiotic", "Target antibiotic", flow$node_text, fixed = TRUE)
flow$node_text <- gsub("target-antibiotic", "target antibiotic", flow$node_text, fixed = TRUE)
flow$node_text <- gsub("T0", "index time", flow$node_text, fixed = TRUE)
flow$node_text <- gsub("target-MDRO", "target MDRO", flow$node_text, fixed = TRUE)
arrow_data <- do.call(rbind, lapply(split(flow, flow$trial), function(d) {
  d <- d[order(d$step_order), , drop = FALSE]
  if (nrow(d) < 2L) return(NULL)
  data.frame(trial_label = d$trial_label[-nrow(d)], y = d$y[-nrow(d)] - 0.22, yend = d$y[-1L] + 0.22)
}))
mf01 <- ggplot2::ggplot(flow, ggplot2::aes(x = 1, y = y)) +
  ggplot2::geom_segment(
    data = arrow_data, ggplot2::aes(x = 1, xend = 1, y = y, yend = yend), inherit.aes = FALSE,
    linewidth = 0.45, colour = "#555555", arrow = grid::arrow(length = grid::unit(0.07, "inches"))
  ) +
  ggplot2::geom_label(
    ggplot2::aes(label = node_text), family = publication_font_family, size = 2.55,
    linewidth = 0.25, label.padding = grid::unit(0.16, "lines"), fill = "white", colour = "black"
  ) +
  ggplot2::facet_wrap(~trial_label, nrow = 1) +
  ggplot2::coord_cartesian(xlim = c(0.45, 1.55), clip = "off") +
  ggplot2::labs(title = "Cohort construction for the two emulated target trials") +
  ggplot2::theme_void(base_family = publication_font_family, base_size = 9) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, size = 11),
    strip.text = ggplot2::element_text(face = "bold", size = 10),
    panel.spacing = grid::unit(0.45, "inches"), plot.margin = ggplot2::margin(8, 18, 8, 18)
  )
register_labels("MF-01", "node", c(levels(flow$trial_label), flow$node_text))
generated[["MF-01"]] <- save_figure(mf01, "MF-01", "MF-01_cohort_flow_diagram", 12, 13)

# MF-02: weighted cumulative incidence of AKI.
curve <- read_locked_csv("09_输出结果/07_图表/controlled_reporting/MF-02_weighted_aki_cumulative_incidence_controlled_data.csv")
annotations <- read_locked_csv("09_输出结果/07_图表/controlled_reporting/MF-02_weighted_aki_cumulative_incidence_controlled_annotations.csv")
curve$strategy_label <- factor(curve$strategy_label, levels = names(publication_strategy_colors))
curve$trial_label <- factor(curve$trial, levels = names(trial_labels), labels = paste0(c("A. ", "B. "), unname(trial_labels)))
annotations$trial_label <- factor(annotations$trial, levels = names(trial_labels), labels = paste0(c("A. ", "B. "), unname(trial_labels)))
mf02 <- ggplot2::ggplot(curve, ggplot2::aes(time_days, aki_risk, colour = strategy_label, fill = strategy_label, linetype = strategy_label)) +
  ggplot2::geom_ribbon(ggplot2::aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.14, colour = NA) +
  ggplot2::geom_step(linewidth = 0.85, direction = "hv") +
  ggplot2::geom_text(
    data = annotations, ggplot2::aes(x = x, y = y, label = annotation), inherit.aes = FALSE,
    family = publication_font_family, size = 3.0, hjust = 0, vjust = 1
  ) +
  ggplot2::facet_wrap(~trial_label, ncol = 1) +
  ggplot2::scale_colour_manual(values = publication_strategy_colors, drop = FALSE) +
  ggplot2::scale_fill_manual(values = publication_strategy_colors, drop = FALSE) +
  ggplot2::scale_linetype_manual(values = publication_strategy_linetypes, drop = FALSE) +
  ggplot2::scale_x_continuous(breaks = 0:7) +
  ggplot2::scale_y_continuous(labels = function(x) paste0(round(100 * x), "%"), limits = c(0, 0.62)) +
  ggplot2::labs(
    title = "Weighted cumulative incidence of 7-day AKI", x = "Days since index time",
    y = "Cumulative incidence", colour = NULL, fill = NULL, linetype = NULL
  ) + publication_theme(10) +
  ggplot2::theme(panel.grid.major.y = ggplot2::element_line(colour = "#E6E6E6", linewidth = 0.25))
register_labels("MF-02", "figure", c(levels(curve$strategy_label), levels(curve$trial_label), "Weighted cumulative incidence of 7-day AKI"))
generated[["MF-02"]] <- save_figure(mf02, "MF-02", "MF-02_weighted_aki_cumulative_incidence", 8.5, 8)

# SF-01: dual DAG with confirmed nodes, edges, and coordinates.
dag_dir <- file.path(workspace_root, "project-meta", "phase-3", "artifacts", "3.4")
dag_nodes <- utils::read.csv(file.path(dag_dir, "dag_batch03_node_dictionary_v1.0_confirmed.csv"), check.names = FALSE)
dag_full_edges <- utils::read.csv(file.path(dag_dir, "dag_batch03_full_causal_edges_v1.0_confirmed.csv"), check.names = FALSE)
dag_measured_edges <- utils::read.csv(file.path(dag_dir, "dag_batch03_measured_operational_edges_v1.0_confirmed.csv"), check.names = FALSE)
dag_coordinates <- utils::read.csv(file.path(dag_dir, "dag_batch03_layout_coordinates_v1.1_confirmed.csv"), check.names = FALSE)
dag_labels <- c(
  AGE = "Age", SEX = "Male sex", RACE_W = "Race", SNF_SRC = "Skilled nursing facility admission source",
  CKD = "Chronic kidney disease", DM = "Diabetes mellitus", CHF = "Congestive heart failure",
  COPD = "Chronic obstructive pulmonary disease", HEM_MAL = "Hematologic malignancy",
  SOLID_MAL = "Solid malignancy", LIVER_DIS = "Liver disease severity", IMMUNOSUPP = "Immunosuppression",
  SCR_BASE = "Baseline serum creatinine", ONSET_SETTING = "Sepsis onset setting", INF_SITE = "Infection site",
  VIRAL_INF = "Acute viral infection status", ABX_TIME = "Time from sepsis onset to target antibiotic initiation",
  TARGET_CLASS = "Initial target antibiotic class", OTHER_BSA = "Other broad-spectrum antibiotic exposure",
  SOURCE_CTRL = "Source-control procedure", NEPHRO_DRUG = "Nephrotoxic drug exposure",
  TARGET_EXPOSURE = "Target antibiotic exposure before index time", CLIN_TRAJ = "Nonrenal SOFA trajectory",
  SCR_DELTA = "Pre-index serum creatinine trajectory", UO_RATE = "Urine-output trajectory",
  MICRO_STATE = "Microbiology status at index time", ALT_ACTIVE = "Alternative active coverage at index time",
  DEESC = "De-escalation strategy", AKI7 = "7-day incident AKI",
  U_INFECTION_CONTROL = "Unmeasured infection control", U_CLINICAL_JUDGMENT = "Unmeasured clinical judgment",
  U_STEWARDSHIP_ACCESS = "Unmeasured stewardship access", U_GOALS_FUNCTION = "Unmeasured goals and function"
)
if (!setequal(names(dag_labels), dag_nodes$node_id)) stop("SF-01 label coverage does not match the confirmed DAG.")
make_dag_panel <- function(edges, included_nodes, title) {
  nodes <- dag_coordinates[dag_coordinates$node_id %in% included_nodes, , drop = FALSE]
  nodes$label <- unname(dag_labels[nodes$node_id])
  nodes$label <- vapply(nodes$label, function(x) paste(strwrap(x, width = 24L), collapse = "\n"), character(1))
  nodes$class <- ifelse(grepl("^U_", nodes$node_id), "Unmeasured", ifelse(nodes$node_id == "DEESC", "Exposure", ifelse(nodes$node_id == "AKI7", "Outcome", "Measured")))
  from <- nodes[match(edges$from, nodes$node_id), c("x", "y")]
  to <- nodes[match(edges$to, nodes$node_id), c("x", "y")]
  if (anyNA(from) || anyNA(to)) stop("SF-01 edge endpoint is missing from the confirmed coordinates.")
  edge_data <- data.frame(x = from$x, y = from$y, xend = to$x, yend = to$y)
  ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = edge_data, ggplot2::aes(x, y, xend = xend, yend = yend), colour = "#888888", linewidth = 0.28,
      arrow = grid::arrow(length = grid::unit(0.035, "inches"), type = "closed")
    ) +
    ggplot2::geom_point(data = nodes, ggplot2::aes(x, y, fill = class, shape = class), size = 3.1, colour = "black", stroke = 0.4) +
    ggrepel::geom_label_repel(
      data = nodes, ggplot2::aes(x, y, label = label), family = publication_font_family,
      size = 2.15, seed = 20260808, box.padding = 0.45, point.padding = 0.25,
      min.segment.length = 0, max.overlaps = Inf, max.iter = 20000, max.time = 3,
      force = 1.25, force_pull = 0.15, segment.colour = "#B3B3B3",
      fill = grDevices::adjustcolor("white", alpha.f = 0.88), linewidth = 0.12
    ) +
    ggplot2::scale_fill_manual(values = c(Measured = "#D9D9D9", Exposure = "#D55E00", Outcome = "#0072B2", Unmeasured = "white")) +
    ggplot2::scale_shape_manual(values = c(Measured = 21, Exposure = 22, Outcome = 23, Unmeasured = 24)) +
    ggplot2::coord_cartesian(xlim = range(dag_coordinates$x) + c(-0.6, 0.6), ylim = range(dag_coordinates$y) + c(-0.7, 0.8), clip = "off") +
    ggplot2::labs(title = title, fill = NULL, shape = NULL) +
    ggplot2::theme_void(base_family = publication_font_family, base_size = 8) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 10), legend.position = "bottom", plot.margin = ggplot2::margin(6, 12, 6, 12))
}
full_panel <- make_dag_panel(dag_full_edges, dag_nodes$node_id[dag_nodes$included_in_full_dag], "A. Full causal DAG")
measured_panel <- make_dag_panel(dag_measured_edges, dag_nodes$node_id[dag_nodes$included_in_measured_dag], "B. Measured operational DAG")
sf01 <- with_font_metric_guard(gridExtra::arrangeGrob(full_panel, measured_panel, ncol = 1))
register_labels("SF-01", "node", c(unname(dag_labels), "Full causal DAG", "Measured operational DAG"))
generated[["SF-01"]] <- save_figure(sf01, "SF-01", "SF-01_dual_DAG", 13, 16)

# SF-02: missing data.
missing_long <- read_locked_csv("09_输出结果/09_tables/controlled_reporting/ST-05_missing_data_long_audit.csv")
mapped <- unname(display_by_field[missing_long$field])
if (anyNA(mapped)) stop("SF-02 contains an unmapped variable.")
missing_long$variable_label <- mapped
missing_long$trial_label <- factor(missing_long$trial, levels = names(trial_labels), labels = unname(trial_labels))
missing_long$stratum_label <- factor(missing_long$stratum, levels = c("overall", "continued", "deescalated"), labels = c("Overall", "Continuation", "De-escalation"))
variable_order <- unique(missing_long$variable_label[order(missing_long$display_order)])
missing_long$variable_label <- factor(missing_long$variable_label, levels = rev(variable_order))
sf02 <- ggplot2::ggplot(missing_long, ggplot2::aes(stratum_label, variable_label, fill = missing_pct)) +
  ggplot2::geom_tile(colour = "white", linewidth = 0.25) +
  ggplot2::geom_text(ggplot2::aes(label = ifelse(missing_pct == 0, "", sprintf("%.1f", missing_pct))), family = publication_font_family, size = 2.1) +
  ggplot2::facet_wrap(~trial_label, nrow = 1) +
  ggplot2::scale_fill_gradient(low = "white", high = "#0072B2", name = "Missing, %") +
  ggplot2::labs(title = "Missing data by trial and treatment strategy", x = NULL, y = NULL) + publication_theme(8.5) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1), panel.grid = ggplot2::element_blank())
register_labels("SF-02", "variable", c(as.character(variable_order), levels(missing_long$trial_label), levels(missing_long$stratum_label)))
generated[["SF-02"]] <- save_figure(sf02, "SF-02", "SF-02_missing_data", 12, 12)

# SF-03: continuous-variable correlation heat map.
corr <- read_locked_csv("09_输出结果/09_tables/controlled_reporting/SF-03_continuous_correlation_data.csv")
label1 <- unname(display_by_field[corr$variable_1]); label2 <- unname(display_by_field[corr$variable_2])
if (anyNA(label1) || anyNA(label2)) stop("SF-03 contains an unmapped variable.")
corr$label1 <- label1; corr$label2 <- label2
corr$trial_label <- factor(corr$trial, levels = names(trial_labels), labels = unname(trial_labels))
cor_order <- unique(label1[corr$trial == "anti_mrsa"])
corr$label1 <- factor(corr$label1, levels = cor_order)
corr$label2 <- factor(corr$label2, levels = rev(cor_order))
sf03 <- ggplot2::ggplot(corr, ggplot2::aes(label1, label2, fill = spearman_rho)) +
  ggplot2::geom_tile(colour = "white", linewidth = 0.4) +
  ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", spearman_rho)), family = publication_font_family, size = 2.35) +
  ggplot2::facet_wrap(~trial_label, nrow = 1) +
  ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, limits = c(-1, 1), name = "Spearman ρ") +
  ggplot2::labs(title = "Spearman correlations among continuous variables", x = NULL, y = NULL) + publication_theme(8.5) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1), panel.grid = ggplot2::element_blank())
register_labels("SF-03", "variable", c(as.character(cor_order), levels(corr$trial_label)))
generated[["SF-03"]] <- save_figure(sf03, "SF-03", "SF-03_continuous_correlation_heatmap", 12, 7)

# SF-04: MICE chain means.
chain <- rbind(
  read_locked_csv("09_输出结果/09_tables/controlled_reporting/SF-04_SF-05_mice_chain_statistics_anti_mrsa.csv"),
  read_locked_csv("09_输出结果/09_tables/controlled_reporting/SF-04_SF-05_mice_chain_statistics_anti_psa.csv")
)
chain$variable_label <- unname(display_by_field[chain$variable])
if (anyNA(chain$variable_label)) stop("SF-04/SF-05 contains an unmapped variable.")
chain$trial_label <- factor(chain$trial, levels = names(trial_labels), labels = unname(trial_labels))
chain$chain_group <- interaction(chain$trial, chain$variable, chain$chain, drop = TRUE)
make_chain_plot <- function(value, title, y_title) {
  ggplot2::ggplot(chain, ggplot2::aes(x = iteration, y = .data[[value]], group = chain_group)) +
    ggplot2::geom_line(colour = "#4D4D4D", alpha = 0.30, linewidth = 0.25, na.rm = TRUE) +
    ggplot2::facet_grid(trial_label ~ variable_label, scales = "free_y") +
    ggplot2::labs(title = title, x = "Iteration", y = y_title) + publication_theme(7.5) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(size = 6), axis.text.y = ggplot2::element_text(size = 6),
      strip.text = ggplot2::element_text(size = 6.5), panel.grid.major = ggplot2::element_line(colour = "#EEEEEE", linewidth = 0.2),
      legend.position = "none"
    )
}
sf04 <- make_chain_plot("chain_mean", "MICE chain means across 20 imputations", "Chain mean")
register_labels("SF-04", "variable", c(unique(chain$variable_label), levels(chain$trial_label)))
generated[["SF-04"]] <- save_figure(sf04, "SF-04", "SF-04_mice_chain_mean", 15, 7.5)

# SF-05: covariate balance.
love <- read_locked_csv("09_输出结果/09_tables/controlled_reporting/SF-06_love_plot_data_controlled.csv")
term_dictionary <- utils::read.csv(file.path(project_root, "02_配置与变量字典", "publication_model_term_labels_v1.0.csv"), check.names = FALSE, stringsAsFactors = FALSE)
term_map <- setNames(term_dictionary$display_label, term_dictionary$term)
love$display_label <- unname(term_map[love$term])
if (anyNA(love$display_label)) stop("SF-05 contains an unmapped balance term.")
love$trial_label <- factor(love$trial, levels = names(trial_labels), labels = unname(trial_labels))
love$stage <- factor(love$stage, levels = c("Unweighted", "Weighted"))
ordered_terms <- term_dictionary$display_label[order(term_dictionary$display_order)]
love$display_label <- factor(love$display_label, levels = rev(ordered_terms))
sf06 <- ggplot2::ggplot(love, ggplot2::aes(median, display_label, colour = stage, shape = stage)) +
  ggplot2::geom_vline(xintercept = 0.10, linetype = "dashed", colour = "#666666", linewidth = 0.45) +
  ggplot2::geom_segment(ggplot2::aes(x = minimum, xend = maximum, yend = display_label), linewidth = 0.55) +
  ggplot2::geom_point(size = 2) +
  ggplot2::facet_wrap(~trial_label, ncol = 1, scales = "free_y") +
  ggplot2::scale_colour_manual(values = c(Unweighted = "#D55E00", Weighted = "#0072B2")) +
  ggplot2::scale_shape_manual(values = c(Unweighted = 1, Weighted = 16)) +
  ggplot2::labs(
    title = "Covariate balance before and after CBPS-ATE weighting",
    subtitle = "Points are medians and horizontal lines are ranges across 20 imputed datasets",
    x = "Absolute standardized mean difference", y = NULL, colour = NULL, shape = NULL,
    caption = "Dashed line indicates the prespecified absolute SMD threshold of 0.10."
  ) + publication_theme(8.5) +
  ggplot2::theme(panel.grid.major.y = ggplot2::element_blank(), legend.position = "top")
register_labels("SF-05", "variable", c(as.character(ordered_terms), levels(love$trial_label), levels(love$stage)))
generated[["SF-05"]] <- save_figure(sf06, "SF-05", "SF-05_love_plot_cbps_ate", 10.5, 12.5)

# SF-07 and SF-08: read locked weight objects; derive plot-only patient summaries.
weight_dir <- file.path(project_root, "09_输出结果", "03_权重")
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
        propensity_score = as.numeric(fit$propensity_score), weight = as.numeric(fit$weight), stringsAsFactors = FALSE
      )
    }
  }
}
weight_long <- do.call(rbind, weight_rows)
if (any(!is.finite(weight_long$propensity_score)) || any(weight_long$propensity_score <= 0 | weight_long$propensity_score >= 1)) stop("Invalid locked propensity score.")
if (any(!is.finite(weight_long$weight)) || any(weight_long$weight <= 0)) stop("Invalid locked weight.")
cbps <- weight_long[weight_long$model == "cbps_ate", , drop = FALSE]
ps_summary <- stats::aggregate(propensity_score ~ trial + patient_key + strategy, cbps, mean)
ps_summary$trial_label <- factor(ps_summary$trial, levels = names(trial_labels), labels = unname(trial_labels))
ps_summary$strategy_label <- factor(strategy_labels[as.character(ps_summary$strategy)], levels = names(publication_strategy_colors))
sf07 <- ggplot2::ggplot(ps_summary, ggplot2::aes(propensity_score, colour = strategy_label, fill = strategy_label, linetype = strategy_label)) +
  ggplot2::geom_density(alpha = 0.14, linewidth = 0.8) +
  ggplot2::facet_wrap(~trial_label, ncol = 1) +
  ggplot2::scale_colour_manual(values = publication_strategy_colors) +
  ggplot2::scale_fill_manual(values = publication_strategy_colors) +
  ggplot2::scale_linetype_manual(values = publication_strategy_linetypes) +
  ggplot2::scale_x_continuous(limits = c(0, 1)) +
  ggplot2::labs(
    title = "Propensity-score overlap", subtitle = "Patient-level mean CBPS propensity score across 20 imputed datasets",
    x = "Propensity score for de-escalation", y = "Density", colour = NULL, fill = NULL, linetype = NULL
  ) + publication_theme(10)
register_labels("SF-06", "figure", c(levels(ps_summary$trial_label), levels(ps_summary$strategy_label)))
generated[["SF-06"]] <- save_figure(sf07, "SF-06", "SF-06_propensity_score_overlap", 8.5, 7.5)

weight_long$trial_label <- factor(weight_long$trial, levels = names(trial_labels), labels = unname(trial_labels))
weight_long$model_label <- factor(weight_long$model, levels = c("cbps_ate", "ato_brlogit"), labels = c("CBPS-ATE", "Overlap weighting (ATO)"))
weight_long$strategy_label <- factor(strategy_labels[as.character(weight_long$strategy)], levels = names(publication_strategy_colors))
sf08 <- ggplot2::ggplot(weight_long, ggplot2::aes(log10(weight), colour = strategy_label, fill = strategy_label, linetype = strategy_label)) +
  ggplot2::geom_density(alpha = 0.12, linewidth = 0.65) +
  ggplot2::facet_grid(trial_label ~ model_label, scales = "free_y") +
  ggplot2::scale_colour_manual(values = publication_strategy_colors) +
  ggplot2::scale_fill_manual(values = publication_strategy_colors) +
  ggplot2::scale_linetype_manual(values = publication_strategy_linetypes) +
  ggplot2::scale_x_continuous(breaks = log10(c(0.01, 0.03, 0.1, 0.3, 1, 3, 10)), labels = c("0.01", "0.03", "0.10", "0.30", "1", "3", "10")) +
  ggplot2::labs(
    title = "Distribution of analysis weights", subtitle = "All 20 imputed datasets; original untruncated weights",
    x = "Weight (log10 scale)", y = "Density", colour = NULL, fill = NULL, linetype = NULL
  ) + publication_theme(9)
register_labels("SF-07", "figure", c(levels(weight_long$trial_label), levels(weight_long$model_label), levels(weight_long$strategy_label)))
generated[["SF-07"]] <- save_figure(sf08, "SF-07", "SF-07_weight_distribution", 11.5, 8)

# SF-09: sensitivity analyses for the AKI risk difference.
sensitivity <- read_locked_csv("09_输出结果/06_敏感性分析/controlled_reporting/ST-08_sensitivity_percentile_ci.csv")
sensitivity <- sensitivity[sensitivity$estimand_field == "rd_aki", , drop = FALSE]
sensitivity$specification_label <- unname(c(
  main_cbps_ate = "Main analysis: CBPS-ATE",
  overlap_weight_ato_brlogit = "Overlap weighting with bias-reduced logistic regression",
  aki_cross_t0_uo_boundary = "Urine-output window allowed to cross index time",
  death_calendar_boundary = "Alternative death-date boundary",
  weight_truncation_1_99 = "ATE weights truncated at the 1st and 99th percentiles",
  weight_truncation_5_95 = "ATE weights truncated at the 5th and 95th percentiles",
  creatinine_only_phenotype = "Creatinine-only AKI phenotype",
  NEW_MS03 = "Alternative DAG adjustment set including Race",
  ABX_TIME_proxy = "Adjustment including the antibiotic-timing proxy",
  MICRO_STATE_proxy = "Microbiology-status proxy perturbation",
  complete_case = "Complete-case analysis"
)[sensitivity$specification])
if (anyNA(sensitivity$specification_label)) stop("SF-08 contains an unmapped sensitivity specification.")
sensitivity$trial_label <- factor(sensitivity$trial, levels = names(trial_labels), labels = unname(trial_labels))
spec_order <- unique(sensitivity$specification_label[sensitivity$trial == "anti_mrsa"])
sensitivity$specification_label <- factor(sensitivity$specification_label, levels = rev(spec_order))
sf09 <- ggplot2::ggplot(sensitivity, ggplot2::aes(point_estimate, specification_label, xmin = ci_lower, xmax = ci_upper)) +
  ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "#666666") +
  ggplot2::geom_errorbar(width = 0.16, orientation = "y", linewidth = 0.55, colour = "#333333") +
  ggplot2::geom_point(size = 2.3, shape = 18, colour = "#0072B2") +
  ggplot2::facet_wrap(~trial_label, ncol = 1) +
  ggplot2::labs(
    title = "Sensitivity analyses for 7-day AKI", subtitle = "Risk difference: de-escalation minus continuation",
    x = "Risk difference (95% bootstrap percentile interval)", y = NULL
  ) + publication_theme(9) +
  ggplot2::theme(panel.grid.major.y = ggplot2::element_blank())
register_labels("SF-08", "specification", c(as.character(spec_order), levels(sensitivity$trial_label)))
generated[["SF-08"]] <- save_figure(sf09, "SF-08", "SF-08_sensitivity_aki_risk_difference", 9.5, 8.5)

expected_ids <- c("MF-01", "MF-02", sprintf("SF-%02d", 1:8))
stopifnot(setequal(names(generated), expected_ids))
all_paths <- unlist(generated, use.names = FALSE)
if (!all(file.exists(all_paths) & file.info(all_paths)$size > 0)) stop("One or more publication figure files are missing or empty.")
label_audit <- do.call(rbind, label_audit_rows)
if (!all(label_audit$pass)) stop("Publication figure label audit failed.")
manifest <- do.call(rbind, lapply(names(generated), function(id) {
  paths <- generated[[id]]
  data.frame(
    artifact_id = id, format = c("png", "tiff", "pdf"),
    absolute_path = normalizePath(paths, winslash = "/", mustWork = TRUE), size_bytes = file.info(paths)$size,
    sha256 = vapply(paths, digest::digest, character(1), algo = "sha256", file = TRUE, serialize = FALSE),
    stringsAsFactors = FALSE
  )
}))
utils::write.csv(label_audit, file.path(audit_dir, "publication_figure_label_audit.csv"), row.names = FALSE, na = "", fileEncoding = "UTF-8")
utils::write.csv(manifest, file.path(audit_dir, "publication_figure_manifest.csv"), row.names = FALSE, na = "", fileEncoding = "UTF-8")
cat("PUBLICATION_FIGURE_IDS=", length(expected_ids), "\n", sep = "")
cat("PUBLICATION_FIGURE_FILES=", nrow(manifest), "\n", sep = "")
cat("PUBLICATION_FIGURE_GENERATION=PASS\n")
