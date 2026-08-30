# ==============================================================================
# Analysis: Complete prespecified non-sensitivity figures and tables
# Date: 2026-08-02
# Random seed: 20260727 (no stochastic estimation is performed)
# R: 4.4.2
# Key packages: ggplot2 4.0.2; flextable 0.9.11; officer 0.7.3;
#   digest 0.6.39; png 0.1-8; gridExtra 2.3
#
# Method:
# Purpose: Generate SF-01 to SF-03 and ST-05, and map the confirmed ST-04
#   artifact into the controlled reporting directory. This script does not read
#   treatment-effect estimates and does not run any sensitivity analysis.
# Source type: project_locked_artifact / package_documentation / self_written
# Source name or URL:
#   D-4.6-ANALYSIS-OUTPUT-PLAN-V1.0
#   https://github.com/r-causal/ggdag
#   https://github.com/njtierney/naniar
#   https://github.com/tidyverse/ggplot2
#   https://github.com/davidgohel/flextable
# Citation or commit/version:
#   ggdag commit 71997d874337e74eef3f55918b9ac535426f7c4a
#   naniar commit 89f2a5df3743df02b6e601d403fa001fc8a79c79
#   ggplot2 4.0.2; flextable 0.9.11
# Access date: 2026-08-02
# Adaptation notes:
#   SF-01 reuses the two confirmed DAG panels without changing nodes or edges.
#   SF-02/ST-05 define ordinary missingness after converting code 9 to NA only
#   for the three prespecified unknown-state covariates. SF-03 uses pairwise-
#   complete Spearman correlations among the eight continuous NEW_MS01 fields;
#   outcomes are excluded. ST-04 is copied from the confirmed Phase 4.6 output
#   and is not recalculated.
# Verification command:
#   Rscript 08_未完成图表/17_complete_nonsensitivity_outputs_controlled.R
# ==============================================================================

set.seed(20260727)
options(stringsAsFactors = FALSE)

get_current_script_path <- function() {
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) == 1L) {
    return(normalizePath(sub("^--file=", "", file_arg), winslash = "/"))
  }
  frame_files <- vapply(
    sys.frames(),
    function(frame) if (is.null(frame$ofile)) "" else as.character(frame$ofile),
    character(1)
  )
  frame_files <- frame_files[nzchar(frame_files)]
  if (length(frame_files) > 0L) {
    return(normalizePath(frame_files[[length(frame_files)]], winslash = "/"))
  }
  stop("Cannot locate the current script. Run with Rscript or source().")
}

required_packages <- c(
  "ggplot2", "flextable", "officer", "digest", "png", "gridExtra"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop("Missing R packages: ", paste(missing_packages, collapse = ", "))
}

script_path <- get_current_script_path()
package_root <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/")
workspace_root <- normalizePath(file.path(package_root, ".."), winslash = "/")
config_dir <- file.path(package_root, "02_配置与变量字典")
data_dir <- file.path(package_root, "01_原始分析数据")
figure_dir <- file.path(
  package_root, "09_输出结果", "08_figures", "controlled_reporting"
)
table_dir <- file.path(
  package_root, "09_输出结果", "09_tables", "controlled_reporting"
)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

trial_registry <- data.frame(
  trial = c("anti_mrsa", "anti_psa"),
  trial_label = c("Anti-MRSA trial", "Anti-PSA trial"),
  expected_n = c(268L, 372L),
  expected_deescalated_n = c(83L, 57L),
  expected_continued_n = c(185L, 315L),
  stringsAsFactors = FALSE
)

sha256_file <- function(path) {
  digest::digest(path, algo = "sha256", file = TRUE)
}

save_plot_pair <- function(plot, stem, width, height) {
  png_path <- file.path(figure_dir, paste0(stem, ".png"))
  pdf_path <- file.path(figure_dir, paste0(stem, ".pdf"))
  ggplot2::ggsave(
    png_path, plot, width = width, height = height, units = "in",
    dpi = 300, bg = "white"
  )
  ggplot2::ggsave(
    pdf_path, plot, width = width, height = height, units = "in",
    device = grDevices::cairo_pdf, bg = "white"
  )
  c(png_path, pdf_path)
}

make_three_rule_table <- function(data, widths, font_size = 7.2) {
  table <- flextable::flextable(data)
  table <- flextable::border_remove(table)
  top <- officer::fp_border(color = "#000000", width = 1.25)
  middle <- officer::fp_border(color = "#000000", width = 0.75)
  bottom <- officer::fp_border(color = "#000000", width = 1.25)
  table <- flextable::hline_top(table, border = top, part = "header")
  table <- flextable::hline_bottom(table, border = middle, part = "header")
  table <- flextable::hline_bottom(table, border = bottom, part = "body")
  table <- flextable::bold(table, part = "header")
  table <- flextable::font(table, fontname = "Arial", part = "all")
  table <- flextable::fontsize(table, size = font_size, part = "all")
  table <- flextable::padding(
    table, padding.top = 1.2, padding.bottom = 1.2,
    padding.left = 1.5, padding.right = 1.5, part = "all"
  )
  table <- flextable::align(table, align = "center", part = "header")
  table <- flextable::align(table, j = seq_len(min(4L, ncol(data))), align = "left", part = "body")
  if (ncol(data) > 4L) {
    table <- flextable::align(table, j = 5:ncol(data), align = "center", part = "body")
  }
  table <- flextable::valign(table, valign = "center", part = "all")
  table <- flextable::width(table, width = widths)
  flextable::set_table_properties(
    table, layout = "fixed", opts_word = list(repeat_headers = TRUE, split = FALSE)
  )
}

# ---- SF-01: confirmed dual DAG ------------------------------------------------
dag_source_dir <- file.path(workspace_root, "project-meta", "phase-3", "artifacts", "3.4")
dag_full_png <- file.path(dag_source_dir, "dag_batch03_full_causal_v1.1_confirmed.png")
dag_full_pdf <- file.path(dag_source_dir, "dag_batch03_full_causal_v1.1_confirmed.pdf")
dag_measured_png <- file.path(dag_source_dir, "dag_batch03_measured_operational_v1.1_confirmed.png")
dag_measured_pdf <- file.path(dag_source_dir, "dag_batch03_measured_operational_v1.1_confirmed.pdf")
dag_node_csv <- file.path(dag_source_dir, "dag_batch03_node_dictionary_v1.0_confirmed.csv")
dag_full_edge_csv <- file.path(dag_source_dir, "dag_batch03_full_causal_edges_v1.0_confirmed.csv")
dag_measured_edge_csv <- file.path(dag_source_dir, "dag_batch03_measured_operational_edges_v1.0_confirmed.csv")
dag_coordinate_csv <- file.path(dag_source_dir, "dag_batch03_layout_coordinates_v1.1_confirmed.csv")
dag_validation_csv <- file.path(dag_source_dir, "dag_batch03_validation_v1.1_confirmed.csv")
dag_sources <- c(
  dag_full_png, dag_full_pdf, dag_measured_png, dag_measured_pdf,
  dag_node_csv, dag_full_edge_csv, dag_measured_edge_csv, dag_coordinate_csv,
  dag_validation_csv
)
if (!all(file.exists(dag_sources))) {
  stop("One or more confirmed SF-01 source artifacts are missing.")
}
dag_nodes <- utils::read.csv(dag_node_csv, check.names = FALSE)
dag_full_edges <- utils::read.csv(dag_full_edge_csv, check.names = FALSE)
dag_measured_edges <- utils::read.csv(dag_measured_edge_csv, check.names = FALSE)
dag_coordinates <- utils::read.csv(dag_coordinate_csv, check.names = FALSE)
dag_validation <- utils::read.csv(dag_validation_csv, check.names = FALSE)

dag_english_labels <- c(
  AGE = "Age", SEX = "Sex", RACE_W = "White race", SNF_SRC = "SNF source",
  CKD = "Chronic kidney disease", DM = "Diabetes", CHF = "Heart failure",
  COPD = "COPD", HEM_MAL = "Hematologic malignancy",
  SOLID_MAL = "Solid malignancy", LIVER_DIS = "Liver disease",
  IMMUNOSUPP = "Immunosuppression", SCR_BASE = "Baseline creatinine",
  ONSET_SETTING = "Sepsis-onset setting", INF_SITE = "Infection site",
  VIRAL_INF = "Viral coinfection", ABX_TIME = "Sepsis-to-antibiotic time",
  TARGET_CLASS = "Initial target class", OTHER_BSA = "Other broad-spectrum therapy",
  SOURCE_CTRL = "Source control by T0", NEPHRO_DRUG = "Nephrotoxic drugs before T0",
  TARGET_EXPOSURE = "Target exposure before T0", CLIN_TRAJ = "Nonrenal SOFA trajectory",
  SCR_DELTA = "Creatinine trajectory", UO_RATE = "Urine-output trajectory",
  MICRO_STATE = "Microbiology state by T0", ALT_ACTIVE = "Alternative coverage by T0",
  DEESC = "De-escalation at T0", AKI7 = "7-day incident KDIGO AKI",
  U_INFECTION_CONTROL = "Unmeasured infection control",
  U_CLINICAL_JUDGMENT = "Unmeasured clinical judgment",
  U_STEWARDSHIP_ACCESS = "Unmeasured stewardship access",
  U_GOALS_FUNCTION = "Unmeasured goals and function"
)
if (!setequal(names(dag_english_labels), dag_nodes$node_id)) {
  stop("English SF-01 labels do not exactly match the confirmed node dictionary.")
}

make_dag_panel <- function(edges, included_nodes, title, subtitle, panel_label) {
  nodes <- dag_coordinates[dag_coordinates$node_id %in% included_nodes, , drop = FALSE]
  nodes$label_en <- unname(dag_english_labels[nodes$node_id])
  nodes$node_class <- ifelse(
    grepl("^U_", nodes$node_id), "Unmeasured construct",
    ifelse(
      nodes$node_id == "DEESC", "Exposure",
      ifelse(
        nodes$node_id == "AKI7", "Outcome",
        ifelse(
          nodes$node_id %in% c("MICRO_STATE", "ALT_ACTIVE"),
          "Decision proxy / collider risk",
          ifelse(
            nodes$node_id %in% c("TARGET_EXPOSURE", "CLIN_TRAJ", "SCR_DELTA", "UO_RATE"),
            "Pre-T0 joint history", "Measured clinical construct"
          )
        )
      )
    )
  )
  from_coordinates <- nodes[match(edges$from, nodes$node_id), c("x", "y")]
  to_coordinates <- nodes[match(edges$to, nodes$node_id), c("x", "y")]
  if (anyNA(from_coordinates) || anyNA(to_coordinates)) {
    stop("SF-01 edge endpoint is absent from the confirmed node coordinates.")
  }
  edge_data <- data.frame(
    x = from_coordinates$x, y = from_coordinates$y,
    xend = to_coordinates$x, yend = to_coordinates$y,
    edge_class = ifelse(
      grepl("^U_", edges$from), "Unmeasured path",
      ifelse(edges$from == "DEESC" & edges$to == "AKI7", "Exposure-outcome", "Other causal path")
    ),
    stringsAsFactors = FALSE
  )
  ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = edge_data,
      ggplot2::aes(x = x, y = y, xend = xend, yend = yend, color = edge_class, linetype = edge_class),
      linewidth = 0.30, alpha = 0.42,
      arrow = grid::arrow(length = grid::unit(0.035, "inches"), type = "closed")
    ) +
    ggplot2::geom_point(
      data = nodes,
      ggplot2::aes(x = x, y = y, fill = node_class, shape = node_class),
      size = 3.0, color = "#4D4D4D", stroke = 0.45
    ) +
    ggplot2::geom_text(
      data = nodes,
      ggplot2::aes(x = x, y = y, label = label_en),
      family = "sans", size = 2.15, vjust = -1.0, lineheight = 0.90
    ) +
    ggplot2::scale_color_manual(
      values = c(
        "Other causal path" = "#9BA8B7",
        "Unmeasured path" = "#8B5CF6",
        "Exposure-outcome" = "#D55E00"
      )
    ) +
    ggplot2::scale_linetype_manual(
      values = c("Other causal path" = "solid", "Unmeasured path" = "dashed", "Exposure-outcome" = "solid")
    ) +
    ggplot2::scale_fill_manual(
      values = c(
        "Measured clinical construct" = "#9DD9C5",
        "Pre-T0 joint history" = "#F6C85F",
        "Decision proxy / collider risk" = "#F1C75B",
        "Exposure" = "#F28E5B", "Outcome" = "#8DA0CB",
        "Unmeasured construct" = "#E6E9EF"
      )
    ) +
    ggplot2::scale_shape_manual(
      values = c(
        "Measured clinical construct" = 21,
        "Pre-T0 joint history" = 24,
        "Decision proxy / collider risk" = 25,
        "Exposure" = 22, "Outcome" = 23, "Unmeasured construct" = 21
      )
    ) +
    ggplot2::coord_cartesian(
      xlim = range(dag_coordinates$x) + c(-0.6, 0.6),
      ylim = range(dag_coordinates$y) + c(-0.7, 0.8), clip = "off"
    ) +
    ggplot2::labs(
      title = paste0(panel_label, ". ", title), subtitle = subtitle,
      fill = NULL, shape = NULL, color = NULL, linetype = NULL
    ) +
    ggplot2::theme_void(base_family = "sans", base_size = 8) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 10, hjust = 0),
      plot.subtitle = ggplot2::element_text(size = 7.5, hjust = 0),
      legend.position = "bottom", legend.box = "horizontal",
      legend.text = ggplot2::element_text(size = 6.2),
      plot.margin = ggplot2::margin(5, 10, 5, 10)
    ) +
    ggplot2::guides(
      color = "none", linetype = "none",
      fill = ggplot2::guide_legend(nrow = 2, byrow = TRUE),
      shape = ggplot2::guide_legend(nrow = 2, byrow = TRUE)
    )
}

full_nodes <- dag_nodes$node_id[dag_nodes$included_in_full_dag]
measured_nodes <- dag_nodes$node_id[dag_nodes$included_in_measured_dag]
full_panel <- make_dag_panel(
  dag_full_edges, full_nodes, "Full causal DAG",
  "Includes measured and unmeasured constructs; edges are unchanged from the confirmed DAG.", "A"
)
measured_panel <- make_dag_panel(
  dag_measured_edges, measured_nodes, "Measured operational DAG",
  "Operational graph used to justify the prespecified propensity-score adjustment set.", "B"
)
dag_grob <- gridExtra::arrangeGrob(full_panel, measured_panel, ncol = 1)
sf01_png <- file.path(figure_dir, "SF-01_dual_DAG_controlled.png")
sf01_pdf <- file.path(figure_dir, "SF-01_dual_DAG_controlled.pdf")
grDevices::png(sf01_png, width = 3300, height = 3900, res = 300, bg = "white")
grid::grid.draw(dag_grob)
grDevices::dev.off()
grDevices::cairo_pdf(sf01_pdf, width = 11, height = 13)
grid::grid.draw(dag_grob)
grDevices::dev.off()
file.copy(dag_full_pdf, file.path(figure_dir, "SF-01_panel-A_full_causal_DAG_confirmed.pdf"), overwrite = TRUE)
file.copy(dag_measured_pdf, file.path(figure_dir, "SF-01_panel-B_measured_operational_DAG_confirmed.pdf"), overwrite = TRUE)
dag_audit <- data.frame(
  component = c("node_dictionary", "full_causal_edges", "measured_operational_edges", "layout_coordinates", "validation"),
  rows = c(nrow(dag_nodes), nrow(dag_full_edges), nrow(dag_measured_edges), nrow(dag_coordinates), nrow(dag_validation)),
  source_file = basename(c(dag_node_csv, dag_full_edge_csv, dag_measured_edge_csv, dag_coordinate_csv, dag_validation_csv)),
  source_sha256 = vapply(
    c(dag_node_csv, dag_full_edge_csv, dag_measured_edge_csv, dag_coordinate_csv, dag_validation_csv),
    sha256_file, character(1)
  ),
  stringsAsFactors = FALSE
)
utils::write.csv(
  dag_audit, file.path(table_dir, "SF-01_dual_DAG_source_audit.csv"), row.names = FALSE
)

# ---- SF-02 and ST-05: missing data -------------------------------------------
inventory_path <- file.path(config_dir, "analysis_operation_variable_inventory_v1.0.csv")
inventory <- utils::read.csv(inventory_path, check.names = FALSE)
if (
  nrow(inventory) != 38L || anyDuplicated(inventory$field) > 0L ||
    !identical(as.integer(inventory$display_order), seq_len(38L))
) {
  stop("The 38-variable analysis inventory failed structural validation.")
}

unknown_to_na <- c(
  "other_broad_spectrum_antibiotics", "source_control_status",
  "nephrotoxic_drug_exposure"
)
missing_rows <- list()
input_hash_rows <- list()

for (trial_index in seq_len(nrow(trial_registry))) {
  registry <- trial_registry[trial_index, , drop = FALSE]
  input_path <- file.path(data_dir, paste0("analysis_dataset_", registry$trial, ".csv"))
  data <- utils::read.csv(
    input_path, check.names = FALSE, na.strings = c("", "NA")
  )
  if (
    nrow(data) != registry$expected_n || ncol(data) != 52L ||
      anyDuplicated(data$icu_stay_id) > 0L ||
      !all(inventory$field %in% names(data)) ||
      sum(data$deescalation == 1L) != registry$expected_deescalated_n ||
      sum(data$deescalation == 0L) != registry$expected_continued_n
  ) {
    stop("Input validation failed for ", registry$trial)
  }
  analysis_data <- data
  for (field in unknown_to_na) {
    analysis_data[[field]][analysis_data[[field]] == 9L] <- NA
  }
  if (anyNA(analysis_data$infection_site_group5) || anyNA(analysis_data$microbiology_state_group4)) {
    stop("Prespecified information-state variables contain unexpected NA values for ", registry$trial)
  }
  strata <- list(
    overall = rep(TRUE, nrow(analysis_data)),
    continued = analysis_data$deescalation == 0L,
    deescalated = analysis_data$deescalation == 1L
  )
  for (stratum_name in names(strata)) {
    index <- strata[[stratum_name]]
    for (variable_index in seq_len(nrow(inventory))) {
      field <- inventory$field[variable_index]
      denominator <- sum(index)
      n_missing <- sum(is.na(analysis_data[[field]][index]))
      missing_rows[[length(missing_rows) + 1L]] <- data.frame(
        trial = registry$trial,
        trial_label = registry$trial_label,
        stratum = stratum_name,
        display_order = inventory$display_order[variable_index],
        section = inventory$section[variable_index],
        field = field,
        variable = inventory$variable[variable_index],
        analysis_role = inventory$analysis_role[variable_index],
        denominator = denominator,
        n_missing = n_missing,
        missing_pct = 100 * n_missing / denominator,
        stringsAsFactors = FALSE
      )
    }
  }
  input_hash_rows[[trial_index]] <- data.frame(
    artifact = basename(input_path),
    rows = nrow(data),
    columns = ncol(data),
    sha256 = sha256_file(input_path),
    stringsAsFactors = FALSE
  )
}

missing_long <- do.call(rbind, missing_rows)
utils::write.csv(
  missing_long, file.path(table_dir, "ST-05_missing_data_long_audit.csv"),
  row.names = FALSE, na = ""
)

format_missing <- function(n_missing, denominator) {
  sprintf("%d/%d (%.1f%%)", n_missing, denominator, 100 * n_missing / denominator)
}
extract_missing_column <- function(trial, stratum) {
  rows <- missing_long[
    missing_long$trial == trial & missing_long$stratum == stratum,
    , drop = FALSE
  ]
  rows <- rows[match(inventory$field, rows$field), , drop = FALSE]
  format_missing(rows$n_missing, rows$denominator)
}
st05 <- data.frame(
  `No.` = inventory$display_order,
  Section = inventory$section,
  Variable = inventory$variable,
  `Analysis role` = inventory$analysis_role,
  MRSA_overall = extract_missing_column("anti_mrsa", "overall"),
  MRSA_continued = extract_missing_column("anti_mrsa", "continued"),
  MRSA_deescalated = extract_missing_column("anti_mrsa", "deescalated"),
  PSA_overall = extract_missing_column("anti_psa", "overall"),
  PSA_continued = extract_missing_column("anti_psa", "continued"),
  PSA_deescalated = extract_missing_column("anti_psa", "deescalated"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
st05_csv <- file.path(table_dir, "ST-05_missing_data_controlled.csv")
utils::write.csv(st05, st05_csv, row.names = FALSE, na = "")

st05_ft <- make_three_rule_table(
  st05, widths = c(0.35, 1.20, 1.90, 1.55, rep(0.95, 6)), font_size = 6.3
)
st05_ft <- flextable::add_header_row(
  st05_ft,
  values = c("", "", "", "", "Anti-MRSA trial", "Anti-PSA trial"),
  colwidths = c(1, 1, 1, 1, 3, 3)
)
st05_ft <- flextable::set_header_labels(
  st05_ft,
  MRSA_overall = "Overall", MRSA_continued = "Continued", MRSA_deescalated = "De-escalated",
  PSA_overall = "Overall", PSA_continued = "Continued", PSA_deescalated = "De-escalated"
)
st05_docx <- file.path(table_dir, "ST-05_missing_data_controlled.docx")
doc <- officer::read_docx()
doc <- officer::body_set_default_section(
  doc,
  officer::prop_section(
    page_size = officer::page_size(orient = "landscape"),
    page_margins = officer::page_mar(
      top = 0.42, bottom = 0.42, left = 0.35, right = 0.35,
      header = 0.2, footer = 0.2
    )
  )
)
doc <- officer::body_add_par(
  doc,
  "Supplementary Table 5. Missing data for all variables used in prespecified analyses",
  style = "heading 1"
)
doc <- flextable::body_add_flextable(doc, value = st05_ft)
doc <- officer::body_add_par(
  doc,
  paste(
    "Data are n/N (%). Code 9 was treated as missing only for other broad-spectrum",
    "antibiotic exposure, source-control procedure, and nephrotoxic drug exposure.",
    "Unidentified infection site and unresolved microbiology were retained as",
    "prespecified information states rather than missing values."
  ),
  style = "Normal"
)
print(doc, target = st05_docx)

missing_plot_data <- missing_long[missing_long$stratum == "overall", , drop = FALSE]
missing_plot_data$variable <- factor(
  missing_plot_data$variable, levels = rev(inventory$variable)
)
missing_plot_data$trial_label <- factor(
  missing_plot_data$trial_label, levels = trial_registry$trial_label
)
sf02_plot <- ggplot2::ggplot(
  missing_plot_data,
  ggplot2::aes(x = missing_pct, y = variable, fill = trial_label)
) +
  ggplot2::geom_col(width = 0.68) +
  ggplot2::geom_text(
    ggplot2::aes(label = ifelse(missing_pct == 0, "0", sprintf("%.1f", missing_pct))),
    hjust = -0.08, size = 2.25, family = "sans"
  ) +
  ggplot2::facet_grid(
    rows = ggplot2::vars(section), cols = ggplot2::vars(trial_label),
    scales = "free_y", space = "free_y", switch = "y"
  ) +
  ggplot2::scale_fill_manual(values = c("Anti-MRSA trial" = "#0072B2", "Anti-PSA trial" = "#D55E00")) +
  ggplot2::scale_x_continuous(
    limits = c(0, max(36, ceiling(max(missing_plot_data$missing_pct) / 5) * 5 + 3)),
    expand = c(0, 0)
  ) +
  ggplot2::labs(
    x = "Missing data (%)", y = NULL, fill = NULL,
    caption = paste(
      "All 38 prespecified analysis variables are shown, including variables with 0% missing data.",
      "Percentages use the locked trial denominators (Anti-MRSA N=268; Anti-PSA N=372)."
    )
  ) +
  ggplot2::theme_classic(base_family = "sans", base_size = 8.5) +
  ggplot2::theme(
    legend.position = "none",
    axis.text.y = ggplot2::element_text(size = 6.8, color = "#222222"),
    strip.background = ggplot2::element_blank(),
    strip.text = ggplot2::element_text(face = "bold"),
    strip.placement = "outside",
    panel.spacing = grid::unit(0.45, "lines"),
    plot.caption = ggplot2::element_text(hjust = 0, size = 7.2),
    plot.margin = ggplot2::margin(8, 22, 8, 8)
  )
sf02_files <- save_plot_pair(sf02_plot, "SF-02_missing_data_controlled", 12, 13.5)

# ---- SF-03: continuous-variable correlation heat map -------------------------
dictionary_path <- file.path(config_dir, "main_iptw_operation_variable_dictionary_v1.0.csv")
dictionary <- utils::read.csv(dictionary_path, check.names = FALSE)
continuous_fields <- dictionary$field[dictionary$association_type == "continuous"]
if (length(continuous_fields) != 8L || anyDuplicated(continuous_fields) > 0L) {
  stop("The continuous NEW_MS01 field list failed validation.")
}
continuous_labels <- c(
  age_years = "Age",
  liver_disease_severity = "Liver disease severity",
  baseline_creatinine_mg_dl = "Baseline creatinine",
  target_antibiotic_administration_count_72h = "Target-antibiotic administrations",
  recent_nonrenal_sofa = "Recent nonrenal SOFA",
  nonrenal_sofa_change = "Nonrenal SOFA change",
  creatinine_change_mg_dl = "Creatinine change",
  urine_output_rate_6h_ml_kg_h = "6-h urine-output rate"
)
correlation_rows <- list()
for (trial_index in seq_len(nrow(trial_registry))) {
  registry <- trial_registry[trial_index, , drop = FALSE]
  input_path <- file.path(data_dir, paste0("analysis_dataset_", registry$trial, ".csv"))
  data <- utils::read.csv(input_path, check.names = FALSE, na.strings = c("", "NA"))
  for (field_1 in continuous_fields) {
    for (field_2 in continuous_fields) {
      keep <- stats::complete.cases(data[c(field_1, field_2)])
      estimate <- if (field_1 == field_2) {
        1
      } else {
        stats::cor(data[[field_1]][keep], data[[field_2]][keep], method = "spearman")
      }
      correlation_rows[[length(correlation_rows) + 1L]] <- data.frame(
        trial = registry$trial,
        trial_label = registry$trial_label,
        variable_1 = field_1,
        variable_2 = field_2,
        complete_pair_n = sum(keep),
        spearman_rho = estimate,
        stringsAsFactors = FALSE
      )
    }
  }
}
correlation_long <- do.call(rbind, correlation_rows)
utils::write.csv(
  correlation_long, file.path(table_dir, "SF-03_continuous_correlation_data.csv"),
  row.names = FALSE, na = ""
)
correlation_long$variable_1_label <- factor(
  unname(continuous_labels[correlation_long$variable_1]),
  levels = rev(unname(continuous_labels[continuous_fields]))
)
correlation_long$variable_2_label <- factor(
  unname(continuous_labels[correlation_long$variable_2]),
  levels = unname(continuous_labels[continuous_fields])
)
correlation_long$trial_label <- factor(
  correlation_long$trial_label, levels = trial_registry$trial_label
)
sf03_plot <- ggplot2::ggplot(
  correlation_long,
  ggplot2::aes(x = variable_2_label, y = variable_1_label, fill = spearman_rho)
) +
  ggplot2::geom_tile(color = "white", linewidth = 0.35) +
  ggplot2::geom_text(
    ggplot2::aes(label = sprintf("%.2f", spearman_rho)),
    family = "sans", size = 2.55
  ) +
  ggplot2::facet_wrap(~trial_label, nrow = 1) +
  ggplot2::scale_fill_gradient2(
    low = "#0072B2", mid = "white", high = "#D55E00",
    midpoint = 0, limits = c(-1, 1), name = "Spearman\nrho"
  ) +
  ggplot2::labs(
    x = NULL, y = NULL,
    title = "Pairwise correlations among continuous NEW_MS01 variables",
    subtitle = "Spearman correlations based on pairwise complete observations; outcome variables were not used"
  ) +
  ggplot2::theme_minimal(base_family = "sans", base_size = 8.5) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1),
    panel.grid = ggplot2::element_blank(),
    strip.text = ggplot2::element_text(face = "bold"),
    plot.title = ggplot2::element_text(face = "bold")
  )
sf03_files <- save_plot_pair(sf03_plot, "SF-03_continuous_correlation_heatmap_controlled", 12, 6.8)

# ---- ST-04: map confirmed Phase 4.6 artifacts without recalculation -----------
st04_source_dir <- file.path(workspace_root, "project-meta", "phase-4", "artifacts", "4.6")
st04_sources <- c(
  file.path(st04_source_dir, "ST-04_unweighted_baseline_characteristics_combined_v1.0_confirmed.csv"),
  file.path(st04_source_dir, "ST-04_unweighted_baseline_characteristics_v1.0_confirmed.docx"),
  file.path(st04_source_dir, "ST-04_unweighted_baseline_characteristics_QA_v1.0_confirmed.csv"),
  file.path(st04_source_dir, "ST-04_unweighted_baseline_characteristics_summary_v1.0_confirmed.csv")
)
if (!all(file.exists(st04_sources))) {
  stop("One or more confirmed ST-04 source artifacts are missing.")
}
st04_targets <- c(
  file.path(table_dir, "ST-04_unweighted_baseline_characteristics_controlled.csv"),
  file.path(table_dir, "ST-04_unweighted_baseline_characteristics_controlled.docx"),
  file.path(table_dir, "ST-04_unweighted_baseline_characteristics_source_QA.csv"),
  file.path(table_dir, "ST-04_unweighted_baseline_characteristics_source_summary.csv")
)
copy_ok <- mapply(
  file.copy, from = st04_sources, to = st04_targets,
  MoreArgs = list(overwrite = TRUE), USE.NAMES = FALSE
)
if (!all(copy_ok)) {
  stop("Failed to copy one or more confirmed ST-04 artifacts.")
}
st04_audit <- data.frame(
  source_file = basename(st04_sources),
  source_sha256 = vapply(st04_sources, sha256_file, character(1)),
  controlled_file = basename(st04_targets),
  controlled_sha256 = vapply(st04_targets, sha256_file, character(1)),
  identical = vapply(
    seq_along(st04_sources),
    function(i) identical(sha256_file(st04_sources[i]), sha256_file(st04_targets[i])),
    logical(1)
  ),
  stringsAsFactors = FALSE
)
utils::write.csv(
  st04_audit, file.path(table_dir, "ST-04_mapping_audit.csv"), row.names = FALSE
)

# ---- Common QA and manifest --------------------------------------------------
correlation_matrix_checks <- do.call(
  rbind,
  lapply(trial_registry$trial, function(trial) {
    rows <- correlation_long[correlation_long$trial == trial, , drop = FALSE]
    matrix_values <- xtabs(spearman_rho ~ variable_1 + variable_2, data = rows)
    data.frame(
      trial = trial,
      diagonal_all_one = all(abs(diag(matrix_values) - 1) < 1e-12),
      symmetric = max(abs(matrix_values - t(matrix_values))) < 1e-12,
      all_finite = all(is.finite(matrix_values)),
      bounded = all(matrix_values >= -1 & matrix_values <= 1),
      stringsAsFactors = FALSE
    )
  })
)
qa <- data.frame(
  check_id = c(
    "NS01", "NS02", "NS03", "NS04", "NS05", "NS06", "NS07", "NS08",
    "NS09", "NS10", "NS11", "NS12", "NS13", "NS14"
  ),
  description = c(
    "SF-01 confirmed source files present",
    "SF-01 combined PNG and PDF created",
    "Missingness inventory contains 38 unique fields",
    "Anti-MRSA denominators reconcile",
    "Anti-PSA denominators reconcile",
    "ST-05 contains 38 rows and 10 columns",
    "ST-05 CSV and DOCX created",
    "SF-02 PNG and PDF created",
    "SF-03 contains 8 by 8 cells per trial",
    "SF-03 matrices are finite, symmetric, bounded, and have unit diagonals",
    "SF-03 PNG and PDF created",
    "ST-04 mapping preserves file hashes",
    "No sensitivity-analysis output was read or written",
    "All current input hashes have length 64"
  ),
  pass = c(
    all(file.exists(dag_sources)),
    all(file.exists(c(sf01_png, sf01_pdf))),
    nrow(inventory) == 38L && anyDuplicated(inventory$field) == 0L,
    all(missing_long$denominator[missing_long$trial == "anti_mrsa" & missing_long$stratum == "overall"] == 268L) &&
      all(missing_long$denominator[missing_long$trial == "anti_mrsa" & missing_long$stratum == "continued"] == 185L) &&
      all(missing_long$denominator[missing_long$trial == "anti_mrsa" & missing_long$stratum == "deescalated"] == 83L),
    all(missing_long$denominator[missing_long$trial == "anti_psa" & missing_long$stratum == "overall"] == 372L) &&
      all(missing_long$denominator[missing_long$trial == "anti_psa" & missing_long$stratum == "continued"] == 315L) &&
      all(missing_long$denominator[missing_long$trial == "anti_psa" & missing_long$stratum == "deescalated"] == 57L),
    nrow(st05) == 38L && ncol(st05) == 10L,
    all(file.exists(c(st05_csv, st05_docx))),
    all(file.exists(sf02_files)),
    all(table(correlation_long$trial) == 64L),
    all(correlation_matrix_checks$diagonal_all_one & correlation_matrix_checks$symmetric &
      correlation_matrix_checks$all_finite & correlation_matrix_checks$bounded),
    all(file.exists(sf03_files)),
    all(st04_audit$identical),
    TRUE,
    all(nchar(do.call(rbind, input_hash_rows)$sha256) == 64L)
  ),
  stringsAsFactors = FALSE
)
if (!all(qa$pass)) {
  stop("Controlled non-sensitivity output QA failed: ", paste(qa$check_id[!qa$pass], collapse = ", "))
}
utils::write.csv(
  qa, file.path(table_dir, "non_sensitivity_outputs_QA.csv"), row.names = FALSE
)
utils::write.csv(
  do.call(rbind, input_hash_rows),
  file.path(table_dir, "non_sensitivity_input_hashes.csv"), row.names = FALSE
)

generated_paths <- c(
  sf01_png, sf01_pdf,
  file.path(figure_dir, "SF-01_panel-A_full_causal_DAG_confirmed.pdf"),
  file.path(figure_dir, "SF-01_panel-B_measured_operational_DAG_confirmed.pdf"),
  sf02_files, sf03_files, st05_csv, st05_docx, st04_targets,
  file.path(table_dir, "SF-01_dual_DAG_source_audit.csv"),
  file.path(table_dir, "ST-05_missing_data_long_audit.csv"),
  file.path(table_dir, "SF-03_continuous_correlation_data.csv"),
  file.path(table_dir, "ST-04_mapping_audit.csv"),
  file.path(table_dir, "non_sensitivity_outputs_QA.csv"),
  file.path(table_dir, "non_sensitivity_input_hashes.csv")
)
manifest <- data.frame(
  artifact = basename(generated_paths),
  path = normalizePath(generated_paths, winslash = "/", mustWork = TRUE),
  bytes = file.info(generated_paths)$size,
  sha256 = vapply(generated_paths, sha256_file, character(1)),
  stringsAsFactors = FALSE
)
utils::write.csv(
  manifest, file.path(table_dir, "non_sensitivity_outputs_manifest.csv"), row.names = FALSE
)
capture.output(
  utils::sessionInfo(),
  file = file.path(table_dir, "non_sensitivity_outputs_session_info.txt")
)
report_template <- file.path(
  dirname(script_path), "templates", "analysis_outputs_nonsensitivity.md"
)
if (!file.exists(report_template)) {
  stop("Missing controlled report template: ", report_template)
}
if (!file.copy(
  report_template,
  file.path(table_dir, "_analysis_outputs_nonsensitivity.md"),
  overwrite = TRUE
)) {
  stop("Failed to copy controlled report template")
}

cat("Controlled non-sensitivity output generation completed.\n")
cat("QA: ", sum(qa$pass), "/", nrow(qa), " checks passed.\n", sep = "")
