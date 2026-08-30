# ============================================================
# Supplementary Figure 1: two-page portrait dual DAG
#
# Code source: 02_代码/R/01_研究设计与DAG/67_phase34_external_review_batch03_dual_dag_revision.R
# Citation/version: dag_batch03 v1.1 confirmed node dictionary, edge tables, and layout coordinates
# Access date: 2026-08-30
# Adaptation: preserve the confirmed graph topology, coordinates, roles, nodes and edges;
#             display stable short node codes in the graph and a complete code-to-label key below each panel;
#             export panel A on page 1 and panel B on page 2 in portrait orientation.
# Data/model status: display-only. No participant data, model fitting, bootstrap, or result calculation is run.
# Verification command: D:/software/R4.4.2/R-4.4.2/bin/x64/Rscript.exe manuscript/plot_english_dual_dag.R
# ============================================================

options(stringsAsFactors = FALSE, warn = 1)
suppressPackageStartupMessages({
  library(ggplot2)
  library(ggdag)
  library(dagitty)
  library(gridExtra)
})

# ---- 1. Locate locked inputs and candidate-output directory -----------------
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("Cannot determine the current script path.")
script_path <- normalizePath(sub("^--file=", "", script_arg[[1L]]), winslash = "/", mustWork = TRUE)
project_root <- normalizePath(Sys.getenv("TTE_PROJECT_ROOT", unset = "."), winslash = "/", mustWork = TRUE)
dag_dir <- file.path(project_root, "project-meta", "phase-3", "artifacts", "3.4")
out_dir <- file.path(project_root, "manuscript", "plot_new")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

nodes <- read.csv(file.path(dag_dir, "dag_batch03_node_dictionary_v1.0_confirmed.csv"), check.names = FALSE)
full_edges <- read.csv(file.path(dag_dir, "dag_batch03_full_causal_edges_v1.0_confirmed.csv"), check.names = FALSE)
measured_edges <- read.csv(file.path(dag_dir, "dag_batch03_measured_operational_edges_v1.0_confirmed.csv"), check.names = FALSE)
coordinates <- read.csv(file.path(dag_dir, "dag_batch03_layout_coordinates_v1.1_confirmed.csv"), check.names = FALSE)

# ---- 2. Stable display codes and complete English node dictionary -----------
node_codes <- c(
  AGE = "AGE", SEX = "SEX", RACE_W = "RACE", SNF_SRC = "SNF",
  CKD = "CKD", DM = "DM", CHF = "CHF", COPD = "COPD",
  HEM_MAL = "HEM-MAL", SOLID_MAL = "SOL-MAL", LIVER_DIS = "LIVER",
  IMMUNOSUPP = "IMMUNO", SCR_BASE = "SCr-BAS", ONSET_SETTING = "ONSET",
  INF_SITE = "SITE", VIRAL_INF = "VIRAL", ABX_TIME = "ABX-TIME",
  TARGET_CLASS = "T-CLASS", OTHER_BSA = "O-BSA",
  SOURCE_CTRL = "SRC-CTRL", NEPHRO_DRUG = "NEPHRO",
  TARGET_EXPOSURE = "T-EXP", CLIN_TRAJ = "SOFA-TR",
  SCR_DELTA = "SCr-TR", UO_RATE = "UO-TR", MICRO_STATE = "MICRO",
  ALT_ACTIVE = "ALT-ACT", DEESC = "DE-ESC", AKI7 = "AKI-7D",
  U_INFECTION_CONTROL = "U1", U_CLINICAL_JUDGMENT = "U2",
  U_STEWARDSHIP_ACCESS = "U3", U_GOALS_FUNCTION = "U4"
)

node_labels <- c(
  AGE = "Age", SEX = "Biological sex", RACE_W = "Race (White/non-White)",
  SNF_SRC = "Skilled nursing facility source", CKD = "Chronic kidney disease",
  DM = "Diabetes mellitus", CHF = "Congestive heart failure",
  COPD = "Chronic obstructive pulmonary disease", HEM_MAL = "Hematologic malignancy",
  SOLID_MAL = "Solid or metastatic malignancy", LIVER_DIS = "Liver disease",
  IMMUNOSUPP = "Immunosuppression", SCR_BASE = "Baseline serum creatinine",
  ONSET_SETTING = "Infection or sepsis onset setting", INF_SITE = "Infection site",
  VIRAL_INF = "Viral coinfection", ABX_TIME = "Time from sepsis onset to first target antibiotic",
  TARGET_CLASS = "Initial target-antibiotic class", OTHER_BSA = "Other broad-spectrum antibiotic coverage",
  SOURCE_CTRL = "Source control by T0", NEPHRO_DRUG = "Nephrotoxic co-medication before T0",
  TARGET_EXPOSURE = "Cumulative target-antibiotic exposure before T0",
  CLIN_TRAJ = "Nonrenal SOFA level and change before T0",
  SCR_DELTA = "Serum creatinine trajectory before T0",
  UO_RATE = "Urine-output trajectory before T0",
  MICRO_STATE = "Microbiological status by T0",
  ALT_ACTIVE = "Adequacy of alternative active coverage by T0",
  DEESC = "Removal of target coverage at T0",
  AKI7 = "Incident in-hospital KDIGO-AKI within 7 days after T0",
  U_INFECTION_CONTROL = "Unmeasured: true infection control or persistent inflammation",
  U_CLINICAL_JUDGMENT = "Unmeasured: clinical or prognostic judgement and related care",
  U_STEWARDSHIP_ACCESS = "Unmeasured: access to consultation and antimicrobial stewardship",
  U_GOALS_FUNCTION = "Unmeasured: treatment goals and baseline functional status"
)

role_labels <- c(
  "已测临床构念" = "Measured clinical construct",
  "T0前联合历史组成" = "Pre-T0 history component",
  "决策代理/碰撞风险" = "Decision proxy/collider risk",
  "暴露" = "Exposure",
  "主要结局" = "Primary outcome",
  "未测量构念" = "Unmeasured construct"
)

role_colours <- c(
  "Measured clinical construct" = "#A7D8C5",
  "Pre-T0 history component" = "#F6C85F",
  "Decision proxy/collider risk" = "#E9C46A",
  "Exposure" = "#F68B5B",
  "Primary outcome" = "#8DA0CB",
  "Unmeasured construct" = "#E5E7EB"
)
role_shapes <- c(
  "Measured clinical construct" = 21,
  "Pre-T0 history component" = 24,
  "Decision proxy/collider risk" = 25,
  "Exposure" = 22,
  "Primary outcome" = 23,
  "Unmeasured construct" = 21
)

if (!setequal(names(node_codes), nodes$node_id)) stop("Node-code dictionary does not exactly cover the confirmed node dictionary.")
if (anyDuplicated(unname(node_codes))) stop("Node display codes must be unique.")
if (!setequal(names(node_labels), nodes$node_id)) stop("English node-label dictionary does not exactly cover the confirmed node dictionary.")
if (!setequal(coordinates$node_id, nodes$node_id)) stop("Confirmed layout coordinates do not exactly cover the confirmed node dictionary.")

# ---- 3. Structural hard stops: display changes must never change the DAG -----
assert_graph <- function(edges, included_ids, expected_nodes, expected_edges, graph_name) {
  if (length(included_ids) != expected_nodes || nrow(edges) != expected_edges) stop(sprintf("%s count check failed.", graph_name))
  if (any(!edges$from %in% included_ids) || any(!edges$to %in% included_ids)) stop(sprintf("%s contains an edge endpoint outside its confirmed node set.", graph_name))
  if (!any(edges$from == "DEESC" & edges$to == "AKI7")) stop(sprintf("%s is missing the confirmed DEESC -> AKI7 edge.", graph_name))
  if (any(edges$from == "ABX_TIME" & edges$to == "AKI7")) stop(sprintf("%s contains the prohibited ABX_TIME -> AKI7 edge.", graph_name))
}

full_ids <- nodes$node_id[nodes$included_in_full_dag]
measured_ids <- nodes$node_id[nodes$included_in_measured_dag]
assert_graph(full_edges, full_ids, 33L, 170L, "Full causal DAG")
assert_graph(measured_edges, measured_ids, 29L, 153L, "Measured operational DAG")
if (any(grepl("^U_", measured_ids))) stop("Measured operational DAG must not contain unmeasured nodes.")
if (!all(paste(measured_edges$from, measured_edges$to, sep = "->") %in% paste(full_edges$from, full_edges$to, sep = "->"))) stop("Measured operational edges must be a subset of full causal edges.")

# ---- 4. Draw a graph with short, unique node codes ---------------------------
make_graph <- function(edges, included_ids, title_text, subtitle_text, y_upper) {
  node_data <- coordinates[match(included_ids, coordinates$node_id), c("node_id", "role", "x", "y")]
  node_data$code <- unname(node_codes[node_data$node_id])
  node_data$role_en <- factor(unname(role_labels[node_data$role]), levels = unname(role_labels))
  node_data$is_unmeasured <- grepl("^U_", node_data$node_id)
  node_data$label_y <- ifelse(
    node_data$node_id %in% c("DEESC", "AKI7"),
    node_data$y - 0.58,
    ifelse(node_data$is_unmeasured, node_data$y + 0.60, node_data$y + 0.34)
  )
  node_data$label_x <- node_data$x

  edge_data <- data.frame(from = edges$from, to = edges$to, stringsAsFactors = FALSE)
  edge_data$name <- edge_data$from
  edge_data$x <- node_data$x[match(edge_data$from, node_data$node_id)]
  edge_data$y <- node_data$y[match(edge_data$from, node_data$node_id)]
  edge_data$xend <- node_data$x[match(edge_data$to, node_data$node_id)]
  edge_data$yend <- node_data$y[match(edge_data$to, node_data$node_id)]
  if (anyNA(edge_data[, c("x", "y", "xend", "yend")])) stop("An edge endpoint is missing from the confirmed layout.")
  edge_data$latent_edge <- grepl("^U_", edge_data$from) | grepl("^U_", edge_data$to)
  edge_data$causal_edge <- edge_data$from == "DEESC" & edge_data$to == "AKI7"
  # The focal causal edge is drawn separately as a straight arrow. Shortening
  # both ends prevents the line and arrowhead from being hidden by the source
  # and outcome nodes, while a white halo separates it from dense background
  # edges. This is a display-only transformation; endpoints remain DEESC/AKI7.
  causal_data <- edge_data[edge_data$causal_edge, , drop = FALSE]
  if (nrow(causal_data) != 1L) stop("Exactly one DEESC -> AKI7 edge is required for focal-arrow rendering.")
  causal_dx <- causal_data$xend - causal_data$x
  causal_dy <- causal_data$yend - causal_data$y
  causal_distance <- sqrt(causal_dx^2 + causal_dy^2)
  causal_ux <- causal_dx / causal_distance
  causal_uy <- causal_dy / causal_distance
  causal_data$x_start <- causal_data$x + 0.30 * causal_ux
  causal_data$y_start <- causal_data$y + 0.30 * causal_uy
  causal_data$x_end <- causal_data$xend - 0.38 * causal_ux
  causal_data$y_end <- causal_data$yend - 0.38 * causal_uy

  ggplot() +
    ggdag::geom_dag_edges_fan(
      data = edge_data[!edge_data$latent_edge & !edge_data$causal_edge, , drop = FALSE],
      mapping = aes(x = x, y = y, xend = xend, yend = yend),
      spread = 0.54, n = 80, edge_colour = "#64748B", edge_width = 0.26, edge_alpha = 0.14
    ) +
    ggdag::geom_dag_edges_fan(
      data = edge_data[edge_data$latent_edge & !edge_data$causal_edge, , drop = FALSE],
      mapping = aes(x = x, y = y, xend = xend, yend = yend),
      spread = 0.62, n = 80, edge_colour = "#7C3AED", edge_width = 0.48,
      edge_alpha = 0.42, edge_linetype = "dashed"
    ) +
    geom_segment(
      data = causal_data,
      aes(x = x_start, y = y_start, xend = x_end, yend = y_end),
      inherit.aes = FALSE, colour = "white", linewidth = 2.25,
      lineend = "round"
    ) +
    geom_segment(
      data = causal_data,
      aes(x = x_start, y = y_start, xend = x_end, yend = y_end),
      inherit.aes = FALSE, colour = "#D55E00", linewidth = 1.15,
      lineend = "round",
      arrow = grid::arrow(type = "closed", length = grid::unit(0.12, "in"))
    ) +
    geom_point(
      data = node_data, aes(x = x, y = y, fill = role_en, shape = role_en),
      size = 4.0, colour = "#1F2937", stroke = 0.55
    ) +
    geom_label(
      data = node_data,
      aes(x = label_x, y = label_y, label = code),
      family = "sans", size = 1.85, fontface = "bold", colour = "#111827",
      fill = "#FFFFFFEB", linewidth = 0, label.padding = grid::unit(0.035, "lines"),
      lineheight = 0.90
    ) +
    scale_fill_manual(values = role_colours, drop = TRUE) +
    scale_shape_manual(values = role_shapes, drop = TRUE) +
    coord_fixed(ratio = 0.70, xlim = c(0.0, 13.45), ylim = c(0.0, y_upper), clip = "off") +
    ggdag::theme_dag() +
    guides(
      fill = guide_legend(nrow = 2, byrow = TRUE),
      shape = guide_legend(nrow = 2, byrow = TRUE)
    ) +
    theme(
      text = element_text(family = "sans"),
      legend.position = "bottom", legend.title = element_blank(),
      legend.text = element_text(size = 6.0), legend.key.size = grid::unit(0.14, "in"),
      legend.spacing.x = grid::unit(0.05, "in"),
      plot.title = element_text(face = "bold", size = 11.2, hjust = 0.5, margin = margin(b = 3)),
      plot.subtitle = element_text(size = 7.2, hjust = 0.5, margin = margin(b = 3)),
      plot.margin = margin(10, 18, 2, 18)
    ) +
    labs(
      title = title_text,
      subtitle = paste(strwrap(subtitle_text, width = 88L), collapse = "\n")
    )
}

# ---- 5. Build a three-column code-to-label key for each page ----------------
wrap_one <- function(x, width = 48L) paste(strwrap(x, width = width), collapse = "\n")

make_node_key <- function(included_ids) {
  n_cols <- 2L
  n_rows <- ceiling(length(included_ids) / n_cols)
  key <- data.frame(
    node_id = included_ids,
    order = seq_along(included_ids),
    stringsAsFactors = FALSE
  )
  key$column <- ((key$order - 1L) %/% n_rows) + 1L
  key$row <- ((key$order - 1L) %% n_rows) + 1L
  key$x0 <- c(0.10, 1.60)[key$column]
  key$y <- n_rows - key$row + 1L
  key$code <- unname(node_codes[key$node_id])
  key$label <- vapply(unname(node_labels[key$node_id]), wrap_one, character(1))
  key$role_en <- factor(
    unname(role_labels[nodes$role[match(key$node_id, nodes$node_id)]]),
    levels = unname(role_labels)
  )

  ggplot(key, aes(y = y)) +
    geom_point(aes(x = x0, fill = role_en, shape = role_en), size = 2.8, colour = "#1F2937", stroke = 0.45) +
    geom_text(aes(x = x0 + 0.055, label = code), hjust = 0, family = "sans", fontface = "bold", size = 2.05, colour = "#111827") +
    geom_text(aes(x = x0 + 0.39, label = label), hjust = 0, family = "sans", size = 1.82, lineheight = 0.90, colour = "#374151") +
    annotate("text", x = 0.03, y = n_rows + 0.90, label = "Node key: graph code - full construct name", hjust = 0, family = "sans", fontface = "bold", size = 2.45, colour = "#111827") +
    scale_fill_manual(values = role_colours, drop = TRUE) +
    scale_shape_manual(values = role_shapes, drop = TRUE) +
    coord_cartesian(xlim = c(0.0, 3.0), ylim = c(0.4, n_rows + 1.05), clip = "off") +
    theme_void() +
    theme(plot.margin = margin(0, 18, 0, 18), legend.position = "none")
}

make_footer <- function(text) {
  grid::textGrob(
    paste(strwrap(text, width = 118L), collapse = "\n"),
    x = grid::unit(0.02, "npc"), y = grid::unit(0.50, "npc"),
    just = c("left", "center"),
    gp = grid::gpar(fontfamily = "sans", fontsize = 6.0, col = "#374151")
  )
}

full_graph <- make_graph(
  full_edges, full_ids,
  "A. Full causal DAG",
  "T0 removal of target coverage and incident in-hospital KDIGO-AKI within 7 days",
  13.95
)
measured_graph <- make_graph(
  measured_edges, measured_ids,
  "B. Measured operational DAG",
  "Variables derivable from the structured database for the target-trial emulation",
  12.40
)

full_page <- gridExtra::arrangeGrob(
  full_graph,
  make_node_key(full_ids),
  make_footer("Light-gray nodes and purple dashed edges denote unmeasured constructs. The orange arrow highlights the prespecified DE-ESC to AKI-7D contrast. All graph codes correspond one-to-one to the full labels in the node key."),
  ncol = 1,
  heights = c(0.58, 0.36, 0.06)
)
measured_page <- gridExtra::arrangeGrob(
  measured_graph,
  make_node_key(measured_ids),
  make_footer("The measured operational DAG omits U1-U4 because these constructs are not derivable from structured data. The adjustment set is sufficient only if the assumptions encoded in this operational DAG hold."),
  ncol = 1,
  heights = c(0.61, 0.33, 0.06)
)

# ---- 6. Export a two-page portrait PDF and page-specific raster files --------
output_stem <- file.path(out_dir, "SF01_dual_DAG_english_confirmed_layout")

grDevices::cairo_pdf(paste0(output_stem, ".pdf"), width = 6.7, height = 9.5, onefile = TRUE, bg = "white")
grid::grid.newpage(); grid::grid.draw(full_page)
grid::grid.newpage(); grid::grid.draw(measured_page)
grDevices::dev.off()

ggsave(paste0(output_stem, "_pageA.png"), full_page, width = 6.7, height = 9.5, units = "in", dpi = 320, bg = "white")
ggsave(paste0(output_stem, "_pageB.png"), measured_page, width = 6.7, height = 9.5, units = "in", dpi = 320, bg = "white")
ggsave(paste0(output_stem, "_pageA.tiff"), full_page, width = 6.7, height = 9.5, units = "in", dpi = 600, device = "tiff", compression = "lzw", bg = "white")
ggsave(paste0(output_stem, "_pageB.tiff"), measured_page, width = 6.7, height = 9.5, units = "in", dpi = 600, device = "tiff", compression = "lzw", bg = "white")

# Tall A-over-B previews retain the historic one-file raster names for audit.
dual_preview <- gridExtra::arrangeGrob(full_page, measured_page, ncol = 1, heights = c(1, 1))
ggsave(paste0(output_stem, ".png"), dual_preview, width = 6.7, height = 19.0, units = "in", dpi = 320, bg = "white")
ggsave(paste0(output_stem, ".tiff"), dual_preview, width = 6.7, height = 19.0, units = "in", dpi = 600, device = "tiff", compression = "lzw", bg = "white")

message("Created two-page portrait Supplementary Figure 1 in: ", out_dir)
