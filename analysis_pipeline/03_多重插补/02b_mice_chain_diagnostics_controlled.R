# ==============================================================================
# Analysis: Controlled MICE chain-mean and chain-SD diagnostic figures
# Date: 2026-08-02
# Random seed: 20260727 (no stochastic estimation is performed)
# R: 4.4.2
# Key packages: ggplot2 4.0.2; digest 0.6.39
#
# Method:
# Purpose: Generate the prespecified SF-04 (chain means) and SF-05 (chain
#   standard deviations) directly from the saved mids chain statistics.
# Source type: GitHub package documentation / project_locked_object / self_written
# Source name or URL: https://github.com/amices/mice
# Citation or commit/version:
#   mice commit 61083e667fa41cc4bc655492591b100bf8a7e443;
#   saved project objects: m=20, maxit=50
# Access date: 2026-08-02
# Adaptation notes:
#   No imputation is rerun. Every actually imputed variable with finite saved
#   chain statistics is included, paginated for readability. Chain SD is the
#   square root of the saved chain variance. Output data are retained in CSV.
# Verification command:
#   Rscript 03_多重插补/02b_mice_chain_diagnostics_controlled.R
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

required_packages <- c("ggplot2", "digest")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop("Missing R packages: ", paste(missing_packages, collapse = ", "))
}

script_path <- get_current_script_path()
package_root <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/")
mice_dir <- file.path(package_root, "09_输出结果", "02_MICE")
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
  stringsAsFactors = FALSE
)

array_to_long <- function(chain_mean, chain_var, trial, trial_label) {
  if (!identical(dim(chain_mean), dim(chain_var))) {
    stop("chainMean and chainVar dimensions differ for ", trial)
  }
  if (any(chain_var < -1e-12, na.rm = TRUE)) {
    stop("Negative saved chain variance detected for ", trial)
  }
  dimensions <- dim(chain_mean)
  names_1 <- dimnames(chain_mean)[[1]]
  if (is.null(names_1) || length(names_1) != dimensions[1]) {
    stop("Missing variable names in saved chain statistics for ", trial)
  }
  grid <- expand.grid(
    variable = names_1,
    iteration = seq_len(dimensions[2]),
    chain = seq_len(dimensions[3]),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  grid$chain_mean <- as.vector(chain_mean)
  grid$chain_sd <- sqrt(pmax(as.vector(chain_var), 0))
  grid$trial <- trial
  grid$trial_label <- trial_label
  grid
}

make_trace_plot <- function(data, value_column, y_label, title, subtitle) {
  data$chain_factor <- factor(data$chain)
  variable_levels <- unique(gsub("_", " ", data$variable, fixed = TRUE))
  data$variable_label <- factor(
    gsub("_", " ", data$variable, fixed = TRUE), levels = variable_levels
  )
  plot <- ggplot2::ggplot(
    data,
    ggplot2::aes(
      x = iteration, y = .data[[value_column]],
      group = chain_factor, color = chain_factor
    )
  ) +
    ggplot2::geom_line(linewidth = 0.30, alpha = 0.60, na.rm = TRUE) +
    ggplot2::facet_wrap(~variable_label, scales = "free_y", ncol = 2, drop = FALSE) +
    ggplot2::scale_color_manual(values = grDevices::hcl.colors(20, palette = "Dynamic")) +
    ggplot2::scale_x_continuous(breaks = c(1, 10, 20, 30, 40, 50)) +
    ggplot2::labs(
      x = "Iteration", y = y_label, title = title, subtitle = subtitle,
      caption = "Each colored trajectory represents one of 20 imputation chains."
    ) +
    ggplot2::theme_classic(base_family = "Arial", base_size = 8.5) +
    ggplot2::theme(
      legend.position = "none",
      strip.background = ggplot2::element_rect(fill = "#F2F2F2", color = NA),
      strip.text = ggplot2::element_text(face = "bold", size = 7.5),
      plot.title = ggplot2::element_text(face = "bold"),
      plot.caption = ggplot2::element_text(hjust = 0, size = 7.2)
    )
  if (value_column == "chain_sd") {
    estimable <- tapply(is.finite(data$chain_sd), data$variable_label, any)
    nonestimable_levels <- names(estimable)[!estimable]
    if (length(nonestimable_levels) > 0L) {
      annotations <- data.frame(
        variable_label = factor(nonestimable_levels, levels = variable_levels),
        iteration = 25,
        chain_sd = 0,
        label = "Not estimable:\n<=1 imputed value per chain",
        stringsAsFactors = FALSE
      )
      plot <- plot + ggplot2::geom_text(
        data = annotations,
        ggplot2::aes(x = iteration, y = chain_sd, label = label),
        inherit.aes = FALSE,
        family = "Arial", size = 3, color = "#555555"
      )
    }
  }
  plot
}

all_qa <- list()
manifest_rows <- list()
variables_per_page <- 6L

for (trial_index in seq_len(nrow(trial_registry))) {
  registry <- trial_registry[trial_index, , drop = FALSE]
  object_path <- file.path(mice_dir, paste0("mids_", registry$trial, "_m20_maxit50.rds"))
  if (!file.exists(object_path)) {
    stop("Missing saved MICE object: ", object_path)
  }
  saved <- readRDS(object_path)
  imputation <- saved$imputation
  chain_mean <- imputation$chainMean
  chain_var <- imputation$chainVar
  if (
    !identical(saved$trial, registry$trial) || nrow(saved$raw_data) != registry$expected_n ||
      imputation$m != 20L || imputation$iteration != 50L ||
      !identical(as.integer(dim(chain_mean)), c(30L, 50L, 20L)) ||
      !identical(as.integer(dim(chain_var)), c(30L, 50L, 20L))
  ) {
    stop("Saved MICE object validation failed for ", registry$trial)
  }
  long <- array_to_long(
    chain_mean, chain_var, registry$trial, registry$trial_label
  )
  active_variables <- unique(long$variable[
    is.finite(long$chain_mean) | is.finite(long$chain_sd)
  ])
  method_targets <- names(imputation$method)[nzchar(imputation$method)]
  if (!setequal(active_variables, method_targets)) {
    stop("Saved finite chain variables differ from active imputation targets for ", registry$trial)
  }
  active_long <- long[long$variable %in% active_variables, , drop = FALSE]
  valid_chain_sd <- is.na(active_long$chain_sd) |
    (is.finite(active_long$chain_sd) & active_long$chain_sd >= 0)
  if (!all(is.finite(active_long$chain_mean)) || !all(valid_chain_sd)) {
    stop("Invalid active chain statistics detected for ", registry$trial)
  }
  csv_path <- file.path(
    table_dir, paste0("SF-04_SF-05_mice_chain_statistics_", registry$trial, ".csv")
  )
  utils::write.csv(active_long, csv_path, row.names = FALSE, na = "")

  pages <- split(
    active_variables,
    ceiling(seq_along(active_variables) / variables_per_page)
  )
  for (metric in c("chain_mean", "chain_sd")) {
    figure_id <- if (metric == "chain_mean") "SF-04" else "SF-05"
    y_label <- if (metric == "chain_mean") "Chain mean" else "Chain standard deviation"
    title <- paste0(
      if (metric == "chain_mean") "MICE chain-mean trajectories: " else "MICE chain-SD trajectories: ",
      registry$trial_label
    )
    pdf_path <- file.path(
      figure_dir,
      paste0(figure_id, "_mice_", metric, "_", registry$trial, "_controlled.pdf")
    )
    grDevices::cairo_pdf(pdf_path, width = 8.5, height = 10.5, onefile = TRUE)
    for (page_index in seq_along(pages)) {
      page_variables <- pages[[page_index]]
      page_data <- active_long[active_long$variable %in% page_variables, , drop = FALSE]
      page_data$variable <- factor(page_data$variable, levels = page_variables)
      plot <- make_trace_plot(
        page_data, metric, y_label, title,
        sprintf("Page %d of %d; m=20; 50 iterations", page_index, length(pages))
      )
      print(plot)
      png_path <- file.path(
        figure_dir,
        sprintf(
          "%s_mice_%s_%s_controlled_page%02d.png",
          figure_id, metric, registry$trial, page_index
        )
      )
      ggplot2::ggsave(
        png_path, plot, width = 8.5, height = 10.5, units = "in",
        dpi = 300, bg = "white"
      )
      manifest_rows[[length(manifest_rows) + 1L]] <- data.frame(
        artifact = basename(png_path),
        bytes = file.info(png_path)$size,
        sha256 = digest::digest(png_path, algo = "sha256", file = TRUE),
        stringsAsFactors = FALSE
      )
    }
    grDevices::dev.off()
    manifest_rows[[length(manifest_rows) + 1L]] <- data.frame(
      artifact = basename(pdf_path),
      bytes = file.info(pdf_path)$size,
      sha256 = digest::digest(pdf_path, algo = "sha256", file = TRUE),
      stringsAsFactors = FALSE
    )
  }
  manifest_rows[[length(manifest_rows) + 1L]] <- data.frame(
    artifact = basename(csv_path),
    bytes = file.info(csv_path)$size,
    sha256 = digest::digest(csv_path, algo = "sha256", file = TRUE),
    stringsAsFactors = FALSE
  )
  all_qa[[trial_index]] <- data.frame(
    trial = registry$trial,
    object_sha256 = digest::digest(object_path, algo = "sha256", file = TRUE),
    n = nrow(saved$raw_data),
    m = imputation$m,
    iterations = imputation$iteration,
    active_variable_n = length(active_variables),
    expected_cells = length(active_variables) * 50L * 20L,
    observed_cells = nrow(active_long),
    all_chain_means_finite = all(is.finite(active_long$chain_mean)),
    chain_sd_finite_n = sum(is.finite(active_long$chain_sd)),
    chain_sd_not_estimable_n = sum(is.na(active_long$chain_sd)),
    chain_sd_not_estimable_variable_n = sum(vapply(
      split(active_long$chain_sd, active_long$variable),
      function(x) all(is.na(x)), logical(1)
    )),
    all_chain_sds_valid = all(valid_chain_sd),
    passed = nrow(active_long) == length(active_variables) * 50L * 20L &&
      all(is.finite(active_long$chain_mean)) &&
      all(valid_chain_sd),
    stringsAsFactors = FALSE
  )
}

qa <- do.call(rbind, all_qa)
if (!all(qa$passed)) {
  stop("Controlled MICE chain-diagnostic QA failed.")
}
utils::write.csv(
  qa, file.path(table_dir, "SF-04_SF-05_mice_chain_diagnostics_QA.csv"),
  row.names = FALSE
)
manifest <- do.call(rbind, manifest_rows)
utils::write.csv(
  manifest, file.path(table_dir, "SF-04_SF-05_mice_chain_diagnostics_manifest.csv"),
  row.names = FALSE
)
capture.output(
  utils::sessionInfo(),
  file = file.path(table_dir, "SF-04_SF-05_mice_chain_diagnostics_session_info.txt")
)

cat("Controlled SF-04 and SF-05 generation completed.\n")
cat("QA: ", sum(qa$passed), "/", nrow(qa), " trials passed.\n", sep = "")
