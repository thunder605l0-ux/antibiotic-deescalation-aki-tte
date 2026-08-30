# Publication presentation helpers only. No statistical estimation is performed here.
# Date: 2026-08-08
# Source: project-locked display dictionary; ggplot2 and flextable documented APIs.
# Verification: source this file, load_publication_dictionary(), publication_three_rule_table().

publication_font_family <- "Times New Roman"
publication_strategy_colors <- c("Continuation" = "#0072B2", "De-escalation" = "#D55E00")
publication_strategy_linetypes <- c("Continuation" = "dashed", "De-escalation" = "solid")
publication_strategy_shapes <- c("Continuation" = 1L, "De-escalation" = 16L)

publication_project_root <- function() {
  current <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  while (!file.exists(file.path(current, "02_配置与变量字典", "variable_coding_table_en.csv"))) {
    parent <- dirname(current)
    if (identical(parent, current)) stop("Cannot locate project root.")
    current <- parent
  }
  current
}

load_publication_dictionary <- function(project_root = publication_project_root()) {
  path <- file.path(project_root, "02_配置与变量字典", "publication_display_labels_v1.0.csv")
  if (!file.exists(path)) stop("Missing publication display dictionary: ", path)
  dictionary <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE, fileEncoding = "UTF-8")
  if (nrow(dictionary) != 52L || anyDuplicated(dictionary$clean_field) || any(!nzchar(dictionary$display_label))) {
    stop("Publication display dictionary failed structural validation.")
  }
  dictionary
}

publication_label <- function(fields, dictionary = load_publication_dictionary(), with_unit = FALSE) {
  index <- match(fields, dictionary$clean_field)
  if (anyNA(index)) stop("Unmapped publication field(s): ", paste(unique(fields[is.na(index)]), collapse = ", "))
  labels <- dictionary$display_label[index]
  if (isTRUE(with_unit)) {
    unit <- dictionary$unit[index]
    labels[nzchar(unit)] <- paste0(labels[nzchar(unit)], ", ", unit[nzchar(unit)])
  }
  unname(labels)
}

publication_term_label <- function(terms, project_root = publication_project_root()) {
  path <- file.path(project_root, "02_配置与变量字典", "publication_model_term_labels_v1.0.csv")
  dictionary <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE, fileEncoding = "UTF-8")
  index <- match(terms, dictionary$term)
  if (anyNA(index)) stop("Unmapped publication model term(s): ", paste(unique(terms[is.na(index)]), collapse = ", "))
  unname(dictionary$display_label[index])
}

publication_theme <- function(base_size = 9) {
  ggplot2::theme_classic(base_size = base_size, base_family = publication_font_family) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = base_size + 1),
      plot.subtitle = ggplot2::element_text(size = base_size),
      axis.title = ggplot2::element_text(face = "plain"),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold"),
      legend.title = ggplot2::element_blank(),
      legend.position = "bottom",
      panel.spacing = grid::unit(1.1, "lines")
    )
}

publication_apply_strategy_scales <- function(plot) {
  plot +
    ggplot2::scale_colour_manual(values = publication_strategy_colors, drop = FALSE) +
    ggplot2::scale_fill_manual(values = publication_strategy_colors, drop = FALSE) +
    ggplot2::scale_linetype_manual(values = publication_strategy_linetypes, drop = FALSE) +
    ggplot2::scale_shape_manual(values = publication_strategy_shapes, drop = FALSE)
}

publication_three_rule_table <- function(data, font_size = 9, widths = NULL) {
  if (!requireNamespace("flextable", quietly = TRUE)) stop("Package 'flextable' is required.")
  if (!requireNamespace("officer", quietly = TRUE)) stop("Package 'officer' is required.")
  table <- flextable::flextable(data)
  table <- flextable::border_remove(table)
  top_border <- officer::fp_border(color = "000000", width = 1.25)
  mid_border <- officer::fp_border(color = "000000", width = 0.75)
  bottom_border <- officer::fp_border(color = "000000", width = 1.25)
  table <- flextable::hline_top(table, border = top_border, part = "header")
  table <- flextable::hline_bottom(table, border = mid_border, part = "header")
  if (nrow(data) > 0L) {
    table <- flextable::hline(table, i = nrow(data), border = bottom_border, part = "body")
  }
  table <- flextable::font(table, fontname = publication_font_family, part = "all")
  table <- flextable::fontsize(table, size = font_size, part = "all")
  table <- flextable::bold(table, bold = TRUE, part = "header")
  table <- flextable::align(table, align = "left", part = "all")
  table <- flextable::valign(table, valign = "center", part = "all")
  table <- flextable::padding(table, padding = 2.5, part = "all")
  if (!is.null(widths)) {
    if (length(widths) != ncol(data)) stop("widths must have one value per table column.")
    table <- flextable::width(table, j = seq_along(widths), width = widths)
  } else {
    table <- flextable::autofit(table)
  }
  table
}

publication_validate_reader_labels <- function(x) {
  x <- as.character(x)
  bad <- grepl("[A-Za-z0-9]+_[A-Za-z0-9_]+", x) |
    grepl("factor\\(", x) |
    grepl("race_white", x, fixed = TRUE) |
    grepl("White race", x, fixed = TRUE)
  data.frame(text = x, pass = !bad, stringsAsFactors = FALSE)
}
