#' Flexible and Batch Publication-Ready Bar Plots
#'
#' Generates publication-ready bar plots from a data frame, CSV, or Excel file.
#' Supports frequency counting, summary metric aggregation (mean, sum, identity),
#' grouped/stacked arrangements, automated data labelling, axis flipping, and
#' multi-column batch plotting.
#'
#' @param data A \code{data.frame}, \code{matrix}, or a character string path to a \code{.csv},
#'   \code{.xlsx}, or \code{.xls} file.
#' @param x Character vector. Name(s) of the categorical grouping variable(s).
#'   If multiple variables are passed without \code{y}, a list of frequency plots is returned.
#' @param y Optional character vector. Numeric variable(s) to aggregate or plot directly.
#'   If multiple variables are passed, a list of plots is returned.
#' @param group Optional character. Secondary categorical variable for grouped/stacked bars.
#' @param stat Character. Aggregation statistic when \code{y} is provided:
#'   \code{"mean"} (default), \code{"sum"}, or \code{"identity"} (values plotted as-is).
#' @param position Character. Bar positioning: \code{"dodge"} (default for comparisons),
#'   \code{"stack"}, or \code{"fill"} (proportional 100\% stacked).
#' @param show_labels Logical. If \code{TRUE}, adds value labels on or above the bars.
#' @param horizontal Logical. If \code{TRUE}, flips axes for long taxonomic or station names.
#' @param palette Character. An RColorBrewer palette name (e.g., \code{"Set2"}, \code{"Dark2"}, \code{"Blues"}).
#' @param bar_width Numeric. Width of bars (default: 0.7).
#' @param xlab Optional character string for x-axis title.
#' @param ylab Optional character string for y-axis title.
#' @param title_prefix Optional character string prefixed to the plot title.
#'
#' @return A single \code{ggplot} object, or a named list of \code{ggplot} objects
#'   if batch plotting across multiple columns.
#' @export
#'
#' @import ggplot2
#' @importFrom utils read.csv
#' @importFrom stats aggregate na.omit
#' @importFrom tools file_ext toTitleCase
#' @importFrom rlang .data
plot_bar <- function(data,
                     x,
                     y = NULL,
                     group = NULL,
                     stat = c("mean", "sum", "identity"),
                     position = c("dodge", "stack", "fill"),
                     show_labels = FALSE,
                     horizontal = FALSE,
                     palette = "Set2",
                     bar_width = 0.7,
                     xlab = NULL,
                     ylab = NULL,
                     title_prefix = "") {

  stat <- match.arg(stat)
  position <- match.arg(position)

  # 1. Ingest Data: data.frame, matrix, CSV, or Excel
  if (is.character(data) && length(data) == 1) {
    if (!file.exists(data)) stop(sprintf("The file '%s' does not exist.", data))
    ext <- tolower(tools::file_ext(data))
    if (ext == "csv") {
      df <- utils::read.csv(data, stringsAsFactors = FALSE, check.names = FALSE)
    } else if (ext %in% c("xlsx", "xls")) {
      if (!requireNamespace("readxl", quietly = TRUE)) {
        stop("Package 'readxl' is required to read Excel files. Please install it with install.packages('readxl').")
      }
      df <- as.data.frame(readxl::read_excel(data))
    } else {
      stop("Unsupported file type. Please provide a data frame, matrix, a .csv file, or an .xlsx/.xls file.")
    }
  } else if (is.matrix(data)) {
    df <- as.data.frame(data)
  } else if (is.data.frame(data)) {
    df <- as.data.frame(data)
  } else {
    stop("'data' must be a data.frame, matrix, or a file path to a .csv or .xlsx file.")
  }

  # 2. Handle Multi-Variable Batch Operations
  if (!is.null(y) && length(y) > 1) {
    plot_list <- vector("list", length(y))
    names(plot_list) <- y
    for (y_col in y) {
      plot_list[[y_col]] <- plot_bar(
        data = df, x = x[1], y = y_col, group = group, stat = stat,
        position = position, show_labels = show_labels, horizontal = horizontal,
        palette = palette, bar_width = bar_width, xlab = xlab, ylab = ylab,
        title_prefix = title_prefix
      )
    }
    return(plot_list)
  }

  if (is.null(y) && length(x) > 1) {
    plot_list <- vector("list", length(x))
    names(plot_list) <- x
    for (x_col in x) {
      plot_list[[x_col]] <- plot_bar(
        data = df, x = x_col, y = NULL, group = group, stat = stat,
        position = position, show_labels = show_labels, horizontal = horizontal,
        palette = palette, bar_width = bar_width, xlab = xlab, ylab = ylab,
        title_prefix = title_prefix
      )
    }
    return(plot_list)
  }

  # 3. Single Plot Input Verification
  target_x <- x[1]
  target_y <- if (!is.null(y)) y[1] else NULL

  if (!target_x %in% names(df)) stop(sprintf("Column '%s' was not found in the dataset.", target_x))
  if (!is.null(target_y) && !target_y %in% names(df)) stop(sprintf("Column '%s' was not found in the dataset.", target_y))
  if (!is.null(group) && !group %in% names(df)) stop(sprintf("Column '%s' was not found in the dataset.", group))

  needed_cols <- c(target_x, target_y, group)
  clean_df <- stats::na.omit(df[, needed_cols, drop = FALSE])
  clean_df[[target_x]] <- as.factor(clean_df[[target_x]])
  if (!is.null(group)) clean_df[[group]] <- as.factor(clean_df[[group]])

  # 4. Aggregate or Compute Frequencies
  if (is.null(target_y)) {
    # Frequency/Count workflow
    if (is.null(group)) {
      calc_tbl <- as.data.frame(table(clean_df[[target_x]]))
      colnames(calc_tbl) <- c(target_x, "Metric_Val")
      fill_var <- target_x
    } else {
      calc_tbl <- as.data.frame(table(clean_df[[target_x]], clean_df[[group]]))
      colnames(calc_tbl) <- c(target_x, group, "Metric_Val")
      fill_var <- group
    }
    display_y <- "Frequency (Count)"
  } else {
    # Continuous variable workflow
    if (!is.numeric(clean_df[[target_y]])) stop(sprintf("Column '%s' must be numeric.", target_y))

    if (stat == "identity") {
      calc_tbl <- clean_df
      calc_tbl$Metric_Val <- calc_tbl[[target_y]]
    } else {
      grp_list <- if (is.null(group)) list(clean_df[[target_x]]) else list(clean_df[[target_x]], clean_df[[group]])
      names(grp_list) <- if (is.null(group)) target_x else c(target_x, group)
      run_fn <- if (stat == "mean") mean else sum
      calc_tbl <- stats::aggregate(clean_df[[target_y]], by = grp_list, FUN = run_fn)
      colnames(calc_tbl)[ncol(calc_tbl)] <- "Metric_Val"
      calc_tbl$Metric_Val <- round(calc_tbl$Metric_Val, 2)
    }
    fill_var <- if (is.null(group)) target_x else group
    display_y <- if (stat == "identity") target_y else paste(tools::toTitleCase(stat), "of", target_y)
  }

  # 5. Build ggplot Object
  dodge_pos <- ggplot2::position_dodge(width = bar_width + 0.05)
  effective_pos <- if (position == "dodge") dodge_pos else if (position == "fill") "fill" else "stack"

  p <- ggplot2::ggplot(
    calc_tbl,
    ggplot2::aes(x = .data[[target_x]], y = .data[["Metric_Val"]], fill = .data[[fill_var]])
  )

  if (is.null(group)) {
    p <- p + ggplot2::geom_col(width = bar_width, color = "black", linewidth = 0.5, show.legend = FALSE)
  } else {
    p <- p + ggplot2::geom_col(position = effective_pos, width = bar_width, color = "black", linewidth = 0.5)
  }

  # Optional Numeric Labels
  if (show_labels && position != "fill") {
    if (position == "dodge" || is.null(group)) {
      v_offset <- if (horizontal) 0.5 else -0.4
      h_offset <- if (horizontal) -0.2 else 0.5
      p <- p + ggplot2::geom_text(
        ggplot2::aes(label = .data[["Metric_Val"]]),
        position = if (is.null(group)) ggplot2::position_identity() else dodge_pos,
        vjust = v_offset,
        hjust = h_offset,
        size = 3.3,
        fontface = "bold"
      )
    } else if (position == "stack") {
      p <- p + ggplot2::geom_text(
        ggplot2::aes(label = .data[["Metric_Val"]]),
        position = ggplot2::position_stack(vjust = 0.5),
        size = 3.1,
        color = "white",
        fontface = "bold"
      )
    }
  }

  # 6. Aesthetics and Layout
  p <- p +
    ggplot2::scale_fill_brewer(palette = palette) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.12))) +
    ggplot2::labs(
      x = if (is.null(xlab)) target_x else xlab,
      y = if (is.null(ylab)) display_y else ylab,
      title = if (title_prefix != "") paste(title_prefix, "-", display_y) else NULL,
      fill = group
    ) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(
      axis.line = ggplot2::element_line(linewidth = 0.7, color = "black"),
      axis.text = ggplot2::element_text(color = "black", face = "bold"),
      axis.title = ggplot2::element_text(face = "bold"),
      legend.position = if (is.null(group)) "none" else "top",
      legend.title = ggplot2::element_text(face = "bold")
    )

  if (horizontal) {
    p <- p + ggplot2::coord_flip()
  }

  return(p)
}
