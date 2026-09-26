#' Flexible and Batch Publication-Ready Bar and Lollipop Plots
#'
#' Generates publication-ready bar and lollipop plots from a data frame, CSV,
#' or Excel file. Supports frequency counting, summary metric aggregation (mean,
#' sum, identity), grouped/stacked arrangements, automated data labelling,
#' coordinate flipping, and multi-column batch plotting.
#'
#' @param data A \code{data.frame}, \code{matrix}, or character string file path
#'   to a \code{.csv}, \code{.xlsx}, or \code{.xls} file containing the survey data.
#' @param x Character vector; column name(s) in \code{data} representing the primary categorical
#'   grouping variable(s). If multiple variables are passed without \code{y}, a list of plots is returned.
#' @param y Optional character vector; column name(s) in \code{data} representing numeric response
#'   variable(s) to aggregate or plot directly. If multiple variables are passed, a list of plots is returned. Default is \code{NULL}.
#' @param group Optional character string; column name in \code{data} representing a secondary
#'   categorical factor for grouped or stacked layouts. Default is \code{NULL}.
#' @param type Character string; visual geometry style: \code{"bar"} (default) or \code{"lollipop"}. Default is \code{"bar"}.
#' @param stat Character string; aggregation statistic when \code{y} is provided:
#'   \code{"mean"} (default), \code{"sum"}, or \code{"identity"} (values plotted as-is). Default is \code{"mean"}.
#' @param position Character string; arrangement geometry: \code{"dodge"} (default for comparisons),
#'   \code{"stack"}, or \code{"fill"} (proportional 100\% stacked). For lollipop plots, \code{"dodge"} is always used. Default is \code{"dodge"}.
#' @param show_labels Logical; if \code{TRUE}, displays numeric value labels on or above the geometries. Default is \code{FALSE}.
#' @param horizontal Logical; if \code{TRUE}, flips Cartesian coordinates horizontally for legible category labels. Default is \code{FALSE}.
#' @param palette Character string; a valid \code{RColorBrewer} palette name (e.g., \code{"Set2"}, \code{"Dark2"}). Default is \code{"Set2"}.
#' @param bar_width Numeric; relative width of bars or spacing factor for dodging. Default is \code{0.7}.
#' @param point_size Numeric; size of the marker heads when \code{type = "lollipop"}. Default is \code{3.5}.
#' @param line_size Numeric; width of the stem lines when \code{type = "lollipop"}. Default is \code{0.8}.
#' @param xlab Optional character string; custom x-axis title. Default is \code{NULL}.
#' @param ylab Optional character string; custom y-axis title. Default is \code{NULL}.
#' @param title_prefix Optional character string; text prefix prepended to the plot title. Default is \code{""}.
#' @param ... Additional arguments passed to \code{plot_bar()} when calling \code{plot_lollipop()}.
#'
#' @return A single \code{ggplot2} object, or a named list of \code{ggplot2} objects
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
                     type = c("bar", "lollipop"),
                     stat = c("mean", "sum", "identity"),
                     position = c("dodge", "stack", "fill"),
                     show_labels = FALSE,
                     horizontal = FALSE,
                     palette = "Set2",
                     bar_width = 0.7,
                     point_size = 3.5,
                     line_size = 0.8,
                     xlab = NULL,
                     ylab = NULL,
                     title_prefix = "") {

  type <- match.arg(type)
  stat <- match.arg(stat)
  position <- match.arg(position)

  # Enforce dodge geometry for lollipop plots
  if (type == "lollipop" && position %in% c("stack", "fill")) {
    warning("Stacked and fill positions are not applicable to lollipop plots. Defaulting to 'dodge'.", call. = FALSE)
    position <- "dodge"
  }

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
  } else if (is.matrix(data) || inherits(data, "tbl_df") || inherits(data, "tbl")) {
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
        data = df, x = x[1], y = y_col, group = group, type = type, stat = stat,
        position = position, show_labels = show_labels, horizontal = horizontal,
        palette = palette, bar_width = bar_width, point_size = point_size,
        line_size = line_size, xlab = xlab, ylab = ylab, title_prefix = title_prefix
      )
    }
    return(plot_list)
  }

  if (is.null(y) && length(x) > 1) {
    plot_list <- vector("list", length(x))
    names(plot_list) <- x
    for (x_col in x) {
      plot_list[[x_col]] <- plot_bar(
        data = df, x = x_col, y = NULL, group = group, type = type, stat = stat,
        position = position, show_labels = show_labels, horizontal = horizontal,
        palette = palette, bar_width = bar_width, point_size = point_size,
        line_size = line_size, xlab = xlab, ylab = ylab, title_prefix = title_prefix
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

  # 5. Build ggplot Object & Geometries
  dodge_pos <- ggplot2::position_dodge(width = bar_width)
  effective_pos <- if (position == "dodge") dodge_pos else if (position == "fill") "fill" else "stack"

  p <- ggplot2::ggplot(
    calc_tbl,
    ggplot2::aes(x = .data[[target_x]], y = .data[["Metric_Val"]], fill = .data[[fill_var]])
  )

  if (type == "bar") {
    if (is.null(group)) {
      p <- p + ggplot2::geom_col(width = bar_width, color = "black", linewidth = 0.5, show.legend = FALSE)
    } else {
      p <- p + ggplot2::geom_col(position = effective_pos, width = bar_width, color = "black", linewidth = 0.5)
    }
  } else {
    # Lollipop geometry: uses geom_linerange for robust dodge alignment
    if (is.null(group)) {
      p <- p +
        ggplot2::geom_linerange(
          ggplot2::aes(ymin = 0, ymax = .data[["Metric_Val"]]),
          linewidth = line_size,
          color = "grey40"
        ) +
        ggplot2::geom_point(
          size = point_size,
          shape = 21,
          color = "black",
          stroke = 0.7,
          show.legend = FALSE
        )
    } else {
      p <- p +
        ggplot2::geom_linerange(
          ggplot2::aes(ymin = 0, ymax = .data[["Metric_Val"]], group = .data[[group]]),
          position = dodge_pos,
          linewidth = line_size,
          color = "grey40"
        ) +
        ggplot2::geom_point(
          ggplot2::aes(group = .data[[group]]),
          position = dodge_pos,
          size = point_size,
          shape = 21,
          color = "black",
          stroke = 0.7
        )
    }
  }

  # 6. Optional Numeric Labels
  if (show_labels && position != "fill") {
    if (type == "bar") {
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
    } else {
      # Lollipop labels positioned just beyond the marker head
      v_offset <- if (horizontal) 0.5 else -0.8
      h_offset <- if (horizontal) -0.4 else 0.5
      p <- p + ggplot2::geom_text(
        ggplot2::aes(label = .data[["Metric_Val"]], group = if (!is.null(group)) .data[[group]] else NULL),
        position = if (is.null(group)) ggplot2::position_identity() else dodge_pos,
        vjust = v_offset,
        hjust = h_offset,
        size = 3.2,
        fontface = "bold"
      )
    }
  }

  # 7. Aesthetics and Layout
  p <- p +
    ggplot2::scale_fill_brewer(palette = palette) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, if (show_labels) 0.15 else 0.08))) +
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

#' @rdname plot_bar
#' @export
plot_lollipop <- function(..., point_size = 3.5, line_size = 0.8) {
  plot_bar(..., type = "lollipop", point_size = point_size, line_size = line_size)
}
