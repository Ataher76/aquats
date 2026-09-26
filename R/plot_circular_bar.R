#' Publication-Grade Circular or Radial Bar Plot
#'
#' Builds circular bar charts with concentric reference gridlines, scale labels,
#' and category labels. Best suited for data with six or more categories.
#'
#' @param data A \code{data.frame}, \code{matrix}, or character string file path
#'   to a \code{.csv} or \code{.xlsx} file containing the category observations.
#' @param x Character string; column name in \code{data} representing the primary categorical
#'   grouping factor positioned radially around the circle.
#' @param y Optional character string; column name in \code{data} representing the numeric
#'   variable to aggregate. If \code{NULL}, observation frequencies (counts) are plotted. Default is \code{NULL}.
#' @param group Optional character string; column name in \code{data} representing a secondary
#'   categorical factor for bar fill coloring. Default is \code{NULL}.
#' @param stat Character string; aggregation statistic when \code{y} is provided:
#'   \code{"mean"} (default), \code{"sum"}, or \code{"identity"} (values plotted as-is). Default is \code{"mean"}.
#' @param inner_radius Numeric; fractional proportion of the plot radius reserved for the
#'   central inner void (between 0 and 1). Default is \code{0.35}.
#' @param show_grid Logical; if \code{TRUE}, renders concentric reference grid rings and
#'   numeric scale values. Default is \code{TRUE}.
#' @param palette Character string; a valid \code{RColorBrewer} palette name. Default is \code{"Set2"}.
#' @param title Optional character string; custom plot title. Default is \code{NULL}.
#'
#' @return A publication-ready \code{ggplot2} object.
#' @export
#'
#' @import ggplot2
#' @importFrom stats aggregate na.omit
#' @importFrom tools file_ext
#' @importFrom utils read.csv
#' @importFrom rlang .data
plot_circular_bar <- function(data,
                              x,
                              y = NULL,
                              group = NULL,
                              stat = c("mean", "sum", "identity"),
                              inner_radius = 0.35,
                              show_grid = TRUE,
                              palette = "Set2",
                              title = NULL) {

  stat <- match.arg(stat)

  # 1. Standard Defensive Ingestion Guard
  if (is.character(data) && length(data) == 1) {
    if (!file.exists(data)) stop(sprintf("File '%s' not found.", data))
    ext <- tolower(tools::file_ext(data))
    if (ext == "csv") {
      df <- utils::read.csv(data, stringsAsFactors = FALSE)
    } else if (ext %in% c("xlsx", "xls")) {
      if (!requireNamespace("readxl", quietly = TRUE)) stop("Please install 'readxl'.")
      df <- as.data.frame(readxl::read_excel(data))
    } else {
      stop("Unsupported file format. Provide a .csv or .xlsx file.")
    }
  } else if (is.matrix(data)) {
    df <- as.data.frame(data)
  } else if (is.data.frame(data)) {
    df <- as.data.frame(data)
  } else {
    stop("'data' must be a data frame, matrix, or valid file path.")
  }

  cols <- c(x, y, group)
  clean_df <- stats::na.omit(df[, cols, drop = FALSE])
  clean_df[[x]] <- as.factor(clean_df[[x]])
  if (!is.null(group)) clean_df[[group]] <- as.factor(clean_df[[group]])

  # 2. Aggregations
  fill_col <- if (!is.null(group)) group else x

  if (is.null(y)) {
    grp_by <- if (is.null(group)) list(clean_df[[x]]) else list(clean_df[[x]], clean_df[[group]])
    names(grp_by) <- if (is.null(group)) x else c(x, group)
    pdata <- as.data.frame(table(grp_by))
    colnames(pdata)[ncol(pdata)] <- "Value"
  } else {
    if (!is.numeric(clean_df[[y]])) stop(sprintf("Column '%s' must be numeric.", y))
    if (stat == "identity") {
      pdata <- clean_df
      pdata$Value <- pdata[[y]]
    } else {
      grp_by <- if (is.null(group)) list(clean_df[[x]]) else list(clean_df[[x]], clean_df[[group]])
      names(grp_by) <- if (is.null(group)) x else c(x, group)
      fn <- if (stat == "mean") mean else sum
      pdata <- stats::aggregate(clean_df[[y]], by = grp_by, FUN = fn)
      colnames(pdata)[ncol(pdata)] <- "Value"
    }
  }

  pdata$Value <- round(pdata$Value, 1)
  n_bars <- nrow(pdata)
  pdata$id <- seq_len(n_bars)

  # 3. Label Angles
  angle <- 90 - 360 * (pdata$id - 0.5) / n_bars
  pdata$hjust <- ifelse(angle < -90, 1, 0)
  pdata$angle <- ifelse(angle < -90, angle + 180, angle)

  # 4. Limits and Grid Radii
  max_val <- max(pdata$Value, na.rm = TRUE)
  min_val <- -max_val * inner_radius
  grid_ticks <- pretty(c(0, max_val), n = 3)
  grid_ticks <- grid_ticks[grid_ticks > 0 & grid_ticks <= max_val]

  # 5. Base Plot
  p <- ggplot2::ggplot(pdata, ggplot2::aes(x = .data$id, y = .data$Value, fill = .data[[fill_col]]))

  # Concentric Guide Rings
  if (show_grid) {
    for (tick in grid_ticks) {
      p <- p + ggplot2::geom_hline(
        yintercept = tick,
        color = "grey80",
        linewidth = 0.4,
        linetype = "dashed"
      )
    }
  }

  # Bars
  p <- p +
    ggplot2::geom_col(
      width = 0.72,
      color = "grey30",
      linewidth = 0.35,
      alpha = 0.9
    ) +
    ggplot2::coord_polar(start = 0) +
    ggplot2::ylim(min_val, max_val * 1.32) +
    ggplot2::scale_fill_brewer(palette = palette) +
    ggplot2::theme_void(base_size = 11) +
    ggplot2::theme(
      legend.position = if (is.null(group)) "none" else "bottom",
      legend.title = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 13)
    )

  # Grid Scale Labels
  if (show_grid) {
    p <- p + ggplot2::annotate(
      "text",
      x = rep(0.5, length(grid_ticks)),
      y = grid_ticks,
      label = as.character(grid_ticks),
      color = "grey40",
      size = 2.8,
      hjust = 1
    )
  }

  # Outer Labels
  p <- p + ggplot2::geom_text(
    ggplot2::aes(
      x = .data$id,
      y = .data$Value + (max_val * 0.05),
      label = paste0(.data[[x]], "\n(", .data$Value, ")"),
      hjust = .data$hjust,
      angle = .data$angle
    ),
    size = 3.0,
    fontface = "bold",
    inherit.aes = FALSE
  )

  if (!is.null(title)) p <- p + ggplot2::labs(title = title)

  return(p)
}
