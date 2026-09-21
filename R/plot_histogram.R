#' Publication-Ready Histogram and Density Plot
#'
#' Generates publication-ready frequency histograms for continuous fisheries
#' and environmental variables. Supports optional kernel density curves,
#' theoretical normal distribution overlays, and mean/median reference lines.
#'
#' @param data A \code{data.frame}, \code{matrix}, or file path to a \code{.csv} or \code{.xlsx} file.
#' @param x Character. Name of the continuous numeric variable to plot.
#' @param group Optional character. Categorical grouping variable for faceting or fill coloring.
#' @param bins Integer. Number of bins (default: 30).
#' @param binwidth Optional numeric. Explicit bin width (overrides \code{bins}).
#' @param add_density Logical. If \code{TRUE}, overlays an empirical kernel density curve.
#' @param add_normal Logical. If \code{TRUE}, overlays a theoretical normal distribution curve.
#' @param show_stats Logical. If \code{TRUE}, adds dashed reference lines for the mean and median.
#' @param facet Logical. If \code{TRUE} and \code{group} is provided, facets the plot by group.
#' @param palette Character. RColorBrewer palette name (default: \code{"Blues"}).
#' @param xlab Optional character. Custom x-axis label.
#' @param ylab Optional character. Custom y-axis label (default: \code{"Frequency"}).
#' @param title Optional character. Custom plot title.
#'
#' @return A \code{ggplot} object.
#' @export
#'
#' @import ggplot2
#' @importFrom stats na.omit sd dnorm median density
#' @importFrom tools file_ext
#' @importFrom utils read.csv
#' @importFrom rlang .data
plot_histogram <- function(data,
                           x,
                           group = NULL,
                           bins = 30,
                           binwidth = NULL,
                           add_density = FALSE,
                           add_normal = FALSE,
                           show_stats = TRUE,
                           facet = FALSE,
                           palette = "Blues",
                           xlab = NULL,
                           ylab = "Frequency",
                           title = NULL) {

  # 1. Standard Defensive Ingestion Guard
  if (is.character(data) && length(data) == 1) {
    if (!file.exists(data)) stop(sprintf("File '%s' not found.", data))
    ext <- tolower(tools::file_ext(data))
    if (ext == "csv") {
      df <- utils::read.csv(data, stringsAsFactors = FALSE)
    } else if (ext %in% c("xlsx", "xls")) {
      if (!requireNamespace("readxl", quietly = TRUE)) stop("Install 'readxl' to load Excel files.")
      df <- as.data.frame(readxl::read_excel(data))
    } else {
      stop("Unsupported file type. Use data.frame, matrix, .csv, or .xlsx.")
    }
  } else if (is.matrix(data)) {
    df <- as.data.frame(data)
  } else if (is.data.frame(data)) {
    df <- as.data.frame(data)
  } else {
    stop("'data' must be a data frame, matrix, or valid file path.")
  }

  if (!x %in% names(df)) stop(sprintf("Column '%s' not found.", x))
  if (!is.numeric(df[[x]])) stop(sprintf("Column '%s' must be numeric.", x))
  if (!is.null(group) && !group %in% names(df)) stop(sprintf("Group '%s' not found.", group))

  cols <- c(x, group)
  clean_df <- stats::na.omit(df[, cols, drop = FALSE])
  if (!is.null(group)) clean_df[[group]] <- as.factor(clean_df[[group]])

  # 2. Setup Base Plot (x only, no global y)
  use_density <- isTRUE(add_density) || isTRUE(add_normal)
  y_title <- if (use_density) "Density" else ylab

  p <- ggplot2::ggplot(clean_df, ggplot2::aes(x = .data[[x]]))

  # 3. Add Histogram Bars
  if (is.null(group)) {
    if (use_density) {
      p <- p + ggplot2::geom_histogram(
        ggplot2::aes(y = ggplot2::after_stat(density)),
        bins = bins,
        binwidth = binwidth,
        fill = "#4A90E2",
        color = "black",
        linewidth = 0.4,
        alpha = 0.75
      )
    } else {
      p <- p + ggplot2::geom_histogram(
        bins = bins,
        binwidth = binwidth,
        fill = "#4A90E2",
        color = "black",
        linewidth = 0.4,
        alpha = 0.75
      )
    }
  } else {
    if (use_density) {
      p <- p + ggplot2::geom_histogram(
        ggplot2::aes(y = ggplot2::after_stat(density), fill = .data[[group]]),
        bins = bins,
        binwidth = binwidth,
        color = "black",
        linewidth = 0.4,
        alpha = 0.7,
        position = if (facet) "identity" else "dodge"
      ) +
        ggplot2::scale_fill_brewer(palette = palette)
    } else {
      p <- p + ggplot2::geom_histogram(
        ggplot2::aes(fill = .data[[group]]),
        bins = bins,
        binwidth = binwidth,
        color = "black",
        linewidth = 0.4,
        alpha = 0.7,
        position = if (facet) "identity" else "dodge"
      ) +
        ggplot2::scale_fill_brewer(palette = palette)
    }
  }

  # 4. Add Density & Theoretical Normal Overlay
  if (add_density) {
    if (is.null(group) || !facet) {
      p <- p + ggplot2::geom_density(
        color = "#D9534F",
        linewidth = 1.0
      )
    } else {
      p <- p + ggplot2::geom_density(
        ggplot2::aes(color = .data[[group]]),
        linewidth = 1.0,
        show.legend = FALSE
      )
    }
  }

  if (add_normal && is.null(group)) {
    x_vals <- clean_df[[x]]
    mu <- mean(x_vals, na.rm = TRUE)
    sigma <- stats::sd(x_vals, na.rm = TRUE)
    p <- p + ggplot2::stat_function(
      fun = stats::dnorm,
      args = list(mean = mu, sd = sigma),
      color = "#2E7D32",
      linewidth = 1.0,
      linetype = "dashed"
    )
  }

  # 5. Reference Lines for Mean and Median
  if (show_stats && is.null(group)) {
    m_val <- mean(clean_df[[x]], na.rm = TRUE)
    med_val <- stats::median(clean_df[[x]], na.rm = TRUE)

    p <- p +
      ggplot2::geom_vline(
        xintercept = m_val,
        color = "darkred",
        linewidth = 0.8,
        linetype = "solid"
      ) +
      ggplot2::geom_vline(
        xintercept = med_val,
        color = "darkblue",
        linewidth = 0.8,
        linetype = "dashed"
      )
  }

  # 6. Faceting & Styling
  if (!is.null(group) && facet) {
    p <- p + ggplot2::facet_wrap(ggplot2::vars(.data[[group]]), scales = "free_y")
  }

  p <- p +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.1))) +
    ggplot2::labs(
      x = if (is.null(xlab)) x else xlab,
      y = y_title,
      title = title,
      fill = group
    ) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(
      axis.line = ggplot2::element_line(linewidth = 0.7, color = "black"),
      axis.text = ggplot2::element_text(color = "black", face = "bold"),
      axis.title = ggplot2::element_text(face = "bold"),
      legend.position = "top",
      legend.title = ggplot2::element_text(face = "bold"),
      strip.background = ggplot2::element_rect(fill = "grey92", color = "black"),
      strip.text = ggplot2::element_text(face = "bold")
    )

  return(p)
}
