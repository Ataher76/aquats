#' Publication-Ready Temporal Tracking and Factorial Time-Series Plots
#'
#' Generates publication-ready time-series and line graphs for aquatic,
#' limnological, and fisheries datasets. Supports single-factor trajectories,
#' two-factor comparisons, and three-factor nested layouts with automated
#' replicate aggregation (Mean \eqn{\pm} SE, SD, or 95% CI ribbons) and faceting.
#'
#' @param data A \code{data.frame}, \code{matrix}, or file path to a \code{.csv} or \code{.xlsx} file.
#' @param time_var Character. Column name representing time steps (Date, POSIXct, numeric year, or factor).
#' @param y_var Character. Column name of the continuous response variable.
#' @param color_var Optional character. Categorical variable mapped to line/point colors (2nd factor).
#' @param linetype_var Optional character. Categorical variable mapped to line types (3rd factor).
#' @param facet_var Optional character. Categorical variable used to facet plots into panels.
#' @param ribbon Character. Uncertainty band style when replicates exist per time point:
#'   \code{"ci"} (95% Confidence Interval, default), \code{"se"} (\eqn{\pm 1} Standard Error),
#'   \code{"sd"} (\eqn{\pm 1} Standard Deviation), or \code{"none"}.
#' @param show_points Logical. If \code{TRUE} (default), renders points at each observation/mean time point.
#' @param point_size Numeric. Size of data points (default: 2.5).
#' @param line_width Numeric. Thickness of the trajectory lines (default: 1.0).
#' @param ribbon_alpha Numeric. Opacity of uncertainty ribbons (default: 0.20).
#' @param smooth Logical. If \code{TRUE}, adds a LOESS or spline smoothing trendline (default: \code{FALSE}).
#' @param color_palette Character. RColorBrewer palette name (default: \code{"Dark2"}).
#' @param xlab Optional character. Custom x-axis label.
#' @param ylab Optional character. Custom y-axis label.
#' @param title Optional character. Plot title.
#'
#' @return A list containing:
#'   \item{Plot}{The publication-ready \code{ggplot} object.}
#'   \item{Summary_Data}{A tidy data frame of aggregated Means, SD, SE, and Confidence Intervals.}
#' @export
#'
#' @import ggplot2
#' @importFrom stats aggregate na.omit sd qt as.formula
#' @importFrom tools file_ext
#' @importFrom utils read.csv
#' @importFrom rlang .data
plot_timeseries <- function(data,
                            time_var,
                            y_var,
                            color_var = NULL,
                            linetype_var = NULL,
                            facet_var = NULL,
                            ribbon = c("ci", "se", "sd", "none"),
                            show_points = TRUE,
                            point_size = 2.5,
                            line_width = 1.0,
                            ribbon_alpha = 0.20,
                            smooth = FALSE,
                            color_palette = "Dark2",
                            xlab = NULL,
                            ylab = NULL,
                            title = NULL) {

  ribbon <- match.arg(ribbon)

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
      stop("Unsupported file type. Use .csv, .xlsx, or data.frame.")
    }
  } else if (is.matrix(data)) {
    df <- as.data.frame(data)
  } else if (is.data.frame(data)) {
    df <- as.data.frame(data)
  } else {
    stop("Input 'data' must be a data.frame, matrix, or valid file path.")
  }

  required_cols <- c(time_var, y_var, color_var, linetype_var, facet_var)
  missing_cols <- required_cols[!required_cols %in% names(df)]
  if (length(missing_cols) > 0) {
    stop(sprintf("Specified column(s) not found in data: %s", paste(missing_cols, collapse = ", ")))
  }

  if (!is.numeric(df[[y_var]])) {
    stop(sprintf("Response column '%s' must be numeric.", y_var))
  }

  clean_df <- stats::na.omit(df[, required_cols, drop = FALSE])
  if (nrow(clean_df) < 2) {
    stop("Insufficient complete observations to construct a time-series plot.")
  }

  # Preserve Time Structure
  is_temporal_date <- inherits(clean_df[[time_var]], c("Date", "POSIXt"))
  if (!is_temporal_date && is.character(clean_df[[time_var]])) {
    clean_df[[time_var]] <- factor(clean_df[[time_var]], levels = unique(clean_df[[time_var]]))
  }

  if (!is.null(color_var)) clean_df[[color_var]] <- as.factor(clean_df[[color_var]])
  if (!is.null(linetype_var)) clean_df[[linetype_var]] <- as.factor(clean_df[[linetype_var]])
  if (!is.null(facet_var)) clean_df[[facet_var]] <- as.factor(clean_df[[facet_var]])

  # 2. Replicate Aggregation
  grp_factors <- c(time_var, color_var, linetype_var, facet_var)
  fmla <- stats::as.formula(paste(y_var, "~", paste(grp_factors, collapse = " + ")))

  n_fun <- function(z) sum(!is.na(z))
  sd_fun <- function(z) if (length(z) > 1) stats::sd(z, na.rm = TRUE) else 0
  se_fun <- function(z) {
    n <- sum(!is.na(z))
    if (n > 1) stats::sd(z, na.rm = TRUE) / sqrt(n) else 0
  }
  lci_fun <- function(z) {
    n <- sum(!is.na(z))
    m <- mean(z, na.rm = TRUE)
    if (n > 1) m - stats::qt(0.975, df = n - 1) * (stats::sd(z, na.rm = TRUE) / sqrt(n)) else m
  }
  uci_fun <- function(z) {
    n <- sum(!is.na(z))
    m <- mean(z, na.rm = TRUE)
    if (n > 1) m + stats::qt(0.975, df = n - 1) * (stats::sd(z, na.rm = TRUE) / sqrt(n)) else m
  }

  df_n    <- stats::aggregate(fmla, data = clean_df, FUN = n_fun)
  df_mean <- stats::aggregate(fmla, data = clean_df, FUN = mean)
  df_sd   <- stats::aggregate(fmla, data = clean_df, FUN = sd_fun)
  df_se   <- stats::aggregate(fmla, data = clean_df, FUN = se_fun)
  df_lci  <- stats::aggregate(fmla, data = clean_df, FUN = lci_fun)
  df_uci  <- stats::aggregate(fmla, data = clean_df, FUN = uci_fun)

  summary_df <- df_mean
  val_col_idx <- ncol(summary_df)
  colnames(summary_df)[val_col_idx] <- "Mean"

  summary_df$N   <- df_n[[val_col_idx]]
  summary_df$SD  <- df_sd[[val_col_idx]]
  summary_df$SE  <- df_se[[val_col_idx]]
  summary_df$LCI <- df_lci[[val_col_idx]]
  summary_df$UCI <- df_uci[[val_col_idx]]

  # Re-apply classes
  if (is_temporal_date) {
    summary_df[[time_var]] <- as.Date(summary_df[[time_var]])
  } else if (is.numeric(clean_df[[time_var]])) {
    summary_df[[time_var]] <- as.numeric(summary_df[[time_var]])
  } else {
    summary_df[[time_var]] <- factor(summary_df[[time_var]], levels = levels(clean_df[[time_var]]))
  }

  if (!is.null(color_var)) summary_df[[color_var]] <- factor(summary_df[[color_var]], levels = levels(clean_df[[color_var]]))
  if (!is.null(linetype_var)) summary_df[[linetype_var]] <- factor(summary_df[[linetype_var]], levels = levels(clean_df[[linetype_var]]))
  if (!is.null(facet_var)) summary_df[[facet_var]] <- factor(summary_df[[facet_var]], levels = levels(clean_df[[facet_var]]))

  # Assign ribbon limits
  if (ribbon == "ci") {
    summary_df$ymin <- summary_df$LCI
    summary_df$ymax <- summary_df$UCI
  } else if (ribbon == "se") {
    summary_df$ymin <- summary_df$Mean - summary_df$SE
    summary_df$ymax <- summary_df$Mean + summary_df$SE
  } else if (ribbon == "sd") {
    summary_df$ymin <- summary_df$Mean - summary_df$SD
    summary_df$ymax <- summary_df$Mean + summary_df$SD
  } else {
    summary_df$ymin <- summary_df$Mean
    summary_df$ymax <- summary_df$Mean
  }

  # Build unified grouping string to guarantee connected lines across factors
  extra_factors <- c(color_var, linetype_var, facet_var)
  if (length(extra_factors) == 0) {
    summary_df$Group_Line <- factor(1)
  } else {
    summary_df$Group_Line <- interaction(summary_df[, extra_factors, drop = FALSE], sep = " - ")
  }

  # 3. Base ggplot Assembly
  p <- ggplot2::ggplot(
    summary_df,
    ggplot2::aes(
      x = .data[[time_var]],
      y = .data$Mean,
      group = .data$Group_Line
    )
  )

  has_replicates <- any(summary_df$N > 1)

  # Ribbon Layer
  if (ribbon != "none" && has_replicates) {
    if (is.null(color_var)) {
      p <- p + ggplot2::geom_ribbon(
        ggplot2::aes(ymin = .data$ymin, ymax = .data$ymax),
        fill = "#2A9D8F",
        alpha = ribbon_alpha,
        color = NA
      )
    } else {
      p <- p + ggplot2::geom_ribbon(
        ggplot2::aes(ymin = .data$ymin, ymax = .data$ymax, fill = .data[[color_var]]),
        alpha = ribbon_alpha,
        color = NA
      )
    }
  }

  # Line Layer
  if (is.null(color_var) && is.null(linetype_var)) {
    p <- p + ggplot2::geom_line(linewidth = line_width, color = "#2A9D8F")
  } else if (!is.null(color_var) && is.null(linetype_var)) {
    p <- p + ggplot2::geom_line(ggplot2::aes(color = .data[[color_var]]), linewidth = line_width)
  } else if (is.null(color_var) && !is.null(linetype_var)) {
    p <- p + ggplot2::geom_line(ggplot2::aes(linetype = .data[[linetype_var]]), linewidth = line_width, color = "#2A9D8F")
  } else {
    p <- p + ggplot2::geom_line(ggplot2::aes(color = .data[[color_var]], linetype = .data[[linetype_var]]), linewidth = line_width)
  }

  # Smooth Layer
  if (smooth) {
    p <- p + ggplot2::geom_smooth(
      method = "loess",
      se = FALSE,
      linewidth = 0.8,
      linetype = "dashed",
      color = "grey30"
    )
  }

  # Points Layer
  if (show_points) {
    if (is.null(color_var)) {
      p <- p + ggplot2::geom_point(
        size = point_size,
        shape = 21,
        stroke = 0.6,
        fill = "#2A9D8F",
        color = "black"
      )
    } else {
      p <- p + ggplot2::geom_point(
        ggplot2::aes(fill = .data[[color_var]]),
        size = point_size,
        shape = 21,
        stroke = 0.6,
        color = "black"
      )
    }
  }

  # 4. Faceting & Palette
  if (!is.null(facet_var)) {
    p <- p + ggplot2::facet_wrap(stats::as.formula(paste("~", facet_var)), scales = "free_y")
  }

  if (!is.null(color_var)) {
    p <- p +
      ggplot2::scale_color_brewer(palette = color_palette, name = color_var) +
      ggplot2::scale_fill_brewer(palette = color_palette, name = color_var)
  }

  # 5. Publication Journal Aesthetics
  p <- p +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(
      axis.line = ggplot2::element_line(linewidth = 0.7, color = "black"),
      axis.text = ggplot2::element_text(color = "black", face = "bold"),
      axis.title = ggplot2::element_text(face = "bold"),
      legend.position = if (is.null(color_var) && is.null(linetype_var)) "none" else "top",
      legend.title = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(face = "bold", size = 12, hjust = 0.5),
      plot.subtitle = ggplot2::element_text(size = 9.5, hjust = 0.5, color = "grey35"),
      strip.background = ggplot2::element_rect(fill = "grey92", color = "black"),
      strip.text = ggplot2::element_text(face = "bold")
    ) +
    ggplot2::labs(
      x = if (!is.null(xlab)) xlab else time_var,
      y = if (!is.null(ylab)) ylab else y_var,
      title = if (!is.null(title)) title else paste("Temporal Profile of", y_var),
      subtitle = if (has_replicates && ribbon != "none") paste("Aggregated Mean \u00b1", toupper(ribbon), "Uncertainty Ribbon") else NULL
    )

  return(list(
    Plot = p,
    Summary_Data = summary_df[, c(grp_factors, "N", "Mean", "SD", "SE", "LCI", "UCI")]
  ))
}
