#' Publication-Ready Sankey and Alluvial Flow Diagram
#'
#' Generates publication-ready multi-tier alluvial and Sankey flow diagrams
#' to visualize categorical transitions, trophic pathways, diet shifts, or
#' spatial routing in ecological and environmental data. Automatically pre-aggregates
#' frequencies, guards against missing strata, dynamically adapts color palettes,
#' and returns a fully customizable \code{ggplot2} object.
#'
#' @param data A \code{data.frame}, \code{matrix}, or character string file path
#'   to a \code{.csv} or \code{.xlsx} file containing the observations.
#' @param axes Character vector; column names in \code{data} specifying two or more
#'   categorical variables to serve as consecutive vertical axes (strata) from left to right.
#' @param value Optional character string; column name in \code{data} representing
#'   the flow magnitude, frequency, or biomass. If \code{NULL}, observation frequencies
#'   (row counts) are calculated automatically. Default is \code{NULL}.
#' @param fill_axis Character string or integer; column name or index in \code{axes}
#'   used to color the alluvial ribbons. Default is \code{1}.
#' @param alpha Numeric; transparency of the alluvial flow ribbons (between 0 and 1).
#'   Default is \code{0.65}.
#' @param width Numeric; relative horizontal width of stratum bars (between 0 and 1).
#'   Default is \code{0.25}.
#' @param stratum_fill Character string; fill color of the stratum bounding boxes.
#'   Default is \code{"grey92"}.
#' @param stratum_color Character string; border color of the stratum bounding boxes.
#'   Default is \code{"grey40"}.
#' @param label_size Numeric; font size for stratum category labels. Default is \code{3.5}.
#' @param palette Character string; a valid \code{RColorBrewer} palette name. Default is \code{"Set2"}.
#' @param title Optional character string; custom plot title. Default is \code{NULL}.
#'
#' @return A publication-ready \code{ggplot2} alluvial diagram object.
#' @export
#'
#' @import ggplot2
#' @importFrom ggalluvial geom_alluvium geom_stratum StatStratum
#' @importFrom stats aggregate na.omit
#' @importFrom tools file_ext
#' @importFrom utils read.csv
#' @importFrom grDevices colorRampPalette
#' @importFrom RColorBrewer brewer.pal
#' @importFrom rlang sym .data
plot_sankey <- function(
    data,
    axes,
    value = NULL,
    fill_axis = 1,
    alpha = 0.65,
    width = 0.25,
    stratum_fill = "grey92",
    stratum_color = "grey40",
    label_size = 3.5,
    palette = "Set2",
    title = NULL
) {

  # 1. Defensive Ingestion Guard
  if (is.character(data) && length(data) == 1) {
    if (!file.exists(data)) {
      stop(sprintf("File does not exist: '%s'", data))
    }
    ext <- tolower(tools::file_ext(data))
    if (ext == "csv") {
      data <- utils::read.csv(data, stringsAsFactors = FALSE)
    } else if (ext %in% c("xlsx", "xls")) {
      if (!requireNamespace("readxl", quietly = TRUE)) {
        stop("Package 'readxl' is required to read Excel files. Install it or supply a data.frame.")
      }
      data <- as.data.frame(readxl::read_excel(data))
    } else {
      stop("Unsupported file format. Provide a .csv, .xlsx file, or a data.frame.")
    }
  } else if (is.matrix(data) || inherits(data, "tbl_df") || inherits(data, "tbl")) {
    data <- as.data.frame(data)
  } else if (!is.data.frame(data)) {
    stop("Input 'data' must be a data.frame, tibble, matrix, or file path.")
  }

  data <- as.data.frame(data)

  # 2. Argument and Axis Validation
  if (!is.character(axes) || length(axes) < 2) {
    stop("'axes' must be a character vector specifying at least 2 column names.")
  }

  target_cols <- axes
  if (!is.null(value)) {
    if (!is.character(value) || length(value) != 1) {
      stop("'value' must be a single character string representing a numeric column.")
    }
    target_cols <- c(target_cols, value)
  }

  missing_cols <- setdiff(target_cols, names(data))
  if (length(missing_cols) > 0) {
    stop(sprintf("Column(s) not found in 'data': %s", paste(missing_cols, collapse = ", ")))
  }

  # Validate fill axis
  if (is.numeric(fill_axis)) {
    if (fill_axis < 1 || fill_axis > length(axes)) {
      stop(sprintf("'fill_axis' index (%d) is out of bounds (must be between 1 and %d).", fill_axis, length(axes)))
    }
    fill_var <- axes[as.integer(fill_axis)]
  } else if (is.character(fill_axis)) {
    if (!fill_axis %in% axes) {
      stop(sprintf("'fill_axis' ('%s') must match one of the specified 'axes': %s", fill_axis, paste(axes, collapse = ", ")))
    }
    fill_var <- fill_axis
  } else {
    stop("'fill_axis' must be a column name or an integer index referencing 'axes'.")
  }

  # 3. Clean Missing Cases and Validate Numerical Flow
  sub_df <- data[, target_cols, drop = FALSE]
  n_initial <- nrow(sub_df)
  sub_df <- stats::na.omit(sub_df)
  n_dropped <- n_initial - nrow(sub_df)

  if (n_dropped > 0) {
    warning(sprintf("Removed %d observation(s) containing missing values across selected axes.", n_dropped), call. = FALSE)
  }

  if (nrow(sub_df) == 0) {
    stop("No valid observations remaining after removing missing values.")
  }

  if (!is.null(value)) {
    if (!is.numeric(sub_df[[value]])) {
      stop(sprintf("Flow magnitude variable '%s' must be numeric.", value))
    }
    if (any(sub_df[[value]] < 0)) {
      stop(sprintf("Flow magnitude variable '%s' contains negative values. All flows must be >= 0.", value))
    }
  }

  # Convert all axes to factors to preserve strata integrity
  for (ax in axes) {
    sub_df[[ax]] <- as.factor(sub_df[[ax]])
  }

  # 4. Aggregate Flow Frequencies
  if (is.null(value)) {
    agg_df <- stats::aggregate(
      list(Weight = rep(1, nrow(sub_df))),
      by = sub_df[, axes, drop = FALSE],
      FUN = length
    )
    val_col <- "Weight"
  } else {
    agg_df <- stats::aggregate(
      sub_df[[value]],
      by = sub_df[, axes, drop = FALSE],
      FUN = sum
    )
    names(agg_df)[ncol(agg_df)] <- "Weight"
    val_col <- "Weight"
  }

  # Exclude zero-weight flows
  agg_df <- agg_df[agg_df[[val_col]] > 0, , drop = FALSE]
  if (nrow(agg_df) == 0) {
    stop("All flow aggregates are zero. Unable to construct diagram.")
  }

  # 5. Build Dynamic Mapping for ggalluvial
  aes_args <- list()
  aes_args[["y"]] <- rlang::sym(val_col)
  for (i in seq_along(axes)) {
    aes_args[[paste0("axis", i)]] <- rlang::sym(axes[i])
  }
  main_aes <- do.call(ggplot2::aes, aes_args)

  # 6. Construct Publication Graphic
  p <- ggplot2::ggplot(agg_df, main_aes) +
    ggalluvial::geom_alluvium(
      ggplot2::aes(fill = .data[[fill_var]]),
      width = width,
      alpha = alpha,
      curve_type = "sine",
      knot.pos = 0.4
    ) +
    ggalluvial::geom_stratum(
      width = width,
      fill = stratum_fill,
      color = stratum_color,
      linewidth = 0.45
    ) +
    ggplot2::geom_text(
      stat = ggalluvial::StatStratum,
      ggplot2::aes(label = ggplot2::after_stat(stratum)),
      size = label_size,
      fontface = "bold",
      color = "black"
    ) +
    ggplot2::scale_x_discrete(limits = axes, expand = c(0.08, 0.08)) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.02, 0.06)))

  # 7. Adaptive Palette Handling
  fill_levels <- unique(agg_df[[fill_var]])
  n_levels <- length(fill_levels)

  pal_colors <- tryCatch({
    RColorBrewer::brewer.pal(max(3, n_levels), palette)
  }, error = function(e) NULL)

  if (!is.null(pal_colors)) {
    if (n_levels > length(pal_colors)) {
      extended_pal <- grDevices::colorRampPalette(pal_colors)(n_levels)
      p <- p + ggplot2::scale_fill_manual(values = extended_pal, name = fill_var)
    } else {
      p <- p + ggplot2::scale_fill_brewer(palette = palette, name = fill_var)
    }
  } else {
    p <- p + ggplot2::scale_fill_discrete(name = fill_var)
  }

  # 8. Publication Theming
  default_title <- sprintf("Alluvial Transitions: %s", paste(axes, collapse = " \u2192 "))
  subtitle_text <- sprintf(
    "Stratified flows weighted by %s (Ribbons colored by '%s')",
    if (is.null(value)) "observation counts" else value,
    fill_var
  )

  p <- p +
    ggplot2::labs(
      title = if (!is.null(title)) title else default_title,
      subtitle = subtitle_text,
      x = "Stratum Axis",
      y = if (is.null(value)) "Cumulative Count" else sprintf("Total %s", value)
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(face = "bold", size = 11, color = "black"),
      axis.text.y = ggplot2::element_text(color = "grey30"),
      axis.title = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(face = "bold", size = 13, hjust = 0.5),
      plot.subtitle = ggplot2::element_text(color = "grey35", size = 10, hjust = 0.5),
      legend.position = "right",
      legend.title = ggplot2::element_text(face = "bold")
    )

  return(p)
}

