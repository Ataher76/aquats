#' Publication-Ready Single, Double, and Triple-Layered Donut Plots
#'
#' Generates publication-ready single or nested multi-layered (sunburst)
#' donut charts. Supports up to three concentric hierarchical tiers
#' with automatic alignment and percentage threshold labelling.
#'
#' @param data A \code{data.frame}, \code{matrix}, or file path to a \code{.csv} or \code{.xlsx} file.
#' @param levels Character vector of 1 to 3 categorical variable names
#'   ordered from innermost ring to outermost ring.
#' @param value Optional character. Numeric variable to aggregate (sum).
#'   If \code{NULL}, counts frequencies.
#' @param show_labels Logical. If \code{TRUE}, adds percentage labels to segments.
#' @param min_percent Numeric. Minimum percentage (0 to 100) required to display
#'   a text label on a segment (prevents visual clutter; default: 4).
#' @param palette Character. RColorBrewer palette name (default: \code{"Set2"}).
#' @param title Optional character string for plot title.
#'
#' @return A \code{ggplot} object.
#' @export
#'
#' @import ggplot2
#' @importFrom stats aggregate na.omit
#' @importFrom tools file_ext
#' @importFrom utils read.csv
#' @importFrom rlang .data
plot_donut <- function(data,
                       levels,
                       value = NULL,
                       show_labels = TRUE,
                       min_percent = 4,
                       palette = "Set2",
                       title = NULL) {

  # 1. Standard Defensive Ingestion Guard
  if (is.character(data) && length(data) == 1) {
    if (!file.exists(data)) stop(sprintf("File '%s' not found.", data))
    ext <- tolower(tools::file_ext(data))
    if (ext == "csv") {
      df <- utils::read.csv(data, stringsAsFactors = FALSE)
    } else if (ext %in% c("xlsx", "xls")) {
      if (!requireNamespace("readxl", quietly = TRUE)) {
        stop("Please install 'readxl' to load Excel files.")
      }
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

  n_tiers <- length(levels)
  if (n_tiers < 1 || n_tiers > 3) {
    stop("'levels' must contain 1, 2, or 3 categorical variable names.")
  }

  missing_cols <- levels[!levels %in% names(df)]
  if (length(missing_cols) > 0) {
    stop(sprintf("Columns not found: %s", paste(missing_cols, collapse = ", ")))
  }
  if (!is.null(value) && !value %in% names(df)) {
    stop(sprintf("Value column '%s' not found.", value))
  }

  cols_needed <- c(levels, value)
  clean_df <- stats::na.omit(df[, cols_needed, drop = FALSE])

  for (lvl in levels) {
    clean_df[[lvl]] <- as.factor(clean_df[[lvl]])
  }

  # 2. Define Concentric Ring Radii
  radii <- list(
    tier1 = list(r_in = 1.8, r_out = 2.8),
    tier2 = list(r_in = c(1.5, 2.35), r_out = c(2.25, 3.1)),
    tier3 = list(r_in = c(1.3, 1.95, 2.6), r_out = c(1.85, 2.5, 3.15))
  )
  cfg <- radii[[paste0("tier", n_tiers)]]

  # Compute Global Total
  grand_total <- if (is.null(value)) nrow(clean_df) else sum(clean_df[[value]], na.rm = TRUE)

  # 3. Calculate Geometry for Each Ring
  ring_frames <- list()

  for (k in seq_len(n_tiers)) {
    current_grp <- levels[1:k]

    if (is.null(value)) {
      agg <- stats::aggregate(
        clean_df[[current_grp[1]]],
        by = clean_df[, current_grp, drop = FALSE],
        FUN = length
      )
      colnames(agg)[ncol(agg)] <- "weight"
    } else {
      agg <- stats::aggregate(
        clean_df[[value]],
        by = clean_df[, current_grp, drop = FALSE],
        FUN = sum
      )
      colnames(agg)[ncol(agg)] <- "weight"
    }

    sort_order <- do.call(order, agg[, current_grp, drop = FALSE])
    agg <- agg[sort_order, , drop = FALSE]

    # Proportions & Angular Coordinates
    agg$pct <- (agg$weight / grand_total) * 100
    agg$prop <- agg$weight / grand_total
    agg$ymax <- cumsum(agg$prop)
    agg$ymin <- agg$ymax - agg$prop
    agg$ymid <- (agg$ymin + agg$ymax) / 2

    # Radial Coordinates
    agg$xmin <- cfg$r_in[k]
    agg$xmax <- cfg$r_out[k]
    agg$xmid <- (agg$xmin + agg$xmax) / 2
    agg$ring_level <- k
    agg$label_name <- as.character(agg[[levels[k]]])

    ring_frames[[k]] <- agg
  }

  common_cols <- c("xmin", "xmax", "xmid", "ymin", "ymax", "ymid", "pct", "ring_level", "label_name")
  standardized <- lapply(ring_frames, function(d) d[, common_cols])
  plot_data <- do.call(rbind, standardized)
  plot_data$label_name <- factor(plot_data$label_name, levels = unique(plot_data$label_name))

  # 4. Construct Concentric Donut Plot
  p <- ggplot2::ggplot(plot_data) +
    ggplot2::geom_rect(
      ggplot2::aes(
        xmin = .data$xmin,
        xmax = .data$xmax,
        ymin = .data$ymin,
        ymax = .data$ymax,
        fill = .data$label_name
      ),
      color = "white",
      linewidth = 0.5,
      alpha = 0.9
    ) +
    ggplot2::coord_polar(theta = "y", start = 0) +
    ggplot2::xlim(c(0, max(cfg$r_out) + 0.3)) +
    ggplot2::ylim(c(0, 1)) +
    ggplot2::scale_fill_brewer(palette = palette) +
    ggplot2::theme_void(base_size = 12) +
    ggplot2::theme(
      legend.position = "right",
      legend.title = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, size = 13)
    )

  # 5. Segment Labels
  if (show_labels) {
    lbl_data <- plot_data[plot_data$pct >= min_percent, ]
    p <- p + ggplot2::geom_text(
      data = lbl_data,
      ggplot2::aes(
        x = .data$xmid,
        y = .data$ymid,
        label = paste0(round(.data$pct), "%")
      ),
      size = 3.2,
      fontface = "bold",
      color = "black",
      inherit.aes = FALSE
    )
  }

  if (!is.null(title)) {
    p <- p + ggplot2::labs(title = title)
  }

  return(p)
}
