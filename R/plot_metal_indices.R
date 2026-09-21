#' Visualize Heavy Metal Pollution Indices with Confidence Intervals
#'
#' Generates publication-ready plots for heavy metal pollution metrics
#' (HPI, HEI, Cd, PLI, or Nemerow PN) across monitoring stations, featuring
#' 95\% confidence interval whiskers, regulatory benchmark lines, and intuitive status fills.
#'
#' @param data Output list from \code{calc_metal_indices()} or a raw \code{data.frame}.
#' @param metric Character. Index to plot: \code{"HPI"} (default), \code{"HEI"},
#'   \code{"Cd"}, \code{"PLI"}, or \code{"Nemerow_PN"}.
#' @param id_col Character. Station identification column (default: \code{"Station"}).
#' @param type Character. Style: \code{"bar"} (default) or \code{"lollipop"}.
#' @param show_ci Logical. If \code{TRUE}, displays 95\% CI whiskers when replicates exist.
#' @param show_threshold Logical. If \code{TRUE}, draws dashed regulatory thresholds.
#' @param palette Character. Status color theme: \code{"auto"} for standard ecological colors
#'   (green = suitable/clean, red = polluted), or an RColorBrewer palette name.
#' @param title Optional character. Custom plot title.
#' @param ... Additional arguments passed to \code{calc_metal_indices()}.
#'
#' @return A \code{ggplot} object.
#' @export
#'
#' @import ggplot2
#' @importFrom rlang .data
plot_metal_indices <- function(
    data,
    metric = c("HPI", "HEI", "Cd", "PLI", "Nemerow_PN"),
    id_col = "Station",
    type = c("bar", "lollipop"),
    show_ci = TRUE,
    show_threshold = TRUE,
    palette = "auto",
    title = NULL,
    ...
) {

  metric <- match.arg(metric)
  type   <- match.arg(type)

  # Standard Defensive Ingestion Guard
  if (is.matrix(data)) {
    data <- as.data.frame(data)
  }

  # 1. Ingest Data
  if (is.list(data) && "Summary_Data" %in% names(data)) {
    plot_df <- data$Summary_Data
  } else if (is.data.frame(data) && paste0(metric, "_Mean") %in% names(data)) {
    plot_df <- data
  } else {
    calc_res <- calc_metal_indices(data = data, id_col = id_col, ...)
    plot_df <- calc_res$Summary_Data
  }

  plot_df[[id_col]] <- factor(plot_df[[id_col]], levels = unique(plot_df[[id_col]]))

  # 2. Metric Configurations
  meta_map <- list(
    HPI = list(
      mean_col = "HPI_Mean", lci_col = "HPI_LCI", uci_col = "HPI_UCI",
      status_col = "HPI_Status", threshold = 100,
      thresh_label = "Critical Threshold (HPI = 100)",
      ylab = "Heavy Metal Pollution Index (HPI)"
    ),
    HEI = list(
      mean_col = "HEI_Mean", lci_col = "HEI_LCI", uci_col = "HEI_UCI",
      status_col = "HEI_Status", threshold = c(10, 20),
      thresh_label = "Evaluation Limits (10 = Low, 20 = High)",
      ylab = "Heavy Metal Evaluation Index (HEI)"
    ),
    Cd = list(
      mean_col = "Cd_Mean", lci_col = "Cd_LCI", uci_col = "Cd_UCI",
      status_col = "Cd_Status", threshold = c(6, 12),
      thresh_label = "Contamination Limits (6 = Moderate, 12 = High)",
      ylab = "Degree of Contamination (Cd)"
    ),
    PLI = list(
      mean_col = "PLI_Mean", lci_col = "PLI_LCI", uci_col = "PLI_UCI",
      status_col = "PLI_Status", threshold = 1.0,
      thresh_label = "Baseline Limit (PLI = 1.0; > 1 polluted)",
      ylab = "Pollution Load Index (PLI)"
    ),
    Nemerow_PN = list(
      mean_col = "Nemerow_PN_Mean", lci_col = "Nemerow_PN_LCI", uci_col = "Nemerow_PN_UCI",
      status_col = "Nemerow_Status", threshold = c(0.7, 1.0, 2.0, 3.0),
      thresh_label = "Warning (0.7), Slight (1.0), Moderate (2.0), Heavy (3.0)",
      ylab = "Nemerow Pollution Index (PN)"
    )
  )

  cfg <- meta_map[[metric]]
  plot_df$Val <- plot_df[[cfg$mean_col]]
  plot_df$LCI <- pmax(0, plot_df[[cfg$lci_col]])
  plot_df$UCI <- plot_df[[cfg$uci_col]]
  plot_df$Status <- plot_df[[cfg$status_col]]

  has_ci <- show_ci && any(plot_df$N > 1)
  upper_bound <- max(c(plot_df$UCI, cfg$threshold), na.rm = TRUE) * 1.20

  # 3. Base Plot Layout
  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = .data[[id_col]], y = .data$Val)
  )

  if (show_threshold) {
    p <- p + ggplot2::geom_hline(
      yintercept = cfg$threshold,
      linetype = "dashed",
      color = "#D9534F",
      linewidth = 0.8
    )
  }

  # 4. Geometries & Whiskers
  if (type == "bar") {
    p <- p +
      ggplot2::geom_col(
        ggplot2::aes(fill = .data$Status),
        width = 0.6,
        color = "black",
        linewidth = 0.35,
        alpha = 0.88
      )
  } else {
    p <- p +
      ggplot2::geom_segment(
        ggplot2::aes(x = .data[[id_col]], xend = .data[[id_col]], y = 0, yend = .data$Val),
        color = "grey40",
        linewidth = 0.8
      ) +
      ggplot2::geom_point(
        ggplot2::aes(fill = .data$Status),
        size = 4.5,
        shape = 21,
        color = "black"
      )
  }

  if (has_ci) {
    p <- p + ggplot2::geom_errorbar(
      ggplot2::aes(ymin = .data$LCI, ymax = .data$UCI),
      width = 0.22,
      linewidth = 0.7,
      color = "black"
    )
  }

  # Value labels above error bars/bars
  p <- p +
    ggplot2::geom_text(
      ggplot2::aes(
        y = if (has_ci) .data$UCI else .data$Val,
        label = sprintf("%.2f", .data$Val)
      ),
      vjust = -0.7,
      fontface = "bold",
      size = 3.2
    )

  # 5. Ecological Color Mapping
  metal_colors <- c(
    "Suitable"               = "#2A9D8F",
    "Polluted/Critical"      = "#D9534F",
    "Low"                    = "#2A9D8F",
    "Medium"                 = "#E9C46A",
    "High"                   = "#D9534F",
    "Low Contamination"      = "#2A9D8F",
    "Moderate Contamination" = "#E9C46A",
    "High Contamination"     = "#D9534F",
    "Unpolluted"             = "#2A9D8F",
    "Polluted"               = "#D9534F",
    "Clean"                  = "#2A9D8F",
    "Warning Limit"          = "#4575B4",
    "Slightly Polluted"      = "#E9C46A",
    "Moderately Polluted"    = "#F4A261",
    "Heavily Polluted"       = "#D9534F"
  )

  if (identical(palette, "auto")) {
    p <- p + ggplot2::scale_fill_manual(values = metal_colors, drop = FALSE)
  } else if (is.character(palette) && length(palette) == 1 && palette %in% c("Set1", "Set2", "Set3", "Dark2", "Paired", "Accent", "Spectral")) {
    p <- p + ggplot2::scale_fill_brewer(palette = palette)
  } else {
    p <- p + ggplot2::scale_fill_manual(values = metal_colors, drop = FALSE)
  }

  # 6. Theme and Layout
  need_tilt <- length(unique(plot_df[[id_col]])) > 3 || any(nchar(as.character(plot_df[[id_col]])) > 8)

  p <- p +
    ggplot2::scale_y_continuous(limits = c(0, upper_bound), expand = c(0, 0)) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(
      axis.line = ggplot2::element_line(linewidth = 0.7, color = "black"),
      axis.title = ggplot2::element_text(face = "bold"),
      axis.text = ggplot2::element_text(color = "black", face = "bold"),
      axis.text.x = ggplot2::element_text(
        color = "black",
        face = "bold",
        angle = if (need_tilt) 35 else 0,
        hjust = if (need_tilt) 1 else 0.5
      ),
      legend.position = "top",
      legend.title = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(face = "bold", size = 11, hjust = 0.5),
      plot.subtitle = ggplot2::element_text(size = 9.5, hjust = 0.5, color = "grey35")
    ) +
    ggplot2::labs(
      x = "Monitoring Station",
      y = cfg$ylab,
      fill = "Status",
      title = if (!is.null(title)) title else paste(cfg$ylab, "(Mean \u00b1 95% CI)"),
      subtitle = if (show_threshold) cfg$thresh_label else NULL
    )

  return(p)
}
