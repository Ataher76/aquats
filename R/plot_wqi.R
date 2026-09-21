#' Visualize Water Quality Index (WQI) with Confidence Intervals
#'
#' Generates publication-ready WQI plots across monitoring stations with
#' 95\% confidence interval whiskers, status fills, and background
#' benchmark bands for WAWQI classification tiers.
#'
#' @param data Output list from \code{calc_wqi()} or a raw \code{data.frame}.
#' @param id_col Character. Column identifying stations (default: \code{"Station"}).
#' @param type Character. Plot style: \code{"bar"} (default) or \code{"lollipop"}.
#' @param show_ci Logical. If \code{TRUE}, displays 95\% CI whiskers when replicates exist.
#' @param show_bands Logical. If \code{TRUE}, adds shaded background benchmark tiers.
#' @param title Optional character. Custom plot title.
#' @param ... Additional arguments passed to \code{calc_wqi()}.
#'
#' @return A \code{ggplot} object.
#' @export
#'
#' @import ggplot2
#' @importFrom rlang .data
plot_wqi <- function(
    data,
    id_col = "Station",
    type = c("bar", "lollipop"),
    show_ci = TRUE,
    show_bands = TRUE,
    title = NULL,
    ...
) {

  type <- match.arg(type)

  # Standard Defensive Ingestion Guard
  if (is.matrix(data)) {
    data <- as.data.frame(data)
  }

  # 1. Ingest Data
  if (is.list(data) && "Summary_Data" %in% names(data)) {
    plot_df <- data$Summary_Data
  } else if (is.data.frame(data) && "WQI_Mean" %in% names(data)) {
    plot_df <- data
  } else {
    calc_res <- calc_wqi(data = data, id_col = id_col, ...)
    plot_df <- calc_res$Summary_Data
  }

  plot_df[[id_col]] <- factor(plot_df[[id_col]], levels = unique(plot_df[[id_col]]))

  plot_df$Val <- plot_df$WQI_Mean
  plot_df$LCI <- pmax(0, plot_df$WQI_LCI)
  plot_df$UCI <- plot_df$WQI_UCI

  has_ci <- show_ci && any(plot_df$N > 1)
  max_val <- if (has_ci) max(plot_df$UCI, na.rm = TRUE) else max(plot_df$Val, na.rm = TRUE)
  upper_bound <- max(c(110, max_val * 1.15))

  # 2. Color Palette for Ecological Ratings
  wqi_colors <- c(
    "Excellent"  = "#2A9D8F",
    "Good"       = "#4575B4",
    "Poor"       = "#E9C46A",
    "Very Poor"  = "#F4A261",
    "Unsuitable" = "#D9534F"
  )

  # 3. Base Plot Layout
  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = .data[[id_col]], y = .data$Val)
  )

  # 4. Optional Shaded Benchmark Bands
  if (show_bands) {
    p <- p +
      ggplot2::annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0,   ymax = 25,  fill = "#2A9D8F", alpha = 0.08) +
      ggplot2::annotate("rect", xmin = -Inf, xmax = Inf, ymin = 25,  ymax = 50,  fill = "#4575B4", alpha = 0.08) +
      ggplot2::annotate("rect", xmin = -Inf, xmax = Inf, ymin = 50,  ymax = 75,  fill = "#E9C46A", alpha = 0.10) +
      ggplot2::annotate("rect", xmin = -Inf, xmax = Inf, ymin = 75,  ymax = 100, fill = "#F4A261", alpha = 0.12) +
      ggplot2::annotate("rect", xmin = -Inf, xmax = Inf, ymin = 100, ymax = upper_bound, fill = "#D9534F", alpha = 0.12) +
      ggplot2::geom_hline(yintercept = c(25, 50, 75, 100), linetype = "dashed", color = "grey65", linewidth = 0.4)
  }

  # 5. Geometries & Whiskers
  if (type == "bar") {
    p <- p +
      ggplot2::geom_col(
        ggplot2::aes(fill = .data$Status),
        width = 0.6,
        color = "black",
        linewidth = 0.35,
        alpha = 0.9
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

  # Value labels
  p <- p +
    ggplot2::geom_text(
      ggplot2::aes(
        y = if (has_ci) .data$UCI else .data$Val,
        label = sprintf("%.1f", .data$Val)
      ),
      vjust = -0.7,
      fontface = "bold",
      size = 3.2
    )

  # 6. Theme and Layout
  p <- p +
    ggplot2::scale_y_continuous(limits = c(0, upper_bound), expand = c(0, 0)) +
    ggplot2::scale_fill_manual(values = wqi_colors, drop = FALSE) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(
      axis.line = ggplot2::element_line(linewidth = 0.7, color = "black"),
      axis.title = ggplot2::element_text(face = "bold"),
      axis.text = ggplot2::element_text(color = "black", face = "bold"),
      axis.text.x = ggplot2::element_text(
        angle = if (length(unique(plot_df[[id_col]])) > 3) 35 else 0,
        hjust = if (length(unique(plot_df[[id_col]])) > 3) 1 else 0.5
      ),
      legend.position = "top",
      legend.title = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(face = "bold", size = 12, hjust = 0.5)
    ) +
    ggplot2::labs(
      x = "Monitoring Station",
      y = "Water Quality Index (WQI)",
      fill = "Status",
      title = if (!is.null(title)) title else "Water Quality Index (WAWQI) (Mean \u00b1 95% CI)"
    )

  return(p)
}
