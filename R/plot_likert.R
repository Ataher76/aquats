#' Publication-Ready Likert Scale Visualization
#'
#' Generates modern diverging or 100% stacked horizontal bar charts
#' for survey items and perception data.
#'
#' @param data A \code{data.frame}, \code{matrix}, or file path to a .csv or .xlsx file.
#' @param items Character vector. Column names of Likert question items.
#' @param levels Optional character vector of ordered response levels.
#' @param type Character. Plot style: "diverging" (default) or "stacked".
#' @param palette Character. RColorBrewer palette name (default: "RdYlBu").
#' @param show_labels Logical. If TRUE, displays percentage numbers on segments.
#' @param clean_labels Logical. If TRUE, replaces underscores with spaces in item names.
#' @param title Optional character. Plot title.
#'
#' @return A ggplot object.
#' @export
#'
#' @import ggplot2
#' @importFrom stats na.omit
#' @importFrom tools file_ext
#' @importFrom utils read.csv
#' @importFrom rlang .data
plot_likert <- function(data,
                        items,
                        levels = NULL,
                        type = c("diverging", "stacked"),
                        palette = "RdYlBu",
                        show_labels = TRUE,
                        clean_labels = TRUE,
                        title = NULL) {

  type <- match.arg(type)

  # 1. Standard Defensive Ingestion Guard
  if (is.character(data) && length(data) == 1) {
    if (!file.exists(data)) stop(sprintf("File '%s' not found.", data))
    ext <- tolower(tools::file_ext(data))
    if (ext == "csv") {
      df <- utils::read.csv(data, stringsAsFactors = FALSE, check.names = FALSE)
    } else if (ext %in% c("xlsx", "xls")) {
      if (!requireNamespace("readxl", quietly = TRUE)) stop("Install 'readxl' to read Excel files.")
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

  missing_cols <- items[!items %in% names(df)]
  if (length(missing_cols) > 0) {
    stop(sprintf("Columns not found: %s", paste(missing_cols, collapse = ", ")))
  }

  display_names <- if (clean_labels) gsub("_", " ", items) else items
  names(display_names) <- items

  # 2. Ordered Levels
  all_vals <- stats::na.omit(unlist(df[, items, drop = FALSE]))
  lvl_order <- if (is.null(levels)) sort(unique(as.character(all_vals))) else levels
  n_levels <- length(lvl_order)

  # 3. Compute Percentages
  res_list <- list()
  for (it in items) {
    vals <- stats::na.omit(df[[it]])
    tbl <- table(factor(vals, levels = lvl_order))
    pcts <- (tbl / length(vals)) * 100
    res_list[[it]] <- data.frame(
      Item = display_names[[it]],
      Response = names(pcts),
      Percent = as.numeric(pcts),
      stringsAsFactors = FALSE
    )
  }
  calc_df <- do.call(rbind, res_list)
  calc_df$Response <- factor(calc_df$Response, levels = lvl_order)
  calc_df$Item <- factor(calc_df$Item, levels = rev(display_names))

  # 4. Layout
  if (type == "diverging") {
    mid_idx <- ceiling(n_levels / 2)
    has_neutral <- (n_levels %% 2 == 1)

    plot_rows <- list()
    for (it_display in display_names) {
      sub <- calc_df[calc_df$Item == it_display, ]
      pcts <- sub$Percent

      left_sum <- if (mid_idx > 1) sum(pcts[1:(mid_idx - 1)]) else 0
      neutral_pct <- if (has_neutral) pcts[mid_idx] else 0
      start_x <- -(left_sum + (neutral_pct / 2))

      cur_x <- start_x
      for (i in seq_len(n_levels)) {
        w <- pcts[i]
        plot_rows[[length(plot_rows) + 1]] <- data.frame(
          Item = it_display,
          Response = lvl_order[i],
          Percent = w,
          xmin = cur_x,
          xmax = cur_x + w,
          xmid = cur_x + (w / 2),
          stringsAsFactors = FALSE
        )
        cur_x <- cur_x + w
      }
    }
    div_df <- do.call(rbind, plot_rows)
    div_df$Response <- factor(div_df$Response, levels = lvl_order)
    div_df$Item <- factor(div_df$Item, levels = rev(display_names))

    p <- ggplot2::ggplot(div_df) +
      ggplot2::geom_vline(xintercept = 0, color = "grey50", linewidth = 0.6, linetype = "dashed") +
      ggplot2::geom_rect(
        ggplot2::aes(
          ymin = as.numeric(.data$Item) - 0.32,
          ymax = as.numeric(.data$Item) + 0.32,
          xmin = .data$xmin,
          xmax = .data$xmax,
          fill = .data$Response
        ),
        color = "white",
        linewidth = 0.4
      ) +
      ggplot2::scale_x_continuous(
        labels = function(x) paste0(abs(x), "%"),
        expand = ggplot2::expansion(mult = c(0.08, 0.08))
      ) +
      ggplot2::scale_y_continuous(
        breaks = seq_along(levels(div_df$Item)),
        labels = levels(div_df$Item)
      ) +
      ggplot2::labs(x = "Percentage (%)", y = NULL, fill = NULL, title = title)

    if (show_labels) {
      lbl_df <- div_df[div_df$Percent >= 6, ]
      p <- p + ggplot2::geom_text(
        data = lbl_df,
        ggplot2::aes(
          x = .data$xmid,
          y = as.numeric(.data$Item),
          label = paste0(round(.data$Percent), "%")
        ),
        size = 3.2,
        fontface = "bold",
        color = "black"
      )
    }

  } else {
    p <- ggplot2::ggplot(
      calc_df,
      ggplot2::aes(x = .data$Percent, y = .data$Item, fill = .data$Response)
    ) +
      ggplot2::geom_col(position = "fill", width = 0.65, color = "white", linewidth = 0.4) +
      ggplot2::scale_x_continuous(
        labels = function(x) paste0(round(x * 100), "%"),
        expand = c(0, 0)
      ) +
      ggplot2::labs(x = "Proportion", y = NULL, fill = NULL, title = title)

    if (show_labels) {
      p <- p + ggplot2::geom_text(
        ggplot2::aes(label = paste0(round(.data$Percent), "%")),
        position = ggplot2::position_fill(vjust = 0.5),
        size = 3.2,
        fontface = "bold",
        color = "black"
      )
    }
  }

  # 5. Styling
  p <- p +
    ggplot2::scale_fill_brewer(palette = palette) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(
      axis.line.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(face = "bold", color = "black", size = 11),
      axis.text.x = ggplot2::element_text(color = "black"),
      axis.title.x = ggplot2::element_text(face = "bold"),
      legend.position = "top",
      legend.justification = "center",
      legend.text = ggplot2::element_text(size = 9.5),
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5)
    ) +
    ggplot2::guides(fill = ggplot2::guide_legend(nrow = 1, byrow = TRUE))

  return(p)
}
