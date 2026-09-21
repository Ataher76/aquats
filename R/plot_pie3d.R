#' Publication-Ready 3D Depth Pie Chart
#'
#' Generates an isometric 3D pie chart with perspective tilt, extruded
#' depth ribbons, and realistic shading natively in ggplot2.
#'
#' @param data A \code{data.frame}, \code{matrix}, or file path to a \code{.csv} or \code{.xlsx} file.
#' @param x Character. Categorical variable for the pie slices.
#' @param y Optional character. Numeric variable to aggregate. If \code{NULL}, calculates counts.
#' @param stat Character. Aggregation method: \code{"sum"} (default), \code{"mean"}, or \code{"identity"}.
#' @param tilt Numeric. Perspective vertical compression factor between 0.2 and 0.8 (default: 0.5).
#' @param depth Numeric. Extrusion thickness of the 3D rim (default: 0.25).
#' @param show_labels Logical. If \code{TRUE}, displays percentage labels on slices.
#' @param min_percent Numeric. Minimum percentage required to display slice text (default: 4).
#' @param palette Character. RColorBrewer palette name (default: \code{"Set2"}).
#' @param title Optional character. Plot title.
#'
#' @return A \code{ggplot} object.
#' @export
#'
#' @import ggplot2
#' @importFrom stats aggregate na.omit
#' @importFrom tools file_ext
#' @importFrom utils read.csv
#' @importFrom rlang .data
plot_pie3d <- function(data,
                       x,
                       y = NULL,
                       stat = c("sum", "mean", "identity"),
                       tilt = 0.5,
                       depth = 0.25,
                       show_labels = TRUE,
                       min_percent = 4,
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

  if (!x %in% names(df)) stop(sprintf("Column '%s' not found.", x))
  if (!is.null(y) && !y %in% names(df)) stop(sprintf("Column '%s' not found.", y))

  cols <- c(x, y)
  clean_df <- stats::na.omit(df[, cols, drop = FALSE])
  clean_df[[x]] <- as.factor(clean_df[[x]])

  # 2. Compute Proportions
  if (is.null(y)) {
    agg <- as.data.frame(table(clean_df[[x]]))
    colnames(agg) <- c("Category", "Value")
  } else {
    if (!is.numeric(clean_df[[y]])) stop(sprintf("Column '%s' must be numeric.", y))
    if (stat == "identity") {
      agg <- clean_df
      agg$Category <- agg[[x]]
      agg$Value <- agg[[y]]
    } else {
      fn <- if (stat == "mean") mean else sum
      agg <- stats::aggregate(clean_df[[y]], by = list(clean_df[[x]]), FUN = fn)
      colnames(agg) <- c("Category", "Value")
    }
  }

  agg <- agg[agg$Value > 0, , drop = FALSE]
  if (nrow(agg) == 0) stop("No positive values available to plot.")

  grand_total <- sum(agg$Value)
  agg$prop <- agg$Value / grand_total
  agg$pct <- agg$prop * 100

  # 3. Angular Bounds for Each Slice
  agg$end_angle <- cumsum(agg$prop) * 2 * pi
  agg$start_angle <- c(0, agg$end_angle[-nrow(agg)])
  agg$mid_angle <- (agg$start_angle + agg$end_angle) / 2

  # 4. Generate 3D Polygons (Top Slices & Front Rims)
  radius <- 1.0
  top_poly_list <- list()
  rim_poly_list <- list()

  for (i in seq_len(nrow(agg))) {
    a1 <- agg$start_angle[i]
    a2 <- agg$end_angle[i]
    cat_name <- as.character(agg$Category[i])

    # Top Slice Polygon
    n_pts <- max(12, ceiling(100 * agg$prop[i]))
    theta_top <- seq(a1, a2, length.out = n_pts)
    x_top <- c(0, radius * cos(theta_top), 0)
    y_top <- c(0, radius * sin(theta_top) * tilt, 0)

    top_poly_list[[i]] <- data.frame(
      x = x_top,
      y = y_top,
      Category = cat_name,
      piece = paste0("top_", i),
      stringsAsFactors = FALSE
    )

    # Front Rim Ribbon (overlaps interval pi to 2*pi where y <= 0)
    rim_a1 <- max(a1, pi)
    rim_a2 <- min(a2, 2 * pi)

    if (rim_a1 < rim_a2) {
      n_rim <- max(8, ceiling(80 * (rim_a2 - rim_a1) / (2 * pi)))
      theta_rim <- seq(rim_a1, rim_a2, length.out = n_rim)
      x_rim_top <- radius * cos(theta_rim)
      y_rim_top <- radius * sin(theta_rim) * tilt
      x_rim_bot <- radius * cos(rev(theta_rim))
      y_rim_bot <- (radius * sin(rev(theta_rim)) * tilt) - depth

      rim_poly_list[[length(rim_poly_list) + 1]] <- data.frame(
        x = c(x_rim_top, x_rim_bot),
        y = c(y_rim_top, y_rim_bot),
        Category = cat_name,
        piece = paste0("rim_", i),
        stringsAsFactors = FALSE
      )
    }
  }

  top_df <- do.call(rbind, top_poly_list)
  rim_df <- if (length(rim_poly_list) > 0) do.call(rbind, rim_poly_list) else NULL

  # Slice Labels
  agg$label_x <- 0.65 * radius * cos(agg$mid_angle)
  agg$label_y <- 0.65 * radius * sin(agg$mid_angle) * tilt
  agg$label_text <- paste0(round(agg$pct), "%")

  # 5. Build ggplot Object
  p <- ggplot2::ggplot()

  # Draw 3D rim base
  if (!is.null(rim_df)) {
    p <- p +
      ggplot2::geom_polygon(
        data = rim_df,
        ggplot2::aes(
          x = .data$x,
          y = .data$y,
          group = .data$piece,
          fill = .data$Category
        ),
        color = "black",
        linewidth = 0.35,
        show.legend = FALSE
      ) +
      ggplot2::geom_polygon(
        data = rim_df,
        ggplot2::aes(
          x = .data$x,
          y = .data$y,
          group = .data$piece
        ),
        fill = "black",
        alpha = 0.28,
        color = "black",
        linewidth = 0.35,
        show.legend = FALSE
      )
  }

  # Draw Top Faces
  p <- p +
    ggplot2::geom_polygon(
      data = top_df,
      ggplot2::aes(
        x = .data$x,
        y = .data$y,
        group = .data$piece,
        fill = .data$Category
      ),
      color = "black",
      linewidth = 0.4
    )

  # Labels on Slices
  if (show_labels) {
    lbl_data <- agg[agg$pct >= min_percent, ]
    p <- p + ggplot2::geom_text(
      data = lbl_data,
      ggplot2::aes(
        x = .data$label_x,
        y = .data$label_y,
        label = .data$label_text
      ),
      size = 3.5,
      fontface = "bold",
      color = "black"
    )
  }

  p <- p +
    ggplot2::coord_fixed() +
    ggplot2::scale_fill_brewer(palette = palette) +
    ggplot2::theme_void(base_size = 12) +
    ggplot2::theme(
      legend.position = "right",
      legend.title = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, size = 13)
    )

  if (!is.null(title)) {
    p <- p + ggplot2::labs(title = title)
  }

  return(p)
}
