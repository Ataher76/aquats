#' Flexible Correlation Matrix and Heatmap Generator for Environmental Data
#'
#' Computes correlation matrices (Pearson, Spearman, or Kendall) with exact p-values
#' and generates publication-ready heatmaps. Supports customizable shapes (squares, circles),
#' triangle views (lower, upper, full), optional diagonal display, significance asterisks,
#' and automated contrast text adjustment.
#'
#' @param data A \code{data.frame} or \code{matrix} containing numeric variables.
#' @param numeric_vars Optional character vector specifying numeric columns to include.
#'   If \code{NULL}, all numeric columns are automatically detected.
#' @param method Character. Correlation method: \code{"pearson"} (default),
#'   \code{"spearman"}, or \code{"kendall"}.
#' @param shape Character. Heatmap marker shape: \code{"square"} (default) or \code{"circle"}.
#' @param view Character. Matrix layout view: \code{"lower"} (default), \code{"upper"}, or \code{"full"}.
#' @param diag Logical. Whether to show diagonal self-correlations (default: \code{FALSE}).
#' @param label_type Character. Content displayed inside cells: \code{"both"} (number + stars),
#'   \code{"stars"} (only asterisks, recommended for circles), \code{"values"} (only numbers),
#'   or \code{"none"}. Defaults to \code{"both"} for squares and \code{"stars"} for circles.
#' @param sig_level Optional numeric (e.g., 0.05). If specified, non-significant correlations
#'   are rendered transparent or blanked out.
#' @param color_palette Character. Diverging palette name for \code{scale_fill_distiller}
#'   (default: \code{"RdBu"}). Alternatives: \code{"RdYlBu"}, \code{"Spectral"}, \code{"BrBG"}.
#' @param title Optional character. Plot title.
#'
#' @return A list containing:
#'   \item{Correlation_Matrix}{Symmetric matrix of correlation coefficients.}
#'   \item{P_Value_Matrix}{Matrix of p-values.}
#'   \item{Summary_Table}{Tidy data frame of pairwise correlations, p-values, and significance.}
#'   \item{Plot}{A publication-ready \code{ggplot} object.}
#' @export
#'
#' @import ggplot2
#' @importFrom stats cor cor.test sd na.omit
#' @importFrom tools toTitleCase
#' @importFrom rlang .data
correlation_heatmap <- function(
    data,
    numeric_vars = NULL,
    method = c("pearson", "spearman", "kendall"),
    shape = c("square", "circle"),
    view = c("lower", "upper", "full"),
    diag = FALSE,
    label_type = NULL,
    sig_level = NULL,
    color_palette = "RdBu",
    title = NULL
) {

  # 1. Standard Defensive Ingestion Guard
  if (is.matrix(data)) {
    df <- as.data.frame(data)
  } else if (is.data.frame(data)) {
    df <- as.data.frame(data)
  } else {
    stop("Input 'data' must be a data.frame, tibble, or matrix.")
  }

  method <- match.arg(method)
  shape  <- match.arg(shape)
  view   <- match.arg(view)

  # Default label_type: squares get both, circles get stars
  if (is.null(label_type)) {
    label_type <- if (shape == "square") "both" else "stars"
  }
  label_type <- match.arg(label_type, c("both", "stars", "values", "none"))

  # Select numeric variables
  if (is.null(numeric_vars)) {
    is_num <- vapply(df, is.numeric, logical(1))
    num_data <- df[, is_num, drop = FALSE]
  } else {
    if (!all(numeric_vars %in% names(df))) {
      missing_vars <- numeric_vars[!numeric_vars %in% names(df)]
      stop(sprintf("Specified numeric variables not found: %s", paste(missing_vars, collapse = ", ")))
    }
    num_data <- df[, numeric_vars, drop = FALSE]
  }

  if (ncol(num_data) < 2) {
    stop("At least two numeric variables are required for correlation analysis.")
  }

  # Check zero variance
  zero_var_cols <- names(num_data)[vapply(num_data, function(x) stats::sd(x, na.rm = TRUE) == 0, logical(1))]
  if (length(zero_var_cols) > 0) {
    stop(sprintf("Zero variance detected in column(s): %s. Remove constant variables.", paste(zero_var_cols, collapse = ", ")))
  }

  # Complete-case filtering
  initial_rows <- nrow(num_data)
  num_data <- stats::na.omit(num_data)
  dropped <- initial_rows - nrow(num_data)
  if (dropped > 0) {
    message(sprintf("Note: Removed %d row(s) containing NA values across selected variables.", dropped))
  }

  if (nrow(num_data) < 3) {
    warning("Fewer than 3 complete observations remain for correlation analysis.")
  }

  # 2. Compute Correlation & P-Value Matrices
  n_vars <- ncol(num_data)
  v_names <- names(num_data)
  cor_mat <- stats::cor(num_data, method = method, use = "complete.obs")
  p_mat <- matrix(0, nrow = n_vars, ncol = n_vars, dimnames = list(v_names, v_names))

  for (i in seq_len(n_vars)) {
    for (j in seq_len(n_vars)) {
      if (i != j) {
        test_ij <- stats::cor.test(num_data[[i]], num_data[[j]], method = method)
        p_mat[i, j] <- test_ij$p.value
      }
    }
  }

  # 3. Reshape into Long Format Table
  cor_long <- as.data.frame(as.table(cor_mat), stringsAsFactors = FALSE)
  p_long <- as.data.frame(as.table(p_mat), stringsAsFactors = FALSE)
  names(cor_long) <- c("Var1", "Var2", "Correlation")
  names(p_long)   <- c("Var1", "Var2", "P_Value")

  plot_df <- merge(cor_long, p_long, by = c("Var1", "Var2"))

  # Assign Significance Stars
  plot_df$Significance <- cut(
    plot_df$P_Value,
    breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
    labels = c("***", "**", "*", "ns"),
    right = TRUE
  )
  plot_df$Significance <- as.character(plot_df$Significance)
  plot_df$Significance[plot_df$Var1 == plot_df$Var2] <- ""

  # Generate formatted labels based on label_type
  stars_str <- ifelse(plot_df$Significance %in% c("***", "**", "*"), plot_df$Significance, "")

  if (label_type == "both") {
    plot_df$Label <- ifelse(plot_df$Var1 == plot_df$Var2, "1.00", sprintf("%.2f%s", plot_df$Correlation, stars_str))
  } else if (label_type == "stars") {
    plot_df$Label <- stars_str
  } else if (label_type == "values") {
    plot_df$Label <- sprintf("%.2f", plot_df$Correlation)
  } else {
    plot_df$Label <- ""
  }

  # 4. Handle Matrix Views (lower, upper, full) & Diagonal
  row_idx <- match(plot_df$Var1, v_names)
  col_idx <- match(plot_df$Var2, v_names)

  if (view == "lower") {
    valid_cells <- if (diag) row_idx >= col_idx else row_idx > col_idx
  } else if (view == "upper") {
    valid_cells <- if (diag) row_idx <= col_idx else row_idx < col_idx
  } else {
    valid_cells <- if (diag) rep(TRUE, nrow(plot_df)) else row_idx != col_idx
  }

  plot_df <- plot_df[valid_cells, , drop = FALSE]

  # Optional significance filtering (remove non-significant markers)
  if (!is.null(sig_level)) {
    insig_idx <- plot_df$P_Value > sig_level & plot_df$Var1 != plot_df$Var2
    plot_df$Correlation[insig_idx] <- NA
    plot_df$Label[insig_idx] <- ""
  }

  # Factor level ordering
  plot_df$Var1 <- factor(plot_df$Var1, levels = v_names)
  plot_df$Var2 <- factor(plot_df$Var2, levels = rev(v_names))

  # Dynamic text contrast: white on dark tiles, black on light tiles
  plot_df$TextColor <- ifelse(abs(plot_df$Correlation) > 0.65, "white", "black")
  plot_df$TextColor[is.na(plot_df$TextColor)] <- "black"

  # 5. Build ggplot2 Heatmap
  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = .data$Var1, y = .data$Var2)
  )

  if (shape == "square") {
    p <- p + ggplot2::geom_tile(
      ggplot2::aes(fill = .data$Correlation),
      color = "grey85",
      linewidth = 0.55
    )
    if (label_type != "none") {
      p <- p + ggplot2::geom_text(
        ggplot2::aes(label = .data$Label, color = .data$TextColor),
        fontface = "bold",
        size = if (label_type == "stars") 4.5 else 3.3
      ) +
        ggplot2::scale_color_identity()
    }
  } else {
    p <- p +
      ggplot2::geom_tile(fill = "transparent", color = "grey90", linewidth = 0.3) +
      ggplot2::geom_point(
        ggplot2::aes(
          size = abs(.data$Correlation),
          fill = .data$Correlation
        ),
        shape = 21,
        color = "grey35"
      ) +
      ggplot2::scale_size_continuous(range = c(4, 13), guide = "none")

    if (label_type != "none") {
      p <- p + ggplot2::geom_text(
        ggplot2::aes(label = .data$Label, color = .data$TextColor),
        fontface = "bold",
        size = if (label_type == "stars") 4.5 else 2.8
      ) +
        ggplot2::scale_color_identity()
    }
  }

  subtitle_text <- if (label_type %in% c("both", "stars")) "* p < 0.05; ** p < 0.01; *** p < 0.001" else NULL

  p <- p +
    ggplot2::scale_fill_distiller(
      palette = color_palette,
      limit = c(-1, 1),
      direction = -1,
      na.value = "grey95"
    ) +
    ggplot2::coord_fixed() +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      axis.title = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1, face = "bold", color = "black"),
      axis.text.y = ggplot2::element_text(face = "bold", color = "black"),
      panel.grid = ggplot2::element_blank(),
      legend.position = "right",
      legend.title = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(face = "bold", size = 12, hjust = 0.5),
      plot.subtitle = ggplot2::element_text(size = 9, hjust = 0.5, color = "grey30")
    ) +
    ggplot2::labs(
      fill = "r",
      title = if (!is.null(title)) title else paste0("Correlation Heatmap (", tools::toTitleCase(method), ")"),
      subtitle = subtitle_text
    )

  return(list(
    Correlation_Matrix = cor_mat,
    P_Value_Matrix = p_mat,
    Summary_Table = plot_df[, c("Var1", "Var2", "Correlation", "P_Value", "Significance")],
    Plot = p
  ))
}
