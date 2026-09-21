#' PCA Correlation Circle with cos2 Gradient for Environmental Parameters
#'
#' Generates a publication-grade PCA correlation circle (variables factor map)
#' using base stats and ggplot2. Colors environmental parameters by their
#' representation quality (\eqn{\cos^2}) on the selected principal components.
#'
#' @param data A \code{data.frame}, \code{matrix}, or \code{list} containing
#'   metadata columns and numeric environmental variables.
#' @param exclude_cols Character vector or integer specifying metadata columns
#'   to exclude from PCA (default: \code{1}). If \code{NULL}, all columns are analyzed.
#' @param dim1 Integer. The first principal component to plot (default: 1).
#' @param dim2 Integer. The second principal component to plot (default: 2).
#' @param palette Character vector of colors for the \code{cos2} gradient
#'   (default: \code{c("#2A9D8F", "#E9C46A", "#F4A261", "#E76F51")}).
#' @param title Optional character. Plot title.
#'
#' @return A list containing:
#'   \item{Plot}{A publication-ready \code{ggplot} correlation circle object.}
#'   \item{Variable_Stats}{Data frame of variable coordinates, correlations, and \eqn{\cos^2} values.}
#'   \item{Variance_Summary}{Data frame detailing eigenvalues and percentage variance explained per PC.}
#'   \item{Sample_Scores}{Data frame of sample/station scores along the principal components.}
#'   \item{PCA_Object}{The underlying \code{prcomp} object.}
#' @export
#'
#' @import ggplot2
#' @importFrom stats prcomp sd
#' @importFrom rlang .data
pca_correlation_circle <- function(
    data,
    exclude_cols = 1,
    dim1 = 1,
    dim2 = 2,
    palette = c("#2A9D8F", "#E9C46A", "#F4A261", "#E76F51"),
    title = NULL
) {

  # 1. Standard Defensive Ingestion Guard
  if (is.list(data) && !is.data.frame(data)) {
    df <- as.data.frame(data[[1]])
  } else if (is.matrix(data)) {
    df <- as.data.frame(data)
  } else if (is.data.frame(data)) {
    df <- as.data.frame(data)
  } else {
    stop("'data' must be a data.frame, matrix, or list.")
  }

  # 2. Isolate Pure Environmental Data
  if (!is.null(exclude_cols)) {
    if (is.character(exclude_cols)) {
      if (!all(exclude_cols %in% names(df))) {
        missing_cols <- exclude_cols[!exclude_cols %in% names(df)]
        stop(sprintf("Excluded column(s) not found: %s", paste(missing_cols, collapse = ", ")))
      }
      exc_idx <- match(exclude_cols, names(df))
    } else {
      exc_idx <- as.integer(exclude_cols)
    }
    env_mat <- df[, -exc_idx, drop = FALSE]
  } else {
    env_mat <- df
  }

  # 3. Numeric & NA Validation
  non_num <- names(env_mat)[!vapply(env_mat, is.numeric, logical(1))]
  if (length(non_num) > 0) {
    stop(sprintf(
      "All analyzed columns must be numeric. Non-numeric columns detected: %s. Assign metadata to 'exclude_cols'.",
      paste(non_num, collapse = ", ")
    ))
  }

  if (any(is.na(env_mat))) {
    stop("Environmental data contains NA values. Please impute or remove missing rows before analysis.")
  }

  # Check zero variance
  zero_var_cols <- names(env_mat)[vapply(env_mat, function(x) stats::sd(x, na.rm = TRUE) == 0, logical(1))]
  if (length(zero_var_cols) > 0) {
    stop(sprintf("Zero variance detected in environmental variable(s): %s. Remove constant variables.", paste(zero_var_cols, collapse = ", ")))
  }

  if (nrow(env_mat) < 3) {
    stop("PCA correlation circle analysis requires at least 3 valid observations.")
  }

  # 4. PCA Computation (scaled to unit variance)
  pca_res <- stats::prcomp(env_mat, scale. = TRUE)
  eig <- pca_res$sdev^2
  var_explained <- (eig / sum(eig)) * 100

  max_dim <- length(eig)
  if (dim1 > max_dim || dim2 > max_dim) {
    stop(sprintf("Requested dimensions (dim1 = %d, dim2 = %d) exceed available principal components (%d).", dim1, dim2, max_dim))
  }

  pc1_lab <- sprintf("PC%d (%.1f%%)", dim1, var_explained[dim1])
  pc2_lab <- sprintf("PC%d (%.1f%%)", dim2, var_explained[dim2])

  # Calculate Variable Coordinates (loadings * sdev = correlation with PCs)
  var_coords <- sweep(pca_res$rotation, 2, pca_res$sdev, FUN = "*")
  coords_df <- as.data.frame(var_coords[, c(dim1, dim2), drop = FALSE])
  colnames(coords_df) <- c("x", "y")
  coords_df$Variable <- rownames(coords_df)

  # Calculate cos2 representation quality across selected dimensions
  cos2_mat <- var_coords^2
  coords_df$cos2 <- cos2_mat[, dim1] + cos2_mat[, dim2]

  # Generate unit circle coordinates
  theta <- seq(0, 2 * pi, length.out = 100)
  circle_df <- data.frame(x = cos(theta), y = sin(theta))

  # Variance summary table
  variance_df <- data.frame(
    Dimension = paste0("PC", seq_along(eig)),
    Eigenvalue = eig,
    Variance_Percent = var_explained,
    Cumulative_Variance = cumsum(var_explained),
    stringsAsFactors = FALSE
  )

  sample_scores <- as.data.frame(pca_res$x)

  # 5. Publication-Ready ggplot2 Rendering
  p <- ggplot2::ggplot() +
    ggplot2::geom_path(
      data = circle_df,
      ggplot2::aes(x = .data$x, y = .data$y),
      color = "grey65",
      linewidth = 0.6,
      linetype = "solid"
    ) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.4) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.4) +
    ggplot2::geom_segment(
      data = coords_df,
      ggplot2::aes(x = 0, y = 0, xend = .data$x, yend = .data$y, color = .data$cos2),
      arrow = ggplot2::arrow(length = ggplot2::unit(0.22, "cm"), type = "closed"),
      linewidth = 0.85
    ) +
    ggplot2::geom_text(
      data = coords_df,
      ggplot2::aes(x = .data$x * 1.14, y = .data$y * 1.14, label = .data$Variable, color = .data$cos2),
      size = 4.0,
      fontface = "bold",
      show.legend = FALSE
    ) +
    ggplot2::scale_color_gradientn(colors = palette, name = expression(cos^2)) +
    ggplot2::coord_fixed(ratio = 1, xlim = c(-1.20, 1.20), ylim = c(-1.20, 1.20)) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_line(color = "grey92", linewidth = 0.4),
      panel.grid.minor = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 1),
      axis.title = ggplot2::element_text(face = "bold", size = 12),
      axis.text = ggplot2::element_text(color = "black", face = "bold"),
      legend.position = "right",
      legend.title = ggplot2::element_text(face = "bold", size = 11),
      plot.title = ggplot2::element_text(face = "bold", size = 12, hjust = 0.5),
      plot.margin = ggplot2::margin(10, 10, 10, 10)
    ) +
    ggplot2::labs(
      x = pc1_lab,
      y = pc2_lab,
      title = if (!is.null(title)) title else sprintf("PCA Correlation Circle (PC%d vs PC%d)", dim1, dim2)
    )

  return(list(
    Plot = p,
    Variable_Stats = coords_df,
    Variance_Summary = variance_df,
    Sample_Scores = sample_scores,
    PCA_Object = pca_res
  ))
}
