#' Non-metric Multidimensional Scaling (NMDS) Analysis for Aquatic Communities
#'
#' Performs NMDS ordination with automated data transformations (Hellinger, log, etc.),
#' stress quality checks, and generates a publication-ready ggplot2 scatterplot
#' with 95% confidence ellipses and optional environmental vector overlays.
#'
#' @param data A \code{data.frame}, \code{matrix}, or \code{list} containing community counts.
#' @param group_col Character or integer specifying the grouping metadata column (default: 1).
#' @param env_data Optional \code{data.frame} or \code{matrix} of numeric environmental variables
#'   to fit onto the ordination space via \code{vegan::envfit}.
#' @param dist_method Character. Dissimilarity metric passed to \code{vegan::vegdist}
#'   (default: \code{"bray"}).
#' @param transform Character. Community data transformation via \code{vegan::decostand}:
#'   \code{"none"} (default), \code{"hellinger"} (recommended), \code{"log"}, \code{"pa"}, or \code{"wisconsin"}.
#' @param k Integer. Number of ordination dimensions (default: 2).
#' @param trymax Integer. Maximum number of random starts for \code{metaMDS} (default: 100).
#' @param palette Character. RColorBrewer palette name for group ellipses and points (default: \code{"Dark2"}).
#' @param title Optional character. Plot title.
#'
#' @return A list containing:
#'   \item{Plot}{A publication-ready \code{ggplot} ordination object.}
#'   \item{Plot_Data}{Tidy data frame of NMDS site coordinates and grouping factors.}
#'   \item{Stress}{The final stress value of the ordination.}
#'   \item{NMDS_Object}{The underlying \code{metaMDS} object from vegan.}
#'   \item{Envfit_Results}{Vector fit results if environmental variables were supplied.}
#' @export
#'
#' @import ggplot2
#' @importFrom vegan metaMDS vegdist decostand envfit scores
#' @importFrom stats complete.cases as.formula
#' @importFrom rlang .data
nmds_analysis <- function(
    data,
    group_col = 1,
    env_data = NULL,
    dist_method = "bray",
    transform = c("none", "hellinger", "log", "pa", "wisconsin"),
    k = 2,
    trymax = 100,
    palette = "Dark2",
    title = NULL
) {

  transform <- match.arg(transform)

  # 1. Ingest Data Formats (list, matrix, or data.frame)
  if (is.list(data) && !is.data.frame(data)) {
    df <- as.data.frame(data[[1]])
  } else if (is.matrix(data)) {
    df <- as.data.frame(data)
  } else if (is.data.frame(data)) {
    df <- as.data.frame(data)
  } else {
    stop("'data' must be a data.frame, matrix, or list.")
  }

  # Ingest env_data safely
  if (!is.null(env_data)) {
    if (is.matrix(env_data)) {
      env_data <- as.data.frame(env_data)
    } else if (!is.data.frame(env_data)) {
      stop("'env_data' must be a data.frame, tibble, or matrix.")
    }
  }

  # 2. Separate Metadata and Community Matrix
  if (!is.null(group_col)) {
    if (is.character(group_col)) {
      if (!group_col %in% names(df)) stop("Grouping column not found: ", group_col)
      grp_idx <- match(group_col, names(df))
    } else {
      grp_idx <- as.integer(group_col)
    }
    meta_df <- df[, grp_idx, drop = FALSE]
    comm_mat <- df[, -grp_idx, drop = FALSE]
  } else {
    meta_df <- NULL
    comm_mat <- df
  }

  # 3. Numeric & NA Validation
  non_num <- names(comm_mat)[!vapply(comm_mat, is.numeric, logical(1))]
  if (length(non_num) > 0) {
    stop(sprintf("All species columns must be numeric. Non-numeric columns: %s", paste(non_num, collapse = ", ")))
  }
  if (any(is.na(comm_mat))) {
    stop("Species community data contains NA values. Replace with 0 or clean before analysis.")
  }

  # 4. Filter Empty Samples
  tot_abund <- rowSums(comm_mat)
  if (any(tot_abund == 0)) {
    zero_rows <- which(tot_abund == 0)
    comm_mat <- comm_mat[tot_abund > 0, , drop = FALSE]
    if (!is.null(meta_df)) meta_df <- meta_df[tot_abund > 0, , drop = FALSE]
    if (!is.null(env_data)) env_data <- env_data[tot_abund > 0, , drop = FALSE]
    message(sprintf("Note: Removed %d sample(s) with zero total abundance.", length(zero_rows)))
  }

  if (nrow(comm_mat) < 3) {
    stop("NMDS ordination requires at least 3 valid samples.")
  }

  group_factor <- if (!is.null(meta_df)) as.factor(meta_df[[1]]) else as.factor(rep("Sample", nrow(comm_mat)))
  group_label <- if (!is.null(meta_df)) names(meta_df)[1] else "Group"

  # 5. Ecological Transformations
  if (transform == "hellinger") {
    comm_mat <- vegan::decostand(comm_mat, method = "hellinger")
  } else if (transform == "log") {
    comm_mat <- vegan::decostand(comm_mat, method = "log")
  } else if (transform == "pa") {
    comm_mat <- vegan::decostand(comm_mat, method = "pa")
  } else if (transform == "wisconsin") {
    comm_mat <- vegan::decostand(comm_mat, method = "wisconsin")
  }

  # 6. NMDS Computation
  nmds_res <- vegan::metaMDS(
    comm_mat,
    distance = dist_method,
    k = k,
    trymax = trymax,
    autotransform = FALSE,
    trace = FALSE
  )

  # Stress quality assessment
  stress_val <- nmds_res$stress
  if (stress_val > 0.20) {
    warning(sprintf("High NMDS stress (%.3f > 0.20). Ordination results should be interpreted with caution.", stress_val))
  } else {
    message(sprintf("NMDS converged successfully with stress = %.3f.", stress_val))
  }

  # Extract site coordinates
  site_scores <- vegan::scores(nmds_res, display = "sites")
  plot_data <- data.frame(
    NMDS1 = site_scores[, 1],
    NMDS2 = site_scores[, 2],
    Group = group_factor,
    stringsAsFactors = FALSE
  )
  colnames(plot_data)[3] <- group_label

  # 7. Optional Environmental Vector Fitting (envfit)
  ef_res <- NULL
  ef_df <- NULL
  if (!is.null(env_data)) {
    env_df <- as.data.frame(env_data)
    env_num <- env_df[, vapply(env_df, is.numeric, logical(1)), drop = FALSE]
    if (nrow(env_num) == nrow(plot_data)) {
      ef_res <- vegan::envfit(nmds_res, env_num, permutations = 999)
      ef_vectors <- vegan::scores(ef_res, display = "vectors")
      if (!is.null(ef_vectors) && nrow(ef_vectors) > 0) {
        ef_df <- as.data.frame(ef_vectors)
        ef_df$Variable <- rownames(ef_df)
        p_vals <- ef_res$vectors$pvals
        ef_df$p_value <- p_vals[ef_df$Variable]
        ef_df <- ef_df[ef_df$p_value <= 0.05, , drop = FALSE]
      }
    }
  }

  # 8. Publication-Ready ggplot2 Rendering
  p <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = .data$NMDS1, y = .data$NMDS2, color = .data[[group_label]], fill = .data[[group_label]])
  ) +
    ggplot2::geom_point(size = 3.2, alpha = 0.85) +
    ggplot2::stat_ellipse(type = "norm", level = 0.95, alpha = 0.15, geom = "polygon", show.legend = FALSE) +
    ggplot2::scale_color_brewer(palette = palette) +
    ggplot2::scale_fill_brewer(palette = palette) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(color = "black", linewidth = 1),
      axis.title = ggplot2::element_text(face = "bold"),
      axis.text = ggplot2::element_text(color = "black", face = "bold"),
      legend.position = "right",
      legend.title = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(face = "bold", size = 12, hjust = 0.5)
    ) +
    ggplot2::labs(
      x = "NMDS Axis 1",
      y = "NMDS Axis 2",
      color = group_label,
      fill = group_label,
      title = if (!is.null(title)) title else sprintf("NMDS Ordination (Stress = %.3f)", stress_val)
    )

  # Overlay environmental arrows if significant
  if (!is.null(ef_df) && nrow(ef_df) > 0) {
    p <- p +
      ggplot2::geom_segment(
        data = ef_df,
        ggplot2::aes(x = 0, y = 0, xend = .data$NMDS1, yend = .data$NMDS2),
        inherit.aes = FALSE,
        arrow = ggplot2::arrow(length = ggplot2::unit(0.2, "cm")),
        color = "navyblue",
        linewidth = 0.8
      ) +
      ggplot2::geom_text(
        data = ef_df,
        ggplot2::aes(x = .data$NMDS1 * 1.15, y = .data$NMDS2 * 1.15, label = .data$Variable),
        inherit.aes = FALSE,
        color = "navyblue",
        fontface = "bold",
        size = 3.5
      )
  }

  return(list(
    Plot = p,
    Plot_Data = plot_data,
    Stress = stress_val,
    NMDS_Object = nmds_res,
    Envfit_Results = ef_res
  ))
}
