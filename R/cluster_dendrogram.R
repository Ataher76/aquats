#' Hierarchical Clustering Dendrogram for Aquatic Communities
#'
#' Computes ecological dissimilarity matrices (Bray-Curtis, Jaccard, Euclidean, etc.)
#' with optional multivariate transformations (Hellinger, log, presence/absence, standardization)
#' and performs agglomerative hierarchical clustering (UPGMA, Ward.D2, Complete).
#' Generates a publication-ready ggplot2 dendrogram with color-coded group tips,
#' a cluster cut-off threshold line, and a sample-to-cluster assignment table.
#'
#' @param data A \code{data.frame}, \code{matrix}, or \code{list} containing
#'   community abundance counts or environmental parameters with optional grouping metadata.
#' @param group_col Character or integer specifying the grouping metadata column (default: 1).
#'   If \code{NULL}, clustering runs without group categorization.
#' @param dist_method Character. Dissimilarity metric passed to \code{vegan::vegdist}:
#'   \code{"bray"} (default), \code{"jaccard"}, \code{"euclidean"}, \code{"horn"},
#'   \code{"kulczynski"}, \code{"gower"}, \code{"manhattan"}, or \code{"canberra"}.
#' @param transform Character. Transformation applied via \code{vegan::decostand}:
#'   \code{"none"} (default), \code{"hellinger"} (recommended for abundance data),
#'   \code{"log"} (\eqn{\ln(x + 1)}), \code{"pa"} (presence/absence binary),
#'   \code{"standardize"} (z-score, recommended for water chemistry), or \code{"sqrt"}.
#' @param cluster_method Character. Agglomeration method passed to \code{stats::hclust}:
#'   \code{"average"} (UPGMA, default), \code{"ward.D2"}, \code{"ward.D"},
#'   \code{"complete"}, \code{"single"}, or \code{"centroid"}.
#' @param k Optional integer. Number of clusters to cut into (defaults to number of group levels).
#' @param h Optional numeric. Specific distance height threshold to cut the dendrogram.
#' @param palette Character. RColorBrewer palette name for group leaf tips (default: \code{"Set2"}).
#' @param orientation Character. Tree layout: \code{"vertical"} (default) or \code{"horizontal"}.
#' @param title Optional character. Plot title.
#'
#' @return A list containing:
#'   \item{Plot}{A publication-ready \code{ggplot} dendrogram object.}
#'   \item{Cluster_Assignments}{Tidy data frame showing sample cluster memberships.}
#'   \item{Distance_Matrix}{The computed \code{dist} object.}
#'   \item{HClust_Object}{The underlying \code{hclust} object.}
#'   \item{Transformed_Data}{The processed numeric matrix used for distance computation.}
#' @export
#'
#' @import ggplot2
#' @importFrom vegan vegdist decostand
#' @importFrom stats hclust cutree
#' @importFrom tools toTitleCase
#' @importFrom rlang .data
cluster_dendrogram <- function(
    data,
    group_col = 1,
    dist_method = c("bray", "jaccard", "euclidean", "horn", "kulczynski", "gower", "manhattan", "canberra"),
    transform = c("none", "hellinger", "log", "pa", "standardize", "sqrt"),
    cluster_method = c("average", "ward.D2", "ward.D", "complete", "single", "centroid"),
    k = NULL,
    h = NULL,
    palette = "Set2",
    orientation = c("vertical", "horizontal"),
    title = NULL
) {

  # 1. Standard Defensive Ingestion Guard
  if (is.list(data) && !is.data.frame(data)) {
    grp_name <- names(data)[1]
    message(sprintf("Note: 'data' is a list. Analyzing first element: '%s'.", grp_name))
    df <- as.data.frame(data[[1]])
  } else if (is.matrix(data)) {
    df <- as.data.frame(data)
  } else if (is.data.frame(data)) {
    df <- as.data.frame(data)
  } else {
    stop("'data' must be a data.frame, matrix, or list.")
  }

  dist_method <- match.arg(tolower(dist_method), c("bray", "jaccard", "euclidean", "horn", "kulczynski", "gower", "manhattan", "canberra"))
  transform   <- match.arg(tolower(transform), c("none", "hellinger", "log", "pa", "standardize", "sqrt"))
  orientation <- match.arg(orientation)

  # Map case-sensitive hclust method strings
  hc_methods <- c(
    "average"  = "average",
    "ward.d2"  = "ward.D2",
    "ward.d"   = "ward.D",
    "complete" = "complete",
    "single"   = "single",
    "centroid" = "centroid"
  )
  cluster_method_raw <- match.arg(tolower(cluster_method), names(hc_methods))
  hc_method <- hc_methods[[cluster_method_raw]]

  # 2. Separate Grouping Metadata from Numeric Matrix
  if (!is.null(group_col)) {
    if (is.character(group_col)) {
      if (!all(group_col %in% names(df))) {
        missing_cols <- group_col[!group_col %in% names(df)]
        stop(sprintf("Grouping column(s) not found: %s", paste(missing_cols, collapse = ", ")))
      }
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

  # 3. Numeric and NA Validation
  non_num <- names(comm_mat)[!vapply(comm_mat, is.numeric, logical(1))]
  if (length(non_num) > 0) {
    stop(sprintf(
      "All community columns must be numeric. Non-numeric columns detected: %s. Place metadata in grouping columns.",
      paste(non_num, collapse = ", ")
    ))
  }

  if (any(is.na(comm_mat))) {
    stop("Community matrix contains NA values. Please replace missing entries with 0 before analysis.")
  }

  if (!is.null(meta_df) && any(is.na(meta_df))) {
    stop("Grouping metadata contains NA values. Please remove or impute missing group entries.")
  }

  # 4. Filter Empty Samples
  tot_abund <- rowSums(comm_mat)
  if (any(tot_abund == 0)) {
    zero_rows <- which(tot_abund == 0)
    comm_mat <- comm_mat[tot_abund > 0, , drop = FALSE]
    if (!is.null(meta_df)) meta_df <- meta_df[tot_abund > 0, , drop = FALSE]
    message(sprintf("Note: Removed %d sample(s) with zero total abundance.", length(zero_rows)))
  }

  n_samples <- nrow(comm_mat)
  if (n_samples < 2) {
    stop("Hierarchical clustering requires at least 2 valid samples after data cleaning.")
  }

  # 5. Ecological Data Transformations
  if (transform == "hellinger") {
    comm_mat <- vegan::decostand(comm_mat, method = "hellinger")
  } else if (transform == "log") {
    comm_mat <- vegan::decostand(comm_mat, method = "log")
  } else if (transform == "pa") {
    comm_mat <- vegan::decostand(comm_mat, method = "pa")
  } else if (transform == "standardize") {
    comm_mat <- vegan::decostand(comm_mat, method = "standardize", MARGIN = 2)
  } else if (transform == "sqrt") {
    comm_mat <- sqrt(comm_mat)
  }

  # 6. Distance Matrix & Hierarchical Clustering
  is_binary <- (transform == "pa")
  dist_mat <- vegan::vegdist(comm_mat, method = dist_method, binary = is_binary)
  hc <- stats::hclust(dist_mat, method = hc_method)

  sample_labels <- if (!is.null(rownames(comm_mat))) rownames(comm_mat) else paste0("S_", seq_len(n_samples))
  grp_factor <- if (!is.null(meta_df)) as.factor(meta_df[[1]]) else as.factor(rep("Sample", n_samples))
  grp_label_name <- if (!is.null(meta_df)) names(meta_df)[1] else "Group"

  # 7. Cluster Cut Assignments
  if (is.null(k) && is.null(h)) {
    k <- length(unique(grp_factor))
  }

  if (!is.null(h)) {
    cluster_vec <- stats::cutree(hc, h = h)
    cut_height <- h
  } else if (!is.null(k) && k > 1 && k < n_samples) {
    cluster_vec <- stats::cutree(hc, k = k)
    cut_height <- mean(c(hc$height[n_samples - k + 1], hc$height[n_samples - k]))
  } else if (!is.null(k) && k == n_samples) {
    cluster_vec <- stats::cutree(hc, k = k)
    cut_height <- hc$height[1] / 2
  } else {
    cluster_vec <- rep(1, n_samples)
    cut_height <- NULL
  }

  clusters_df <- data.frame(
    Sample = sample_labels,
    Group = grp_factor,
    Cluster = paste0("Cluster_", cluster_vec),
    stringsAsFactors = FALSE
  )
  colnames(clusters_df)[2] <- grp_label_name

  # 8. Extract Coordinates for ggplot2 Dendrogram
  leaf_x <- integer(n_samples)
  leaf_x[hc$order] <- seq_len(n_samples)

  node_x <- numeric(2 * n_samples - 1)
  node_y <- numeric(2 * n_samples - 1)
  names(node_x) <- as.character(c(-seq_len(n_samples), seq_len(n_samples - 1)))
  names(node_y) <- as.character(c(-seq_len(n_samples), seq_len(n_samples - 1)))

  for (i in seq_len(n_samples)) {
    node_x[as.character(-i)] <- leaf_x[i]
    node_y[as.character(-i)] <- 0
  }

  n_merges <- nrow(hc$merge)
  seg_x <- numeric(3 * n_merges)
  seg_y <- numeric(3 * n_merges)
  seg_xend <- numeric(3 * n_merges)
  seg_yend <- numeric(3 * n_merges)
  idx <- 1

  for (m in seq_len(n_merges)) {
    a <- hc$merge[m, 1]
    b <- hc$merge[m, 2]
    h_val <- hc$height[m]

    xa <- node_x[as.character(a)]
    ya <- node_y[as.character(a)]
    xb <- node_x[as.character(b)]
    yb <- node_y[as.character(b)]

    xc <- (xa + xb) / 2
    node_x[as.character(m)] <- xc
    node_y[as.character(m)] <- h_val

    # Vertical line child a
    seg_x[idx] <- xa; seg_y[idx] <- ya; seg_xend[idx] <- xa; seg_yend[idx] <- h_val; idx <- idx + 1
    # Vertical line child b
    seg_x[idx] <- xb; seg_y[idx] <- yb; seg_xend[idx] <- xb; seg_yend[idx] <- h_val; idx <- idx + 1
    # Horizontal bar
    seg_x[idx] <- xa; seg_y[idx] <- h_val; seg_xend[idx] <- xb; seg_yend[idx] <- h_val; idx <- idx + 1
  }

  segments_df <- data.frame(
    x = seg_x, y = seg_y, xend = seg_xend, yend = seg_yend
  )

  leaves_df <- data.frame(
    x = seq_len(n_samples),
    y = 0,
    Label = sample_labels[hc$order],
    Group = grp_factor[hc$order],
    stringsAsFactors = FALSE
  )

  # 9. Publication-Ready Dendrogram Assembly
  max_height <- max(hc$height)

  y_axis_label <- if (transform != "none") {
    paste0(tools::toTitleCase(dist_method), " Distance (", tools::toTitleCase(transform), ")")
  } else {
    paste(tools::toTitleCase(dist_method), "Distance")
  }

  p <- ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = segments_df,
      ggplot2::aes(x = .data$x, y = .data$y, xend = .data$xend, yend = .data$yend),
      color = "grey25",
      linewidth = 0.65
    )

  if (!is.null(cut_height)) {
    p <- p + ggplot2::geom_hline(
      yintercept = cut_height,
      linetype = "dashed",
      color = "#D9534F",
      linewidth = 0.75
    )
  }

  p <- p +
    ggplot2::geom_point(
      data = leaves_df,
      ggplot2::aes(x = .data$x, y = .data$y, color = .data$Group),
      size = 3.2,
      shape = 16
    ) +
    ggplot2::scale_color_brewer(palette = palette) +
    ggplot2::scale_x_continuous(
      breaks = leaves_df$x,
      labels = leaves_df$Label,
      expand = ggplot2::expansion(add = c(0.8, 0.8))
    ) +
    ggplot2::scale_y_continuous(
      limits = c(0, max_height * 1.12),
      expand = ggplot2::expansion(mult = c(0.02, 0.05))
    ) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(
      axis.line = ggplot2::element_line(linewidth = 0.7, color = "black"),
      axis.title = ggplot2::element_text(face = "bold"),
      axis.text = ggplot2::element_text(color = "black", face = "bold"),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      legend.position = if (!is.null(meta_df)) "top" else "none",
      legend.title = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(face = "bold", size = 12, hjust = 0.5)
    ) +
    ggplot2::labs(
      x = "Sampling Stations",
      y = y_axis_label,
      color = grp_label_name,
      title = if (!is.null(title)) title else paste0("Hierarchical Clustering Dendrogram (", tools::toTitleCase(cluster_method_raw), ")")
    )

  if (orientation == "horizontal") {
    p <- p +
      ggplot2::coord_flip() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5))
  }

  return(list(
    Plot = p,
    Cluster_Assignments = clusters_df,
    Distance_Matrix = dist_mat,
    HClust_Object = hc,
    Transformed_Data = comm_mat
  ))
}
