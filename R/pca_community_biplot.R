#' Principal Component Analysis (PCA) Community Biplot
#'
#' Performs Principal Component Analysis (PCA / RDA) on aquatic community data.
#' Generates a publication-ready ggplot2 biplot showing site coordinates with
#' 95% confidence ellipses and directional species loading arrows.
#'
#' @param data A \code{data.frame}, \code{matrix}, or \code{list} where grouping
#'   metadata is at the start, followed by numeric species counts.
#' @param group_col Character or integer specifying the grouping metadata column (default: 1).
#' @param scale Logical. Should species data be scaled to unit variance (default: \code{FALSE})?
#'   (\code{TRUE} is recommended when species counts span orders of magnitude).
#' @param transform Character. Optional pre-transformation via \code{vegan::decostand}:
#'   \code{"none"} (default), \code{"hellinger"}, \code{"log"}, \code{"pa"}, or \code{"wisconsin"}.
#' @param top_n_species Integer or NULL. If specified (e.g., 5), labels only the top \sQuote{n}
#'   species with the highest vector lengths (loadings) to prevent text overlap.
#' @param color_palette Character. RColorBrewer palette name for site groups (default: \code{"Dark2"}).
#' @param arrow_multiplier Numeric. Scaling factor for species arrow lengths (default: 1).
#' @param title Optional character. Plot title.
#'
#' @return A list containing:
#'   \item{Plot}{A publication-ready \code{ggplot} biplot object.}
#'   \item{Plot_Data}{Tidy data frame of site PC scores and grouping metadata.}
#'   \item{Species_Scores}{Tidy data frame of species loadings and vector coordinates.}
#'   \item{PCA_Object}{The underlying \code{rda} ordination object.}
#' @export
#'
#' @import ggplot2
#' @import vegan
#' @importFrom rlang .data
pca_community_biplot <- function(
    data,
    group_col = 1,
    scale = FALSE,
    transform = c("none", "hellinger", "log", "pa", "wisconsin"),
    top_n_species = NULL,
    color_palette = "Dark2",
    arrow_multiplier = 1,
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
    message(sprintf("Note: Removed %d sample(s) with zero total abundance.", length(zero_rows)))
  }

  if (nrow(comm_mat) < 3) {
    stop("PCA biplot analysis requires at least 3 valid samples.")
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

  # 6. PCA Computation via vegan::rda
  pca_res <- vegan::rda(comm_mat, scale = scale)

  site_scores <- vegan::scores(pca_res, display = "sites")
  plot_data <- data.frame(
    PC1 = site_scores[, 1],
    PC2 = site_scores[, 2],
    Group = group_factor,
    stringsAsFactors = FALSE
  )
  colnames(plot_data)[3] <- group_label

  species_scores <- as.data.frame(vegan::scores(pca_res, display = "species"))
  species_scores$Species <- rownames(species_scores)
  species_scores$PC1 <- species_scores$PC1 * arrow_multiplier
  species_scores$PC2 <- species_scores$PC2 * arrow_multiplier
  species_scores$Vector_Length <- sqrt(species_scores$PC1^2 + species_scores$PC2^2)

  # Filter top N species if requested
  label_df <- species_scores
  if (!is.null(top_n_species) && top_n_species < nrow(species_scores)) {
    species_scores <- species_scores[order(species_scores$Vector_Length, decreasing = TRUE), ]
    label_df <- species_scores[seq_len(top_n_species), ]
  }

  eigenvals <- pca_res$CA$eig
  var_exp <- (eigenvals / sum(eigenvals)) * 100
  pc1_lab <- sprintf("PC1 (%.1f%%)", var_exp[1])
  pc2_lab <- sprintf("PC2 (%.1f%%)", var_exp[2])

  # 7. ggplot2 Biplot Assembly
  p <- ggplot2::ggplot() +
    ggplot2::geom_point(
      data = plot_data,
      ggplot2::aes(x = .data$PC1, y = .data$PC2, color = .data[[group_label]], fill = .data[[group_label]]),
      size = 3.2, alpha = 0.85
    ) +
    ggplot2::stat_ellipse(
      data = plot_data,
      ggplot2::aes(x = .data$PC1, y = .data$PC2, color = .data[[group_label]], fill = .data[[group_label]]),
      type = "norm", level = 0.95, alpha = 0.15, geom = "polygon", show.legend = FALSE
    ) +
    ggplot2::geom_segment(
      data = species_scores,
      ggplot2::aes(x = 0, y = 0, xend = .data$PC1, yend = .data$PC2),
      inherit.aes = FALSE,
      arrow = ggplot2::arrow(length = ggplot2::unit(0.2, "cm")),
      color = "darkred", linewidth = 0.7, alpha = 0.65
    ) +
    ggplot2::geom_text(
      data = label_df,
      ggplot2::aes(x = .data$PC1, y = .data$PC2, label = .data$Species),
      inherit.aes = FALSE,
      color = "darkred", vjust = -0.6, hjust = 0.5, fontface = "bold", size = 3.5
    ) +
    ggplot2::scale_color_brewer(palette = color_palette) +
    ggplot2::scale_fill_brewer(palette = color_palette) +
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
      x = pc1_lab,
      y = pc2_lab,
      color = group_label,
      fill = group_label,
      title = if (!is.null(title)) title else "PCA Community Biplot"
    )

  return(list(
    Plot = p,
    Plot_Data = plot_data,
    Species_Scores = species_scores,
    PCA_Object = pca_res
  ))
}
