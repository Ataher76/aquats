#' Similarity Percentages (SIMPER) Analysis
#'
#' Identifies which species primarily contribute to the differences between groups
#' using Bray-Curtis dissimilarity decomposition (\code{vegan::simper}).
#' Automatically cleans empty samples, validates missing values, supports ecological
#' transformations, and generates a publication-ready ggplot2 bar chart.
#'
#' @param data A \code{data.frame}, \code{matrix}, or \code{list} where grouping
#'   metadata is at the start, followed by numeric species counts.
#' @param group_col Character or integer specifying the grouping metadata column (default: 1).
#' @param top_n Integer. Number of top contributing species to display per comparison (default: 5).
#' @param transform Character. Optional pre-transformation via \code{vegan::decostand}:
#'   \code{"none"} (default), \code{"hellinger"} (recommended), \code{"log"}, \code{"pa"}, or \code{"wisconsin"}.
#' @param permutations Integer. Number of permutations for p-value calculation (default: 999).
#' @param seed Optional integer. Random seed for reproducible permutations (default: 42).
#'
#' @return A list containing:
#'   \item{SIMPER_Object}{The raw \code{simper} object from vegan.}
#'   \item{Summary_Table}{Tidy summary data frame of top contributing species per group comparison.}
#'   \item{Plot}{A publication-ready \code{ggplot} bar chart of top contributors.}
#' @export
#'
#' @import ggplot2
#' @import vegan
#' @importFrom stats reorder
#' @importFrom utils head
#' @importFrom rlang .data
simper_analysis <- function(
    data,
    group_col = 1,
    top_n = 5,
    transform = c("none", "hellinger", "log", "pa", "wisconsin"),
    permutations = 999,
    seed = 42
) {

  transform <- match.arg(transform)
  if (!is.null(seed)) set.seed(seed)

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
    stop("SIMPER analysis requires a grouping column to compare groups.")
  }

  # 3. Numeric & NA Validation
  non_num <- names(comm_mat)[!vapply(comm_mat, is.numeric, logical(1))]
  if (length(non_num) > 0) {
    stop(sprintf("All community columns must be numeric. Non-numeric columns: %s", paste(non_num, collapse = ", ")))
  }

  if (any(is.na(comm_mat))) {
    stop("Species community data contains NA values. Replace with 0 or clean before analysis.")
  }
  if (any(is.na(meta_df))) {
    stop("Grouping metadata contains NA values. Please remove or impute missing group values.")
  }

  # 4. Filter Empty Samples
  tot_abund <- rowSums(comm_mat)
  if (any(tot_abund == 0)) {
    zero_rows <- which(tot_abund == 0)
    comm_mat <- comm_mat[tot_abund > 0, , drop = FALSE]
    meta_df <- meta_df[tot_abund > 0, , drop = FALSE]
    message(sprintf("Note: Removed %d sample(s) with zero total abundance.", length(zero_rows)))
  }

  if (nrow(comm_mat) < 3) {
    stop("SIMPER analysis requires at least 3 valid samples.")
  }

  group_factor <- as.factor(meta_df[[1]])
  if (nlevels(group_factor) < 2) {
    stop("SIMPER analysis requires at least 2 distinct groups in the grouping factor.")
  }

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

  # 6. Run vegan SIMPER
  simp <- vegan::simper(comm_mat, group_factor, permutations = permutations)
  simp_sum <- summary(simp, ordered = TRUE)

  # 7. Extract Top N Contributors per Comparison
  plot_data <- data.frame()
  for (comp in names(simp_sum)) {
    comp_df <- simp_sum[[comp]]
    top_species <- utils::head(comp_df, top_n)

    # Safely resolve cumulative column name across vegan versions ('cusum' vs 'cumsum')
    cum_col_vals <- if (!is.null(top_species$cusum)) top_species$cusum else top_species$cumsum
    if (is.null(cum_col_vals)) cum_col_vals <- rep(NA, nrow(top_species))

    temp_df <- data.frame(
      Comparison = comp,
      Species = rownames(top_species),
      Contribution = top_species$average,
      Cumulative = cum_col_vals,
      p_value = if (!is.null(top_species$p)) top_species$p else NA,
      stringsAsFactors = FALSE
    )
    plot_data <- rbind(plot_data, temp_df)
  }

  if (nrow(plot_data) == 0) {
    warning("No valid SIMPER comparison results extracted.")
  }

  # 8. Publication-Ready ggplot2 Bar Chart
  p <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = stats::reorder(.data$Species, .data$Contribution), y = .data$Contribution, fill = .data$Comparison)
  ) +
    ggplot2::geom_col(width = 0.65, color = "black", linewidth = 0.35, alpha = 0.85) +
    ggplot2::coord_flip() +
    ggplot2::facet_wrap(~ .data$Comparison, scales = "free_y") +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(color = "black", linewidth = 1),
      axis.title = ggplot2::element_text(face = "bold"),
      axis.text = ggplot2::element_text(color = "black", face = "bold"),
      legend.position = "none",
      strip.background = ggplot2::element_rect(fill = "grey90", color = "black"),
      strip.text = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(face = "bold", size = 12, hjust = 0.5)
    ) +
    ggplot2::labs(
      x = "Species / Taxa",
      y = "Average Contribution to Dissimilarity",
      title = "SIMPER: Top Contributing Taxa Across Pairwise Group Comparisons"
    )

  return(list(
    SIMPER_Object = simp,
    Summary_Table = plot_data,
    Plot = p
  ))
}
