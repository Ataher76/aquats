#' Permutational Multivariate Analysis of Variance (PERMANOVA)
#'
#' Performs PERMANOVA (using \code{vegan::adonis2}) to test for statistical differences
#' in community composition between groups. Automatically cleans empty samples, validates
#' numeric matrices, supports ecological transformations, and handles flexible data ingestion.
#'
#' @param data A \code{data.frame} or \code{matrix} containing grouping metadata
#'   and numeric species abundance counts.
#' @param group_col Character string or integer; column name or index in \code{data}
#'   specifying the categorical grouping factor. Default is \code{1}.
#' @param dist_method Character string; dissimilarity or distance metric passed to
#'   \code{vegan::vegdist} (e.g., \code{"bray"}, \code{"jaccard"}, \code{"euclidean"}).
#'   Default is \code{"bray"}.
#' @param transform Character string; pre-transformation applied to community counts via
#'   \code{vegan::decostand}: \code{"none"} (default), \code{"hellinger"} (recommended for
#'   abundance data), \code{"log"}, \code{"pa"} (presence/absence), or \code{"wisconsin"}.
#'   Default is \code{"none"}.
#' @param permutations Numeric integer; number of Monte Carlo permutations for significance
#'   testing. Default is \code{999}.
#' @param seed Optional integer; random seed for reproducible permutation testing.
#'   Default is \code{42}.
#'
#' @return A list containing:
#'   \item{PERMANOVA_Table}{The ANOVA-style permutation results table (\code{anova.cca}) from \code{vegan::adonis2}.}
#'   \item{Cleaned_Data}{A data frame combining the filtered grouping factor and community matrix.}
#'   \item{Distance_Method}{Character string indicating the dissimilarity metric applied.}
#' @export
#'
#' @import vegan
#' @importFrom stats complete.cases as.formula
permanova_analysis <- function(
    data,
    group_col = 1,
    dist_method = "bray",
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
    stop("PERMANOVA requires a grouping column to test for differences between groups.")
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
    stop("PERMANOVA analysis requires at least 3 valid samples.")
  }

  group_factor <- as.factor(meta_df[[1]])
  if (nlevels(group_factor) < 2) {
    stop("PERMANOVA requires at least 2 distinct groups in the grouping factor.")
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

  # 6. PERMANOVA Computation via vegan::adonis2
  design_df <- data.frame(Group = group_factor)
  formula_obj <- stats::as.formula("comm_mat ~ Group")
  perm_res <- vegan::adonis2(
    formula_obj,
    data = design_df,
    method = dist_method,
    permutations = permutations
  )

  return(list(
    PERMANOVA_Table = perm_res,
    Cleaned_Data = cbind(Group = group_factor, comm_mat),
    Distance_Method = dist_method
  ))
}
