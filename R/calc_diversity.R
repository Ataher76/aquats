#' Calculate Comprehensive Alpha Diversity Indices
#'
#' Computes species richness, total abundance, Shannon-Wiener (\eqn{H'}),
#' Gini-Simpson (\eqn{1 - D}), Margalef (\eqn{d}), Menhinick (\eqn{R}),
#' and Pielou's evenness (\eqn{J'}) indices simultaneously from community
#' abundance matrices. Automatically isolates non-numeric metadata columns,
#' guards against mathematical domain errors (e.g., division by zero), and
#' removes empty sample observations.
#'
#' @param data A \code{data.frame}, \code{matrix}, or single-element \code{list}
#'   containing community abundance counts and optional metadata columns.
#' @param group_col Optional character or integer vector; column name(s) or index/indices
#'   in \code{data} specifying metadata variables to retain alongside calculated diversity
#'   indices. Default is \code{NULL}.
#'
#' @return A tidy \code{data.frame} containing sample identifiers, retained metadata
#'   columns, and the following computed diversity metrics:
#'   \item{Sample}{Sample row names or integer index identifiers.}
#'   \item{Abundance}{Total observed individual count/abundance per sample (\eqn{N}).}
#'   \item{Richness}{Observed species richness count (\eqn{S}).}
#'   \item{Shannon}{Shannon-Wiener diversity index (\eqn{H'}, natural log base \eqn{e}).}
#'   \item{Simpson}{Gini-Simpson diversity index (\eqn{1 - D}).}
#'   \item{Margalef}{Margalef richness index (\eqn{d = (S - 1) / \ln(N)}).}
#'   \item{Menhinick}{Menhinick richness index (\eqn{R = S / \sqrt{N}}).}
#'   \item{Pielou}{Pielou's evenness index (\eqn{J' = H' / \ln(S)}).}
#' @export
#'
#' @importFrom vegan diversity specnumber
calc_diversity <- function(
    data,
    group_col = NULL
) {

  # 1. Standard Defensive Ingestion Guard
  if (is.list(data) && !is.data.frame(data)) {
    data <- as.data.frame(data[[1]])
  } else if (is.matrix(data) || inherits(data, "tbl_df") || inherits(data, "tbl")) {
    data <- as.data.frame(data)
  } else if (!is.data.frame(data)) {
    stop("Input 'data' must be a data.frame, tibble, or matrix.")
  }

  # Ensure base data.frame behavior (prevents tibble 1D drop anomalies)
  data <- as.data.frame(data)

  # 2. Extract Metadata vs Community Matrix
  non_num_cols <- names(data)[!vapply(data, is.numeric, logical(1))]

  if (!is.null(group_col)) {
    if (is.numeric(group_col)) {
      invalid_idx <- group_col < 1 | group_col > ncol(data)
      if (any(invalid_idx)) {
        stop("One or more column indices in 'group_col' are out of bounds.")
      }
      grp_names <- names(data)[as.integer(group_col)]
    } else {
      missing_cols <- setdiff(group_col, names(data))
      if (length(missing_cols) > 0) {
        stop(sprintf("Column(s) not found in 'data': %s", paste(missing_cols, collapse = ", ")))
      }
      grp_names <- group_col
    }
    all_meta_cols <- unique(c(grp_names, non_num_cols))
  } else {
    all_meta_cols <- non_num_cols
  }

  if (length(all_meta_cols) > 0) {
    meta_df <- data[, all_meta_cols, drop = FALSE]
    comm_mat <- data[, !names(data) %in% all_meta_cols, drop = FALSE]
  } else {
    meta_df <- NULL
    comm_mat <- data
  }

  # 3. Numeric Validation
  if (ncol(comm_mat) < 2) {
    stop("Species abundance data must contain at least 2 numeric species columns.")
  }

  if (anyNA(comm_mat)) {
    stop("NA values detected in community counts. Replace missing values with 0 before calculating.")
  }

  # 4. Filter Empty Samples
  tot_abundance <- rowSums(comm_mat)
  valid_rows <- tot_abundance > 0

  if (!all(valid_rows)) {
    n_dropped <- sum(!valid_rows)
    warning(sprintf("Removed %d empty sample(s) with zero total abundance.", n_dropped), call. = FALSE)
    comm_mat <- comm_mat[valid_rows, , drop = FALSE]
    if (!is.null(meta_df)) meta_df <- meta_df[valid_rows, , drop = FALSE]
    tot_abundance <- rowSums(comm_mat)
  }

  if (nrow(comm_mat) == 0) {
    stop("No valid samples remaining after filtering zero-abundance rows.")
  }

  # 5. Compute Alpha Diversity Indices
  s_rich  <- vegan::specnumber(comm_mat)
  n_abund <- tot_abundance
  h_shan  <- vegan::diversity(comm_mat, index = "shannon")
  d_simp  <- vegan::diversity(comm_mat, index = "simpson")
  d_marg  <- ifelse(n_abund > 1, (s_rich - 1) / log(n_abund), 0)
  d_menh  <- ifelse(n_abund > 0, s_rich / sqrt(n_abund), 0)
  j_piel  <- ifelse(s_rich > 1, h_shan / log(s_rich), 0)

  # 6. Assemble Output
  # Avoid duplicate 'Sample' column if already supplied in metadata
  has_sample_col <- !is.null(meta_df) && ("Sample" %in% names(meta_df))

  if (!has_sample_col) {
    sample_ids <- if (!is.null(rownames(comm_mat))) rownames(comm_mat) else seq_len(nrow(comm_mat))
    out_df <- data.frame(Sample = sample_ids, stringsAsFactors = FALSE)
    if (!is.null(meta_df)) {
      out_df <- cbind(out_df, meta_df)
    }
  } else {
    out_df <- meta_df
  }

  indices_df <- data.frame(
    Abundance = round(n_abund, 2),
    Richness  = as.integer(s_rich),
    Shannon   = round(h_shan, 3),
    Simpson   = round(d_simp, 3),
    Margalef  = round(d_marg, 3),
    Menhinick = round(d_menh, 3),
    Pielou    = round(j_piel, 3),
    stringsAsFactors = FALSE
  )

  return(cbind(out_df, indices_df))
}
