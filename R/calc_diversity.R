#' Calculate Comprehensive Alpha Diversity Indices
#'
#' Computes Richness, Abundance, Shannon, Simpson, Margalef, Menhinick,
#' and Pielou's Evenness indices simultaneously from community abundance data.
#' Automatically isolates metadata columns to prevent non-numeric errors.
#'
#' @param data A \code{data.frame}, \code{matrix}, or \code{list} of community counts.
#' @param group_col Optional integer or character vector. Explicit metadata column(s) to retain.
#'
#' @return A tidy \code{data.frame} containing metadata and computed diversity indices.
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
  } else if (is.matrix(data)) {
    data <- as.data.frame(data)
  } else if (!is.data.frame(data)) {
    stop("Input 'data' must be a data.frame, tibble, or matrix.")
  }

  # 2. Extract Metadata vs Community Matrix
  # Detect all non-numeric columns automatically
  non_num_cols <- names(data)[!vapply(data, is.numeric, logical(1))]

  if (!is.null(group_col)) {
    if (is.numeric(group_col)) {
      grp_names <- names(data)[as.integer(group_col)]
    } else {
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
    stop("Species abundance data must contain at least 2 species columns.")
  }

  if (any(is.na(comm_mat))) {
    stop("NA values detected. Replace missing values with 0 before calculating.")
  }

  # 4. Filter Empty Samples
  tot_abundance <- rowSums(comm_mat)
  valid_rows <- tot_abundance > 0

  if (!all(valid_rows)) {
    comm_mat <- comm_mat[valid_rows, , drop = FALSE]
    if (!is.null(meta_df)) meta_df <- meta_df[valid_rows, , drop = FALSE]
    tot_abundance <- rowSums(comm_mat)
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
  out_df <- data.frame(
    Sample = if (!is.null(rownames(comm_mat))) rownames(comm_mat) else seq_len(nrow(comm_mat)),
    stringsAsFactors = FALSE
  )

  if (!is.null(meta_df)) {
    out_df <- cbind(out_df, meta_df)
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
