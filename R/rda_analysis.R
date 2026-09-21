#' Redundancy Analysis (RDA) Triplot
#'
#' Performs constrained Redundancy Analysis (RDA) combining community and environmental data.
#' Automatically cleans empty samples from both datasets, runs permutation significance tests,
#' and generates a publication-ready ggplot2 Triplot with ggrepel label adjustments.
#'
#' @param comm_data A data frame or matrix where grouping metadata is at the start,
#'   followed by numeric species counts.
#' @param env_data A data frame or matrix containing environmental variables.
#'   Must have the exact same number of rows as \code{comm_data}.
#' @param group_col Integer or character specifying the grouping column in \code{comm_data}. Default is 1.
#' @param transform Character. Optional pre-transformation: \code{"none"} (default),
#'   \code{"hellinger"}, \code{"log"}, \code{"pa"}, or \code{"wisconsin"}.
#' @param color_palette Character; a valid RColorBrewer palette name (default: "Dark2").
#' @param species_arrow_mult Numeric; scaling factor for species arrows (default: 1).
#' @param env_arrow_mult Numeric; scaling factor for environmental arrows (default: 1).
#'
#' @return A list containing:
#'   \item{RDA_Object}{The underlying \code{rda} ordination object.}
#'   \item{Model_Significance}{Permutation ANOVA test for overall model significance.}
#'   \item{Variable_Significance}{Permutation ANOVA test for marginal term significance.}
#'   \item{Plot}{The publication-ready ggplot2 triplot.}
#' @export
#'
#' @import ggplot2 vegan ggrepel
#' @importFrom stats as.formula anova sd
#' @importFrom rlang .data
rda_analysis <- function(comm_data,
                         env_data,
                         group_col = 1,
                         transform = c("none", "hellinger", "log", "pa", "wisconsin"),
                         color_palette = "Dark2",
                         species_arrow_mult = 1,
                         env_arrow_mult = 1) {

  transform <- match.arg(transform)

  # --- STEP 1: DEFENSIVE INGESTION & MATRIX COERCION ---
  if (is.matrix(comm_data)) {
    comm_data <- as.data.frame(comm_data)
  } else if (!is.data.frame(comm_data)) {
    stop("Input 'comm_data' must be a data.frame, tibble, or matrix.")
  }

  if (is.matrix(env_data)) {
    env_data <- as.data.frame(env_data)
  } else if (!is.data.frame(env_data)) {
    stop("Input 'env_data' must be a data.frame, tibble, or matrix.")
  }

  if (nrow(comm_data) != nrow(env_data)) {
    stop("comm_data and env_data must have the exact same number of rows (samples).")
  }

  # Isolate grouping metadata from community data
  if (is.character(group_col)) {
    if (!group_col %in% names(comm_data)) stop("Grouping column not found in comm_data: ", group_col)
    grp_idx <- match(group_col, names(comm_data))
  } else {
    grp_idx <- as.integer(group_col)
  }

  meta_df <- comm_data[, grp_idx, drop = FALSE]
  comm_mat <- comm_data[, -grp_idx, drop = FALSE]

  # Isolate purely numeric columns from environmental data
  env_mat <- env_data[, vapply(env_data, is.numeric, logical(1)), drop = FALSE]
  if (ncol(env_mat) == 0) stop("No numeric environmental variables found in env_data.")

  if (any(is.na(comm_mat)) || any(is.na(env_mat))) {
    stop("Data contains NA values. Impute or remove missing values before analysis.")
  }

  # Auto-clean zero abundance samples from BOTH matrices to maintain alignment
  n_abund <- rowSums(comm_mat)
  if (any(n_abund == 0)) {
    zero_rows <- which(n_abund == 0)
    comm_mat <- comm_mat[n_abund > 0, , drop = FALSE]
    env_mat  <- env_mat[n_abund > 0, , drop = FALSE]
    meta_df  <- meta_df[n_abund > 0, , drop = FALSE]
    message(sprintf("Note: Removed %d sample(s) with zero abundance from both datasets.", length(zero_rows)))
  }

  if (nrow(comm_mat) < 4) {
    stop("RDA requires at least 4 valid samples.")
  }

  group_factor <- as.factor(meta_df[[1]])
  group_name   <- names(meta_df)[1]

  # --- STEP 2: TRANSFORMATIONS & RDA COMPUTATION ---
  if (transform == "hellinger") {
    comm_mat <- vegan::decostand(comm_mat, method = "hellinger")
  } else if (transform == "log") {
    comm_mat <- vegan::decostand(comm_mat, method = "log")
  } else if (transform == "pa") {
    comm_mat <- vegan::decostand(comm_mat, method = "pa")
  } else if (transform == "wisconsin") {
    comm_mat <- vegan::decostand(comm_mat, method = "wisconsin")
  }

  # Scale environmental variables to unit variance for balanced regression weights
  env_mat_scaled <- as.data.frame(scale(env_mat))

  # Run the constrained RDA model
  rda_res <- vegan::rda(comm_mat ~ ., data = env_mat_scaled)

  # Run permutation tests for overall model and marginal variable significance
  model_sig <- stats::anova(rda_res, permutations = 999)
  terms_sig <- stats::anova(rda_res, by = "margin", permutations = 999)

  # --- STEP 3: SCORES EXTRACTION ---
  eigenvals <- rda_res$CCA$eig
  tot_var   <- rda_res$tot.chi
  var_exp   <- (eigenvals / tot_var) * 100
  rda1_lab  <- sprintf("RDA1 (%.1f%%)", var_exp[1])
  rda2_lab  <- sprintf("RDA2 (%.1f%%)", var_exp[2])

  site_scores <- as.data.frame(vegan::scores(rda_res, display = "sites"))
  site_scores$Group <- group_factor

  species_scores <- as.data.frame(vegan::scores(rda_res, display = "species"))
  if (nrow(species_scores) > 0) {
    species_scores$Species <- rownames(species_scores)
    species_scores$RDA1 <- species_scores$RDA1 * species_arrow_mult
    species_scores$RDA2 <- species_scores$RDA2 * species_arrow_mult
  }

  env_scores <- as.data.frame(vegan::scores(rda_res, display = "bp"))
  env_scores$Variable <- rownames(env_scores)
  env_scores$RDA1 <- env_scores$RDA1 * env_arrow_mult
  env_scores$RDA2 <- env_scores$RDA2 * env_arrow_mult

  # --- STEP 4: PUBLICATION TRIPLOT VISUALIZATION ---
  p <- ggplot2::ggplot() +
    ggplot2::geom_point(
      data = site_scores,
      ggplot2::aes(x = .data$RDA1, y = .data$RDA2, color = .data$Group, fill = .data$Group),
      size = 3.2, alpha = 0.85
    ) +
    ggplot2::stat_ellipse(
      data = site_scores,
      ggplot2::aes(x = .data$RDA1, y = .data$RDA2, color = .data$Group, fill = .data$Group),
      type = "norm", level = 0.95, alpha = 0.15, geom = "polygon", show.legend = FALSE
    ) +
    {if (nrow(species_scores) > 0) list(
      ggplot2::geom_segment(
        data = species_scores,
        ggplot2::aes(x = 0, y = 0, xend = .data$RDA1, yend = .data$RDA2),
        inherit.aes = FALSE,
        arrow = ggplot2::arrow(length = ggplot2::unit(0.16, "cm")),
        color = "darkred", linewidth = 0.55, alpha = 0.45
      ),
      ggrepel::geom_text_repel(
        data = species_scores,
        ggplot2::aes(x = .data$RDA1, y = .data$RDA2, label = .data$Species),
        inherit.aes = FALSE,
        color = "darkred", fontface = "italic", size = 3.6, box.padding = 0.4, max.overlaps = 50
      )
    )} +
    ggplot2::geom_segment(
      data = env_scores,
      ggplot2::aes(x = 0, y = 0, xend = .data$RDA1, yend = .data$RDA2),
      inherit.aes = FALSE,
      arrow = ggplot2::arrow(length = ggplot2::unit(0.22, "cm"), type = "closed"),
      color = "navyblue", linewidth = 0.95
    ) +
    ggrepel::geom_text_repel(
      data = env_scores,
      ggplot2::aes(x = .data$RDA1, y = .data$RDA2, label = .data$Variable),
      inherit.aes = FALSE,
      color = "navyblue", fontface = "bold", size = 4.2, box.padding = 0.5, max.overlaps = Inf
    ) +
    ggplot2::scale_color_brewer(palette = color_palette, name = group_name) +
    ggplot2::scale_fill_brewer(palette = color_palette, name = group_name) +
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
      x = rda1_lab,
      y = rda2_lab,
      title = "Redundancy Analysis (RDA) Triplot"
    )

  return(list(
    RDA_Object = rda_res,
    Model_Significance = model_sig,
    Variable_Significance = terms_sig,
    Plot = p
  ))
}
