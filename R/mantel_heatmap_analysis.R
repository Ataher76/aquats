#' Integrated Mantel Test and Correlation Heatmap
#'
#' Computes internal correlations among environmental variables and performs Mantel tests
#' linking community composition to environmental parameters using base vegan and ggplot2.
#'
#' @param comm_data A data frame containing community or species abundance data.
#' @param env_data A data frame containing numeric environmental variables.
#' @param method Character; correlation method: "pearson" or "spearman". Default is "pearson".
#' @param spec_dist Character; distance metric for community data (e.g., "bray", "euclidean"). Default is "bray".
#' @param env_dist Character; distance metric for environmental data. Default is "euclidean".
#' @param color_palette Character; color palette for heatmap. Default is "RdBu".
#'
#' @return A list containing the correlation matrix, Mantel test results, and the ggplot object.
#'
#' @import ggplot2 dplyr tidyr vegan
#' @importFrom stats cor cor.test
#' @export
mantel_heatmap_analysis <- function(comm_data, env_data, method = "pearson",
                                    spec_dist = "bray", env_dist = "euclidean",
                                    color_palette = "RdBu") {

  # 1. Environmental correlation matrix
  env_mat <- as.data.frame(env_data)
  cor_res <- stats::cor(env_mat, method = method, use = "complete.obs")

  # Melt correlation matrix for ggplot
  cor_df <- as.data.frame(as.table(cor_res))
  colnames(cor_df) <- c("Var1", "Var2", "Correlation")

  # 2. Mantel tests for each environmental variable vs community matrix
  comm_d <- vegan::vegdist(comm_data, method = spec_dist)

  mantel_results <- data.frame(Variable = character(), r = numeric(), p_value = numeric(), stringsAsFactors = FALSE)
  for (col in colnames(env_mat)) {
    single_env <- data.frame(val = env_mat[[col]])
    env_d <- vegan::vegdist(single_env, method = env_dist)
    m_test <- vegan::mantel(comm_d, env_d, method = method, permutations = 999)
    mantel_results <- rbind(mantel_results, data.frame(
      Variable = col,
      r = m_test$statistic,
      p_value = m_test$signif
    ))
  }

  # 3. Plotting correlation heatmap with Mantel bar/points
  p <- ggplot2::ggplot(cor_df, ggplot2::aes(x = Var1, y = Var2, fill = Correlation)) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::scale_fill_gradientn(
      colors = RColorBrewer::brewer.pal(n = 8, name = color_palette),
      limits = c(-1, 1)
    ) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = "Environmental Correlation & Mantel Linkages",
      subtitle = sprintf("Community Distance: %s | Env Distance: %s", spec_dist, env_dist),
      x = NULL, y = NULL, fill = "Pearson r"
    ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, face = "bold"),
      axis.text.y = ggplot2::element_text(face = "bold"),
      panel.grid = ggplot2::element_blank()
    )

  return(list(
    Correlation_Matrix = cor_res,
    Mantel_Results = mantel_results,
    Plot = p
  ))
}
