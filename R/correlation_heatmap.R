#' Flexible Correlation Matrix and Heatmap Generator
#'
#' Computes correlation matrices (Pearson, Spearman, or Kendall) and generates
#' publication-ready heatmaps with customizable shapes (circles, squares), layout
#' views (full, lower, upper), and optional text labels.
#'
#' @param data A data frame containing numeric variables.
#' @param numeric_vars Optional character vector of numeric column names to include. If NULL, all numeric columns are used.
#' @param method Character; correlation method: "pearson", "spearman", or "kendall". Default is "pearson".
#' @param shape Character; marker shape for the heatmap: "circle" or "square". Default is "circle".
#' @param view Character; matrix display layout: "full", "lower", or "upper". Default is "lower".
#' @param show_text Logical; whether to display numeric correlation values inside markers. Default is FALSE for circles and TRUE for squares.
#' @param color_palette Character; name of a diverging color palette from RColorBrewer. Default is "RdBu".
#'
#' @return A list containing the correlation matrix, p-value matrix, and the ggplot object.
#'
#' @import ggplot2 RColorBrewer dplyr tidyr
#' @importFrom stats cor cor.test
#' @importFrom dplyr %>% select where
#' @export
#'
correlation_heatmap <- function(data,
                                numeric_vars = NULL,
                                method = "pearson",
                                shape = "circle",
                                view = "lower",
                                show_text = NULL,
                                color_palette = "RdBu") {

  # Default show_text based on shape if not explicitly provided
  if (is.null(show_text)) {
    show_text <- if (shape == "square") TRUE else FALSE
  }

  # 1. Input Validation & Selection
  if (is.null(numeric_vars)) {
    num_data <- data %>% dplyr::select(where(is.numeric))
  } else {
    if (!all(numeric_vars %in% names(data))) stop("Specified numeric variables not found in data.")
    num_data <- data %>% dplyr::select(dplyr::all_of(numeric_vars))
  }

  if (ncol(num_data) < 2) stop("At least two numeric variables are required for correlation analysis.")

  # 2. Compute Correlation and P-value Matrices
  n_vars <- ncol(num_data)
  cor_mat <- stats::cor(num_data, method = method, use = "complete.obs")

  p_mat <- matrix(NA, nrow = n_vars, ncol = n_vars, dimnames = list(names(num_data), names(num_data)))
  for (i in seq_len(n_vars)) {
    for (j in seq_len(n_vars)) {
      if (i != j) {
        test <- stats::cor.test(num_data[[i]], num_data[[j]], method = method)
        p_mat[i, j] <- test$p.value
      } else {
        p_mat[i, j] <- 0
      }
    }
  }

  # 3. Handle Matrix View (Full, Lower, Upper Triangle)
  if (view == "lower") {
    cor_mat[upper.tri(cor_mat)] <- NA
  } else if (view == "upper") {
    cor_mat[lower.tri(cor_mat)] <- NA
  }

  # 4. Reshape into Long Format for ggplot2
  cor_df <- as.data.frame(as.table(cor_mat))
  names(cor_df) <- c("Var1", "Var2", "Correlation")
  cor_df <- tidyr::drop_na(cor_df, Correlation)

  cor_df$Var1 <- factor(cor_df$Var1, levels = names(num_data))
  cor_df$Var2 <- factor(cor_df$Var2, levels = rev(names(num_data)))

  # 5. Build Heatmap Visualization
  p <- ggplot(cor_df, aes(x = Var1, y = Var2))

  if (shape == "circle") {
    p <- p + geom_point(aes(size = abs(Correlation), fill = Correlation), shape = 21, color = "gray40") +
      scale_size_continuous(range = c(4, 14), guide = "none")

    if (show_text) {
      p <- p + geom_text(aes(label = sprintf("%.2f", Correlation)), color = "black", size = 3)
    }

  } else if (shape == "square") {
    p <- p + geom_tile(aes(fill = Correlation), color = "white", linewidth = 0.5)

    if (show_text) {
      p <- p + geom_text(aes(label = sprintf("%.2f", Correlation)), color = "black", size = 3.5)
    }

  } else {
    stop("Invalid shape. Choose 'circle' or 'square'.")
  }

  p <- p +
    theme_minimal() +
    labs(
      x = "",
      y = "",
      title = paste0("Correlation Matrix Heatmap (", stringr::str_to_title(method), ")"),
      fill = "Corr"
    ) +
    scale_fill_distiller(palette = color_palette, limit = c(-1, 1), direction = -1) +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      panel.grid.major = element_blank()
    ) +
    coord_fixed()

  return(list(
    Correlation_Matrix = cor_mat,
    P_Value_Matrix = p_mat,
    Plot = p
  ))
}
