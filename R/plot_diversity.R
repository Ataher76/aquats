#' Visualize and Test Alpha Diversity Patterns
#'
#' Generates publication-ready boxplots for alpha diversity metrics across
#' single or two-factor environmental groupings. Provides optional parametric
#' (ANOVA) or non-parametric (Kruskal-Wallis) hypothesis testing, with optional faceting.
#'
#' @param data A \code{data.frame}, \code{matrix}, or \code{list}.
#' @param group_col Character or integer vector specifying 1 or 2 grouping variables.
#' @param facet_var Optional character specifying a column to facet by.
#' @param index Character. Metric to plot: \code{"shannon"} (default), \code{"simpson"},
#'   \code{"richness"}, \code{"abundance"}, \code{"margalef"}, \code{"menhinick"}, or \code{"pielou"}.
#' @param test Character. Hypothesis test to perform: \code{"anova"} (default),
#'   \code{"kruskal"}, or \code{"none"}.
#' @param color_palette Character. RColorBrewer palette name (default: \code{"Dark2"}).
#' @param permutations Integer. Bootstrap iterations for confidence intervals (default: 1000).
#' @param file Optional character. Export path for summary tables (\code{.xlsx} or \code{.csv}).
#'
#' @return A list with elements \code{Diversity_Data}, \code{Group_Summary},
#'   \code{Test_Result}, and \code{Plot}.
#' @export
#'
#' @import ggplot2
#' @importFrom stats aov as.formula kruskal.test na.omit quantile
#' @importFrom utils write.csv
#' @importFrom tools file_ext
#' @importFrom rlang .data
plot_diversity <- function(
    data,
    group_col = 1,
    facet_var = NULL,
    index = c("shannon", "simpson", "richness", "abundance", "margalef", "menhinick", "pielou"),
    test = c("anova", "kruskal", "none"),
    color_palette = "Dark2",
    permutations = 1000,
    file = NULL
) {

  if (is.matrix(data)) {
    data <- as.data.frame(data)
  } else if (!is.data.frame(data) && !is.list(data)) {
    stop("Input 'data' must be a data.frame, tibble, matrix, or list.")
  }

  index <- match.arg(tolower(index), c("shannon", "simpson", "richness", "abundance", "margalef", "menhinick", "pielou"))
  test  <- match.arg(tolower(test), c("anova", "kruskal", "none"))

  all_groups <- unique(c(group_col, facet_var))

  # 1. Compute Indices
  div_df <- calc_diversity(
    data = data,
    group_col = all_groups
  )

  if (is.character(group_col)) {
    grp_names <- group_col
  } else {
    grp_names <- names(data)[as.integer(group_col)]
  }

  f1 <- grp_names[1]
  f2 <- if (length(grp_names) >= 2) grp_names[2] else NULL

  # 2. Map Target Index
  metric_map <- c(
    "shannon"   = "Shannon",
    "simpson"   = "Simpson",
    "richness"  = "Richness",
    "abundance" = "Abundance",
    "margalef"  = "Margalef",
    "menhinick" = "Menhinick",
    "pielou"    = "Pielou"
  )
  target_col <- metric_map[[index]]
  div_df$Metric_Value <- div_df[[target_col]]
  div_df[[f1]] <- as.factor(div_df[[f1]])
  if (!is.null(f2)) div_df[[f2]] <- as.factor(div_df[[f2]])
  if (!is.null(facet_var)) div_df[[facet_var]] <- as.factor(div_df[[facet_var]])

  # 3. Hypothesis Testing
  test_res <- NULL
  if (test == "anova") {
    test_terms <- if (!is.null(f2)) paste(f1, "*", f2) else f1
    fmla <- stats::as.formula(paste("Metric_Value ~", test_terms))
    test_res <- summary(stats::aov(fmla, data = div_df))
  } else if (test == "kruskal") {
    if (is.null(f2)) {
      test_res <- stats::kruskal.test(stats::as.formula(paste("Metric_Value ~", f1)), data = div_df)
    } else {
      grp_combined <- interaction(div_df[[f1]], div_df[[f2]], sep = " - ")
      test_res <- stats::kruskal.test(div_df$Metric_Value ~ grp_combined)
    }
  }

  # 4. Bootstrap Summaries
  split_factor <- if (is.null(f2)) div_df[[f1]] else interaction(div_df[[f1]], div_df[[f2]], sep = " - ")
  sub_vals <- split(div_df$Metric_Value, split_factor)

  boot_rows <- lapply(names(sub_vals), function(grp) {
    vals <- stats::na.omit(sub_vals[[grp]])
    n <- length(vals)
    if (n < 2) {
      return(data.frame(Group = grp, Index = target_col, N = n, Mean = mean(vals), LCI = NA, UCI = NA))
    }
    b_means <- replicate(permutations, {
      mean(vals[sample.int(n, size = n, replace = TRUE)])
    })
    data.frame(
      Group = grp,
      Index = target_col,
      N = n,
      Mean = round(mean(b_means), 3),
      LCI = round(stats::quantile(b_means, 0.025, names = FALSE), 3),
      UCI = round(stats::quantile(b_means, 0.975, names = FALSE), 3),
      stringsAsFactors = FALSE
    )
  })
  summary_stats <- do.call(rbind, boot_rows)

  # 5. Publication Visualization
  y_labels <- c(
    "Shannon"   = "Shannon Diversity Index (H')",
    "Simpson"   = "Simpson Diversity Index (1 - D)",
    "Richness"  = "Species Richness (S)",
    "Abundance" = "Total Abundance (N)",
    "Margalef"  = "Margalef's Richness (d)",
    "Menhinick" = "Menhinick's Diversity (D_mn)",
    "Pielou"    = "Pielou's Evenness (J')"
  )

  fill_var <- if (!is.null(f2)) f2 else f1

  p <- ggplot2::ggplot(
    div_df,
    ggplot2::aes(x = .data[[f1]], y = .data[["Metric_Value"]], fill = .data[[fill_var]])
  )

  if (is.null(f2)) {
    p <- p +
      ggplot2::geom_boxplot(width = 0.55, alpha = 0.8, outlier.shape = NA, color = "black") +
      ggplot2::geom_jitter(width = 0.18, size = 2.2, alpha = 0.65, shape = 21, fill = "white", color = "black")
  } else {
    p <- p +
      ggplot2::geom_boxplot(position = ggplot2::position_dodge(0.8), width = 0.65, alpha = 0.8, outlier.shape = NA, color = "black") +
      ggplot2::geom_point(position = ggplot2::position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8), shape = 21, fill = "white", color = "black", size = 2, alpha = 0.7)
  }

  p <- p +
    ggplot2::scale_fill_brewer(palette = color_palette) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(
      axis.line = ggplot2::element_line(linewidth = 0.7, color = "black"),
      axis.text = ggplot2::element_text(color = "black", face = "bold"),
      axis.title = ggplot2::element_text(face = "bold"),
      legend.position = if (is.null(f2)) "none" else "top",
      legend.title = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(face = "bold", size = 11, hjust = 0.5),
      strip.background = ggplot2::element_rect(fill = "grey90", color = "black"),
      strip.text = ggplot2::element_text(face = "bold")
    ) +
    ggplot2::labs(
      x = f1,
      y = y_labels[[target_col]],
      fill = fill_var,
      title = paste("Community Diversity:", y_labels[[target_col]])
    )

  if (!is.null(facet_var)) {
    p <- p + ggplot2::facet_wrap(stats::as.formula(paste("~", facet_var)))
  }

  return(list(
    Diversity_Data = div_df[, !names(div_df) %in% c("Metric_Value")],
    Group_Summary = summary_stats,
    Test_Result = test_res,
    Plot = p
  ))
}
