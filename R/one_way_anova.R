#' One-Way Analysis of Variance (ANOVA) with Faceting & Visualization
#'
#' This function performs a one-way ANOVA, calculates Tukey's HSD post-hoc test,
#' extracts Estimated Marginal Means (EMMs), and generates a publication-ready plot
#' with optional faceting by a secondary categorical variable.
#'
#' @param data A data frame in long format.
#' @param factor_var Character; the name of the primary categorical independent variable.
#' @param numeric_var Character; the name of the continuous response variable.
#' @param facet_var Optional character; the name of a secondary categorical variable to facet the plot by (default: \code{NULL}).
#' @param factor_levels Optional character vector; specifies the exact order of the categorical levels on the x-axis.
#' @param plot_type Character; type of plot: "boxplot", "barplot", or "pointrange". Default is "boxplot".
#' @param error_type Character; error bar type for bar/pointrange plots: "se" (Standard Error) or "sd" (Standard Deviation).
#' @param sig_display Character; display significance via "letters" (compact letter display) or "stars" (p-value brackets).
#' @param y_limits Optional numeric vector of length 2; explicitly sets the Y-axis limits.
#' @param add_jitter Logical; if TRUE, adds jittered points to the boxplot.
#' @param show_mean Logical; if TRUE, displays a mean point inside the boxplot.
#' @param mean_color Character; color for the mean point in the boxplot.
#' @param color_palette Character; name of a color palette from RColorBrewer.
#'
#' @return A list containing the ANOVA summary, Tukey HSD results, EMMs, summary table, and ggplot object.
#'
#' @export
#' @import ggplot2 emmeans multcomp RColorBrewer ggpubr
#' @importFrom stats aov as.formula sd TukeyHSD
#' @importFrom dplyr group_by summarise mutate n
#' @importFrom rlang .data
one_way_anova <- function(data,
                          factor_var,
                          numeric_var,
                          facet_var = NULL,
                          factor_levels = NULL,
                          plot_type = "boxplot",
                          error_type = "se",
                          sig_display = "letters",
                          y_limits = NULL,
                          add_jitter = TRUE,
                          show_mean = TRUE,
                          mean_color = "darkred",
                          color_palette = "Set1") {

  # Standard Defensive Ingestion Guard
  if (is.matrix(data)) {
    data <- as.data.frame(data)
  } else if (!is.data.frame(data)) {
    stop("Input 'data' must be a data.frame, tibble, or matrix.")
  }

  # 1. Input Validation & Defensive NA Cleaning
  required_cols <- c(factor_var, numeric_var)
  if (!is.null(facet_var)) required_cols <- c(required_cols, facet_var)

  if (!all(required_cols %in% names(data))) {
    stop("One or more specified columns not found in data.")
  }
  if (!is.numeric(data[[numeric_var]])) stop(paste(numeric_var, "must be numeric."))

  # Isolate valid rows
  valid_rowsstats <- stats::complete.cases(data[, required_cols, drop = FALSE])
  clean_data <- data[valid_rowsstats, , drop = FALSE]
  na_dropped <- nrow(data) - nrow(clean_data)

  if (na_dropped > 0) {
    message(sprintf("Note: Automatically removed %d row(s) containing NA values in target variables.", na_dropped))
  }

  if (nrow(clean_data) < 5) {
    stop("ANOVA requires at least 5 complete observations.")
  }

  # Factor Ordering / Coercion
  if (!is.null(factor_levels)) {
    clean_data[[factor_var]] <- factor(clean_data[[factor_var]], levels = factor_levels)
  } else {
    clean_data[[factor_var]] <- as.factor(clean_data[[factor_var]])
  }

  if (!is.null(facet_var)) {
    clean_data[[facet_var]] <- as.factor(clean_data[[facet_var]])
  }

  # 2. Statistical Models
  formula_str <- as.formula(paste(numeric_var, "~", factor_var))
  aov_model <- stats::aov(formula_str, data = clean_data)
  aov_summary <- summary(aov_model)
  tukey_res <- stats::TukeyHSD(aov_model)

  # 3. Estimated Marginal Means & Letters
  emmean <- emmeans::emmeans(aov_model, specs = as.formula(paste("~", factor_var)))
  emmean_cld <- multcomp::cld(emmean, Letters = letters)
  emmean_cld$.group <- trimws(emmean_cld$.group)

  # 4. Summary Statistics for Error Bars & Letter Placement
  group_vars <- if (!is.null(facet_var)) c(factor_var, facet_var) else factor_var

  sum_data <- clean_data %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) %>%
    dplyr::summarise(
      mean_val = mean(.data[[numeric_var]], na.rm = TRUE),
      sd_val = stats::sd(.data[[numeric_var]], na.rm = TRUE),
      se_val = sd_val / sqrt(dplyr::n()),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      error_margin = if (error_type == "sd") sd_val else se_val,
      upper_limit = mean_val + error_margin
    )

  # Merge letters if no facet; if faceted, letters apply per facet or main factor
  if (is.null(facet_var)) {
    emmean_cld <- merge(emmean_cld, sum_data, by = factor_var)
  }

  # 5. Base Plot Initialization
  p <- ggplot2::ggplot(
    clean_data,
    ggplot2::aes(x = .data[[factor_var]], y = .data[[numeric_var]], fill = .data[[factor_var]])
  ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::labs(x = factor_var, y = numeric_var, title = "One-Way ANOVA") +
    ggplot2::scale_fill_brewer(palette = color_palette) +
    ggplot2::theme(
      legend.position = "none",
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(color = "black", linewidth = 1),
      axis.title = ggplot2::element_text(face = "bold"),
      axis.text = ggplot2::element_text(color = "black", face = "bold"),
      strip.background = ggplot2::element_rect(fill = "grey90", color = "black"),
      strip.text = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(face = "bold", size = 12, hjust = 0.5)
    )

  # 6. Plot Type Routing
  if (plot_type == "boxplot") {
    p <- p +
      ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.7)

    if (add_jitter) {
      p <- p + ggplot2::geom_jitter(color = "black", size = 1.5, width = 0.2, alpha = 0.5)
    }
    if (show_mean) {
      p <- p + ggplot2::stat_summary(fun = mean, geom = "point", shape = 18, size = 4, color = mean_color)
    }

  } else if (plot_type == "barplot") {
    p <- ggplot2::ggplot(
      sum_data,
      ggplot2::aes(x = .data[[factor_var]], y = .data$mean_val, fill = .data[[factor_var]])
    ) +
      ggplot2::geom_col(alpha = 0.8, color = "black", width = 0.65) +
      ggplot2::geom_errorbar(ggplot2::aes(ymin = .data$mean_val - .data$error_margin, ymax = .data$upper_limit), width = 0.2) +
      ggplot2::theme_bw(base_size = 12) +
      ggplot2::labs(x = factor_var, y = numeric_var, title = paste("One-Way ANOVA (Mean \u00B1", toupper(error_type), ")")) +
      ggplot2::scale_fill_brewer(palette = color_palette) +
      ggplot2::theme(
        legend.position = "none",
        panel.grid.major = ggplot2::element_blank(),
        panel.border = ggplot2::element_rect(color = "black", linewidth = 1),
        axis.title = ggplot2::element_text(face = "bold"),
        axis.text = ggplot2::element_text(color = "black", face = "bold"),
        strip.background = ggplot2::element_rect(fill = "grey90", color = "black"),
        strip.text = ggplot2::element_text(face = "bold"),
        plot.title = ggplot2::element_text(face = "bold", size = 12, hjust = 0.5)
      )

  } else if (plot_type == "pointrange") {
    p <- ggplot2::ggplot(
      sum_data,
      ggplot2::aes(x = .data[[factor_var]], y = .data$mean_val, color = .data[[factor_var]])
    ) +
      ggplot2::geom_point(size = 3) +
      ggplot2::geom_errorbar(ggplot2::aes(ymin = .data$mean_val - .data$error_margin, ymax = .data$upper_limit), width = 0.2) +
      ggplot2::theme_bw(base_size = 12) +
      ggplot2::labs(x = factor_var, y = numeric_var, title = paste("One-Way ANOVA (Mean \u00B1", toupper(error_type), ")")) +
      ggplot2::scale_color_brewer(palette = color_palette) +
      ggplot2::theme(
        legend.position = "none",
        panel.grid.major = ggplot2::element_blank(),
        panel.border = ggplot2::element_rect(color = "black", linewidth = 1),
        axis.title = ggplot2::element_text(face = "bold"),
        axis.text = ggplot2::element_text(color = "black", face = "bold"),
        strip.background = ggplot2::element_rect(fill = "grey90", color = "black"),
        strip.text = ggplot2::element_text(face = "bold"),
        plot.title = ggplot2::element_text(face = "bold", size = 12, hjust = 0.5)
      )
  }

  # 7. Apply Facet Wrap if requested
  if (!is.null(facet_var)) {
    p <- p + ggplot2::facet_wrap(stats::as.formula(paste("~", facet_var)), scales = "free_y")
  }

  # 8. Significance Display (Letters)
  if (sig_display == "letters" && is.null(facet_var)) {
    y_nudge <- max(clean_data[[numeric_var]], na.rm = TRUE) * 0.05
    if (plot_type == "boxplot") {
      max_y <- max(clean_data[[numeric_var]], na.rm = TRUE)
      p <- p + ggplot2::geom_text(data = emmean_cld, ggplot2::aes(x = .data[[factor_var]], y = max_y + y_nudge, label = .data$.group), inherit.aes = FALSE, size = 4.5, fontface = "bold")
    } else {
      p <- p + ggplot2::geom_text(data = emmean_cld, ggplot2::aes(x = .data[[factor_var]], y = .data$upper_limit + y_nudge, label = .data$.group), inherit.aes = FALSE, size = 4.5, fontface = "bold")
    }
  } else if (sig_display == "stars") {
    p <- p + ggpubr::geom_pwc(
      data = clean_data,
      ggplot2::aes(x = .data[[factor_var]], y = .data[[numeric_var]]),
      method = "tukey_hsd",
      label = "p.adj.signif",
      hide.ns = TRUE,
      inherit.aes = FALSE,
      step.increase = 0.1,
      vjust = 0.5
    )
  }

  # 9. Apply Custom Y-Axis Limits (if provided)
  if (!is.null(y_limits)) {
    if (length(y_limits) == 2 && is.numeric(y_limits)) {
      p <- p + ggplot2::coord_cartesian(ylim = y_limits)
    } else {
      warning("y_limits must be a numeric vector of length 2. Ignoring custom limits.")
    }
  }

  return(list(
    ANOVA_Summary = aov_summary,
    TukeyHSD = tukey_res,
    EMMs = emmean,
    Summary_Table = sum_data,
    Plot = p
  ))
}
