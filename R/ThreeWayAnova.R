#' Perform Three-Way ANOVA Analysis
#'
#' This function conducts a three-way ANOVA analysis on the given dataset, calculates Tukey's HSD post-hoc tests,
#' estimates marginal means, and generates a plot visualizing the results. Optional jitter and violin plots can be
#' added for enhanced visualization, and results can be automatically saved.
#'
#' @param data A data frame containing four columns: three factors (categorical variables) and one numeric variable.
#' The first three columns must be factors, and the fourth column must be numeric.
#' @param add_jitter Logical; if TRUE, adds jitter points. Default is TRUE.
#' @param jitter_size Numeric; size of jitter points. Default is 2.
#' @param jitter_width Numeric; width of jitter spread within `position_jitterdodge`. Default is 0.2.
#' @param jitter_alpha Numeric; transparency of jitter points (0 to 1). Default is 0.6.
#' @param add_violin Logical; if TRUE, adds violin plots. Default is FALSE.
#' @param violin_alpha Numeric; transparency of violin plots (0 to 1). Default is 0.2.
#' @param violin_scale Numeric; scaling factor for violin width. Default is 1.
#' @param dodge_width Numeric; width for dodging boxplots, jitter, violin, and text. Default is 1.0.
#' @param save_output Logical; if TRUE, saves plot as PDF and other results as text. Default is TRUE.
#' @param output_dir Character; directory path for saving files. Default is current working directory.
#' @param plot_filename Character; base filename for the PDF plot (e.g., "anova_plot"). Default is "three_way_anova_plot".
#' @param text_filename Character; base filename for the text file (e.g., "anova_results"). Default is "three_way_anova_results".
#' @return A list containing:
#' \item{ANOVA_Summary}{Summary of the ANOVA model, including F-values and p-values.}
#' \item{TukeyHSD}{Results of Tukey's HSD post-hoc test for multiple comparisons.}
#' \item{Estimated_Marginal_Means}{Estimated marginal means for factor combinations.}
#' \item{Compact_Letters_Display}{A compact letter display of group differences for easy interpretation.}
#' \item{Plot}{A ggplot object showing the interaction effect and group comparisons, with optional jitter and violin plots.}
#'
#' @details The function performs a three-way ANOVA, tests for significant interactions, and provides estimated
#' marginal means for all factor combinations. It includes a visualization of results, highlighting significant
#' group differences using compact letter displays. Jitter and violin plots can be enabled with customizable
#' parameters, and the `dodge_width` parameter controls the separation between groups. If `save_output` is TRUE,
#' the plot is saved as a PDF and other results as a text file in the specified directory.
#'
#' @examples
#' # Example dataset
#' data <- data.frame(
#'   Factor1 = rep(c("A", "B"), each = 12),
#'   Factor2 = rep(c("X", "Y"), each = 6, times = 2),
#'   Factor3 = rep(c("P", "Q"), times = 12),
#'   NumericVar = rnorm(24)
#' )
#'
#' # Perform Three-Way ANOVA with default settings and save
#' results <- ThreeWayAnova(data, save_output = TRUE, 
#'                         plot_filename = "my_three_way_plot", 
#'                         text_filename = "my_three_way_results")
#' print(results$Plot)
#' 
#' # With jitter, violin, and custom dodge width
#' results_with_jv <- ThreeWayAnova(data, add_jitter = TRUE, add_violin = TRUE, 
#'                                 jitter_size = 1.5, violin_alpha = 0.3, 
#'                                 dodge_width = 1.2, 
#'                                 save_output = TRUE, output_dir = "output_folder")
#' print(results_with_jv$Plot)
#'
#' @import ggplot2 emmeans multcomp
#' @export
#' 

ThreeWayAnova <- function(data, 
                          add_jitter = TRUE, 
                          jitter_size = 2, 
                          jitter_width = 0.2, 
                          jitter_alpha = 0.6,
                          add_violin = FALSE, 
                          violin_alpha = 0.2, 
                          violin_scale = 1,
                          dodge_width = 1.0,  # New parameter for dodging width
                          save_output = TRUE, 
                          output_dir = getwd(), 
                          plot_filename = "three_way_anova_plot", 
                          text_filename = "three_way_anova_results") {
  # Load necessary packages
  if (!requireNamespace("readxl", quietly = TRUE)) install.packages("readxl")
  if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
  if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
  if (!requireNamespace("emmeans", quietly = TRUE)) install.packages("emmeans")
  if (!requireNamespace("multcomp", quietly = TRUE)) install.packages("multcomp")
  
  library(devtools)
  library(readxl)
  library(ggplot2)
  library(emmeans)
  library(multcomp)
  
  # Enhanced data validation with detailed feedback
  if (!is.data.frame(data)) {
    stop("Input 'data' must be a data frame. Current type: ", class(data)[1])
  }
  if (ncol(data) != 4) {
    cat("Data structure:\n")
    print(str(data))
    stop("Data must have exactly four columns: three factors and a numeric variable. Found ", ncol(data), " columns.")
  }
  factor1 <- names(data)[[1]]
  factor2 <- names(data)[[2]]
  factor3 <- names(data)[[3]]
  numeric_var <- names(data)[[4]]
  
  # Convert character columns to factors if not already
  for (col in c(factor1, factor2, factor3)) {
    if (is.character(data[[col]])) {
      data[[col]] <- as.factor(data[[col]])
      cat("Converted", col, "to factor.\n")
    } else if (!is.factor(data[[col]])) {
      cat("Column", col, "is of type", class(data[[col]])[1], ". Expected factor.\n")
      print(str(data[[col]]))
      stop(paste(col, "must be a factor or character convertible to factor."))
    }
  }
  
  # Check that the numeric variable is numeric with detailed feedback
  if (!is.numeric(data[[numeric_var]])) {
    cat("Column", numeric_var, "is of type", class(data[[numeric_var]])[1], ". Expected numeric.\n")
    print(str(data[[numeric_var]]))
    stop(paste(numeric_var, "must be a numeric variable."))
  }
  
  # Run Three-Way ANOVA
  formula <- as.formula(paste(numeric_var, "~", factor1, "*", factor2, "*", factor3))
  aov_model <- aov(formula, data = data)
  aov_summary <- summary(aov_model) # Summary of the ANOVA model
  
  # Tukey's HSD post-hoc test for interaction and main effects
  TukeySHD <- TukeyHSD(aov_model)
  
  # Check for reference grid availability
  tryCatch({
    emmean <- emmeans(aov_model, specs = as.formula(paste("~", factor1, "*", factor2, "*", factor3)))
    emmean_cld <- cld(emmean, Letters = letters)
  }, error = function(e) {
    stop("Error in calculating estimated marginal means: ", e$message)
  }, warning = function(w) {
    cat("Warning during EMM calculation:", w$message, "\n")
  })
  
  # Adjust CLD for plotting
  cld_data <- as.data.frame(emmean_cld)
  cld_data$group <- paste(cld_data[[factor1]], cld_data[[factor2]], cld_data[[factor3]], sep = "_")
  
  # Interaction Plot with Letters (visualizing Factor1 and Factor2, with Factor3 as fill)
  plot <- ggplot(data, aes_string(x = factor1, y = numeric_var, fill = factor2)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7, position = position_dodge(width = dodge_width)) +
    geom_text(data = cld_data, aes_string(x = factor1, y = "emmean", label = ".group", group = factor2),
              position = position_dodge(width = dodge_width), vjust = -0.5, size = 5) +
    theme_bw() +
    labs(x = factor1, y = numeric_var, title = "Three-Way ANOVA Results") +
    theme(text = element_text(size = 12), plot.title = element_text(hjust = 0.5))
  
  # Add jitter if requested
  if (add_jitter) {
    plot <- plot + geom_jitter(aes_string(color = factor2), size = jitter_size, 
                               alpha = jitter_alpha, 
                               position = position_jitterdodge(jitter.width = jitter_width, dodge.width = dodge_width))
  }
  
  # Add violin if requested
  if (add_violin) {
    plot <- plot + geom_violin(aes_string(fill = factor2), alpha = violin_alpha, scale = violin_scale, 
                               position = position_dodge(width = dodge_width), show.legend = FALSE)
  }
  
  # Save outputs if requested
  if (save_output) {
    # Ensure output directory exists
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
      cat("Created output directory:", output_dir, "\n")
    }
    
    # Save plot as PDF
    plot_pdf <- file.path(output_dir, paste0(plot_filename, ".pdf"))
    ggsave(plot_pdf, plot = plot, device = "pdf", width = 12, height = 6, dpi = 300)  # Adjusted width for three factors
    cat("Plot saved as:", plot_pdf, "\n")
    
    # Save other results as text file
    text_file <- file.path(output_dir, paste0(text_filename, ".txt"))
    sink(text_file)
    cat("ANOVA Summary:\n")
    print(aov_summary)
    cat("\nTukey HSD Results:\n")
    print(TukeySHD)
    cat("\nEstimated Marginal Means:\n")
    print(emmean)
    cat("\nCompact Letter Display:\n")
    print(emmean_cld)
    sink()
    cat("Results saved as:", text_file, "\n")
  }
  
  # Return all relevant information as a list
  return(list(
    ANOVA_Summary = aov_summary,
    TukeyHSD = TukeySHD,
    Estimated_Marginal_Means = emmean,
    Compact_Letters_Display = emmean_cld,
    Plot = plot
  ))
}


data <- data.frame(
  Factor1 = rep(c("A", "B"), each = 12),
  Factor2 = rep(c("X", "Y"), each = 6, times = 2),
  Factor3 = rep(c("P", "Q"), times = 12),
  NumericVar = rnorm(24)
)
results <- ThreeWayAnova(data, save_output = TRUE)
print(results$Plot)




