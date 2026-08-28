#' OneWayAnova: Conduct One-Way Analysis of Variance (ANOVA) with Visualization
#'
#' This function performs one-way ANOVA for a given dataset with one categorical 
#' factor and one numeric response variable. It computes the ANOVA table, 
#' performs Tukey's Honest Significant Difference (HSD) post-hoc test, and 
#' visualizes the results with a boxplot annotated with significant group letters.
#' Optional jitter and violin plots can be added for enhanced visualization.
#' Results are automatically saved as a PDF (plot) and text file (other outputs).
#'
#' @param data A data frame with exactly two columns: 
#'   - The first column must be a factor representing the categorical independent variable.
#'   - The second column must be numeric representing the dependent variable.
#' @param add_jitter Logical; if TRUE, adds jitter points. Default is TRUE.
#' @param jitter_size Numeric; size of jitter points. Default is 2.
#' @param jitter_width Numeric; width of jitter spread. Default is 0.2.
#' @param jitter_alpha Numeric; transparency of jitter points (0 to 1). Default is 0.6.
#' @param add_violin Logical; if TRUE, adds violin plots. Default is FALSE.
#' @param violin_alpha Numeric; transparency of violin plots (0 to 1). Default is 0.2.
#' @param violin_scale Numeric; scaling factor for violin width. Default is 1.
#' @param color_palette Character; name of a color palette from RColorBrewer or a vector of colors. Default is "Set1".
#' @param save_output Logical; if TRUE, saves plot as PDF and other results as text. Default is TRUE.
#' @param output_dir Character; directory path for saving files. Default is current working directory.
#' @param plot_filename Character; base filename for the PDF plot (e.g., "anova_plot"). Default is "anova_plot".
#' @param text_filename Character; base filename for the text file (e.g., "anova_results"). Default is "anova_results".
#' @return A list containing the following elements:
#'   - \code{ANOVA_Summary}: A summary table of the ANOVA model.
#'   - \code{TukeyHSD}: Results from Tukey's HSD post-hoc test.
#'   - \code{Estimated_Marginal_Means}: Estimated marginal means (EMMs).
#'   - \code{Compact_Letters_Display}: Compact letter display for group comparisons.
#'   - \code{Plot}: A \code{ggplot2} boxplot of the results with annotated group letters, 
#'     optionally including jitter and violin plots.
#'
#' @details 
#' The function checks the input dataset for compatibility and performs ANOVA using 
#' the \code{aov()} function. Post-hoc analysis is conducted with Tukey's HSD and 
#' EMMs are estimated for visualization. The boxplot highlights group differences 
#' with distinct letter labels. Jitter and violin plots can be enabled with customizable 
#' parameters, and their colors are mapped to the factor levels using the specified 
#' \code{color_palette}. If \code{save_output} is TRUE, the plot is saved as a PDF and 
#' other results as a text file in the specified directory.
#'
#' @examples 
#' # Example dataset
#' example_data <- data.frame(
#'   Treatment = rep(c("A", "B", "C"), each = 10),
#'   Response = c(rnorm(10, 5, 1), rnorm(10, 6, 1), rnorm(10, 7, 1))
#' )
#' 
#' # Perform One-Way ANOVA with default settings and save
#' result <- OneWayAnova(example_data, save_output = TRUE, 
#'                       plot_filename = "my_anova_plot", 
#'                       text_filename = "my_anova_results")
#' print(result$Plot)
#' 
#' # With jitter, violin, and custom color palette
#' result_with_jv <- OneWayAnova(example_data, add_jitter = TRUE, add_violin = TRUE, 
#'                              jitter_size = 1.5, violin_alpha = 0.3, 
#'                              color_palette = "Dark2", 
#'                              save_output = TRUE, output_dir = "output_folder")
#' print(result_with_jv$Plot)
#'
#' @note Ensure that the input data has the required structure: one factor column 
#' and one numeric column. Use \code{as.factor()} to convert categorical variables if needed.
#' The output directory must exist or be creatable; otherwise, an error will occur.
#' The \code{color_palette} can be any RColorBrewer palette name (e.g., "Set1", "Dark2") 
#' or a vector of color names/hex codes.
#' 
#' @import ggplot2 emmeans multcomp RColorBrewer
#' @export

OneWayAnova <- function(data, 
                        add_jitter = TRUE, 
                        jitter_size = 2, 
                        jitter_width = 0.2, 
                        jitter_alpha = 0.6,
                        add_violin = FALSE, 
                        violin_alpha = 0.2, 
                        violin_scale = 1,
                        color_palette = "Set1",  # New parameter for color customization
                        save_output = TRUE, 
                        output_dir = getwd(), 
                        plot_filename = "anova_plot", 
                        text_filename = "anova_results") {
  # Load necessary packages
  if (!requireNamespace("readxl", quietly = TRUE)) install.packages("readxl")
  if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
  if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
  if (!requireNamespace("emmeans", quietly = TRUE)) install.packages("emmeans")
  if (!requireNamespace("multcomp", quietly = TRUE)) install.packages("multcomp")
  if (!requireNamespace("RColorBrewer", quietly = TRUE)) install.packages("RColorBrewer")
  
  library(devtools)
  library(readxl)
  library(ggplot2)
  library(emmeans)
  library(multcomp)
  library(RColorBrewer)
  
  # Check if data has the correct structure
  if (ncol(data) != 2) {
    stop("Data should have exactly two columns: a factor and a numeric variable.")
  }
  
  # Define column names for easy access
  factor_var <- names(data)[[1]]
  numeric_var <- names(data)[[2]]
  
  # Ensure the first column is a factor
  data[[factor_var]] <- as.factor(data[[factor_var]])
  
  # Check that the numeric variable is numeric
  if (!is.numeric(data[[numeric_var]])) {
    stop(paste(numeric_var, "must be a numeric variable."))
  }
  
  # Run ANOVA
  aov_model <- aov(as.formula(paste(numeric_var, "~", factor_var)), data = data) 
  aov_summary <- summary(aov_model) # Summary of the ANOVA model
  
  # Tukey's HSD post-hoc test
  TukeySHD <- TukeyHSD(aov_model)
  
  # Check for reference grid availability
  tryCatch({
    emmean <- emmeans(aov_model, specs = as.formula(paste("~", factor_var)))
    emmean_cld <- cld(emmean, Letters = letters)
  }, error = function(e) {
    stop("Error in calculating estimated marginal means: ", e$message)
  })
  
  # Plotting ANOVA results
  plot <- ggplot(data, aes_string(x = factor_var, y = numeric_var, fill = factor_var)) +
    geom_boxplot(show.legend = FALSE, outlier.shape = NA, alpha = 0.7) +
    theme_bw() +
    labs(x = factor_var, y = numeric_var, title = "One-Way ANOVA with Tukey's HSD") +
    geom_text(data = emmean_cld, aes_string(x = factor_var, y = "upper.CL", label = ".group"),
              vjust = -0.5, size = 5) +
    scale_fill_brewer(palette = color_palette) +  # Apply color palette to fill
    theme(text = element_text(size = 12), plot.title = element_text(hjust = 0.5))
  
  # Add jitter if requested
  if (add_jitter) {
    plot <- plot + geom_jitter(aes_string(color = factor_var), size = jitter_size, width = jitter_width, 
                               alpha = jitter_alpha, show.legend = TRUE) +
      scale_color_brewer(palette = color_palette)  # Apply color palette to jitter
  }
  
  # Add violin if requested
  if (add_violin) {
    plot <- plot + geom_violin(aes(fill = factor_var), alpha = violin_alpha, scale = violin_scale, 
                               show.legend = TRUE) +
      scale_fill_brewer(palette = color_palette)  # Ensure fill matches palette
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
    ggsave(plot_pdf, plot = plot, device = "pdf", width = 8, height = 6, dpi = 300)
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
