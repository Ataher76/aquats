#' Generalized Linear Models (GLMs) for Fisheries Count and Binary Data
#'
#' This function fits Generalized Linear Models (e.g., Poisson, Quasipoisson, Binomial)
#' to non-normal ecological data, extracts model summaries, incident rate ratios (IRR) or
#' odds ratios, and generates a publication-ready forest plot visualizing effect sizes.
#'
#' @param data A data frame containing the variables.
#' @param response_var Character; the name of the numeric count or binary response variable.
#' @param predictor_vars Character vector; the names of the independent predictor variables.
#' @param family_type Character; distribution family: "poisson", "quasipoisson", or "binomial". Default is "poisson".
#' @param color_palette Character; name of a color palette from RColorBrewer. Default is "Set1".
#'
#' @return A list containing the model summary, coefficient table, and the ggplot object.
#'
#' @import ggplot2 RColorBrewer dplyr
#' @importFrom stats glm as.formula summary.glm poisson quasipoisson binomial coef
#' @importFrom dplyr %>%
#' @export
glm_analysis <- function(data,
                         response_var,
                         predictor_vars,
                         family_type = "poisson",
                         color_palette = "Set1") {

  # 1. Input Validation
  required_cols <- c(response_var, predictor_vars)
  if (!all(required_cols %in% names(data))) stop("Specified columns not found in data.")
  if (!is.numeric(data[[response_var]])) stop(paste(response_var, "must be numeric/integer."))

  # 2. Assign Family Object
  fam_obj <- switch(family_type,
                    "poisson" = stats::poisson(link = "log"),
                    "quasipoisson" = stats::quasipoisson(link = "log"),
                    "binomial" = stats::binomial(link = "logit"),
                    stop("Invalid family_type. Choose 'poisson', 'quasipoisson', or 'binomial'."))

  # 3. Fit GLM Model
  formula_str <- as.formula(paste(response_var, "~", paste(predictor_vars, collapse = " + ")))
  glm_model <- stats::glm(formula_str, family = fam_obj, data = data)
  glm_summary <- summary(glm_model)

  # 4. Extract Coefficients for Forest Plot Visualization
  coef_mat <- stats::coef(glm_summary)
  coef_df <- as.data.frame(coef_mat)
  names(coef_df) <- c("Estimate", "Std_Error", "z_value", "p_value")
  coef_df$Term <- rownames(coef_df)
  coef_df <- coef_df[coef_df$Term != "(Intercept)", ]

  # Calculate exponentiated effect metrics (IRR for Poisson/Quasipoisson, Odds Ratio for Binomial)
  coef_df$Effect_Size <- exp(coef_df$Estimate)
  coef_df$Lower_CI <- exp(coef_df$Estimate - 1.96 * coef_df$Std_Error)
  coef_df$Upper_CI <- exp(coef_df$Estimate + 1.96 * coef_df$Std_Error)

  # 5. Forest Plot Visualization of Model Effects
  p <- ggplot(coef_df, aes(x = stats::reorder(Term, Estimate), y = Effect_Size, color = Term)) +
    geom_pointrange(aes(ymin = Lower_CI, ymax = Upper_CI), size = 1.1) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "gray50", linewidth = 0.8) +
    coord_flip() +
    theme_bw() +
    labs(
      x = "Model Predictors",
      y = if(family_type %in% c("poisson", "quasipoisson")) "Incident Rate Ratio (IRR)" else "Odds Ratio (OR)",
      title = paste0("Generalized Linear Model (", stringr::str_to_title(family_type), " Regression)"),
      subtitle = "Effect Sizes with 95% Confidence Intervals"
    ) +
    scale_color_brewer(palette = color_palette) +
    theme(
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      legend.position = "none"
    )

  return(list(
    Model_Summary = glm_summary,
    Coefficient_Table = coef_df,
    Plot = p
  ))
}
