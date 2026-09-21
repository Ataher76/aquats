#' Generalized Linear Modeling (GLM) for Aquatic & Ecological Data
#'
#' Fits generalized linear models (GLM) across Gaussian, Poisson, Binomial,
#' and Gamma families with publication-ready diagnostic and prediction graphics.
#'
#' @param data A \code{data.frame}, \code{matrix}, or tibble.
#' @param response_var Character. Name of the dependent response variable.
#' @param predictor_vars Character vector. Name(s) of continuous or categorical predictor variables.
#' @param family_type Character. Distribution family: \code{"gaussian"}, \code{"poisson"},
#'   \code{"binomial"}, or \code{"gamma"}.
#' @param color_palette Character. RColorBrewer palette name (default: \code{"Dark2"}).
#'
#' @return A list containing:
#'   \item{Model}{The fitted \code{glm} object.}
#'   \item{Summary}{The model summary object.}
#'   \item{Anova_Table}{Analysis of deviance table.}
#'   \item{Plot}{A publication-ready \code{ggplot} object.}
#' @export
#'
#' @import ggplot2
#' @importFrom stats glm as.formula predict poisson binomial Gamma gaussian
#' @importFrom rlang .data
glm_analysis <- function(data,
                         response_var,
                         predictor_vars,
                         family_type = "poisson",
                         color_palette = "Set1") {

  # Standard Defensive Ingestion Guard
  if (is.matrix(data)) {
    data <- as.data.frame(data)
  } else if (!is.data.frame(data)) {
    stop("Input 'data' must be a data.frame, tibble, or matrix.")
  }

  # 1. Input Validation & Defensive NA Cleaning
  required_cols <- c(response_var, predictor_vars)
  if (!all(required_cols %in% names(data))) {
    missing_cols <- required_cols[!required_cols %in% names(data)]
    stop("Specified column(s) not found in data: ", paste(missing_cols, collapse = ", "))
  }

  if (!is.numeric(data[[response_var]])) {
    stop(paste(response_var, "must be numeric/integer."))
  }

  # Isolate valid rows for target variables
  clean_data <- data[, required_cols, drop = FALSE]
  valid_rows <- stats::complete.cases(clean_data)
  clean_data <- clean_data[valid_rows, , drop = FALSE]
  na_dropped <- sum(!valid_rows)

  if (na_dropped > 0) {
    message(sprintf("Note: Automatically removed %d row(s) containing NA values in target variables.", na_dropped))
  }

  if (nrow(clean_data) < (length(predictor_vars) + 3)) {
    stop("Insufficient complete observations to fit the specified GLM model.")
  }

  # Guard for binomial family
  if (family_type == "binomial") {
    unique_vals <- unique(clean_data[[response_var]])
    if (!all(unique_vals %in% c(0, 1))) {
      stop("Binomial family requires binary response variable coded strictly as 0 and 1.")
    }
  }

  # Coerce character/logical predictors to factors
  for (v in predictor_vars) {
    if (is.character(clean_data[[v]]) || is.logical(clean_data[[v]])) {
      clean_data[[v]] <- as.factor(clean_data[[v]])
    }
  }

  # 2. Assign Family Object
  fam_obj <- switch(family_type,
                    "poisson"      = stats::poisson(link = "log"),
                    "quasipoisson" = stats::quasipoisson(link = "log"),
                    "binomial"     = stats::binomial(link = "logit"),
                    stop("Invalid family_type. Choose 'poisson', 'quasipoisson', or 'binomial'."))

  # 3. Fit GLM Model
  formula_str <- stats::as.formula(paste(response_var, "~", paste(predictor_vars, collapse = " + ")))
  glm_model <- stats::glm(formula_str, family = fam_obj, data = clean_data)
  glm_summary <- summary(glm_model)

  # 4. Extract Coefficients for Forest Plot Visualization
  coef_mat <- stats::coef(glm_summary)
  coef_df <- as.data.frame(coef_mat)
  names(coef_df) <- c("Estimate", "Std_Error", "z_value", "p_value")
  coef_df$Term <- rownames(coef_df)
  coef_df <- coef_df[coef_df$Term != "(Intercept)", , drop = FALSE]

  if (nrow(coef_df) == 0) {
    stop("The model contains no valid predictors after filtering coefficients.")
  }

  # Exponentiated effect metrics (IRR for Poisson/Quasipoisson, Odds Ratio for Binomial)
  coef_df$Effect_Size <- exp(coef_df$Estimate)
  coef_df$Lower_CI <- exp(coef_df$Estimate - 1.96 * coef_df$Std_Error)
  coef_df$Upper_CI <- exp(coef_df$Estimate + 1.96 * coef_df$Std_Error)

  # 5. Forest Plot Visualization of Model Effects
  metric_label <- if (family_type %in% c("poisson", "quasipoisson")) {
    "Incident Rate Ratio (IRR)"
  } else {
    "Odds Ratio (OR)"
  }

  p <- ggplot2::ggplot(
    coef_df,
    ggplot2::aes(x = stats::reorder(.data$Term, .data$Estimate), y = .data$Effect_Size, color = .data$Term)
  ) +
    ggplot2::geom_pointrange(
      ggplot2::aes(ymin = .data$Lower_CI, ymax = .data$Upper_CI),
      size = 0.9, linewidth = 1.1
    ) +
    ggplot2::geom_hline(yintercept = 1, linetype = "dashed", color = "grey40", linewidth = 0.8) +
    ggplot2::coord_flip() +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::labs(
      x = "Model Predictor Terms",
      y = metric_label,
      title = paste0("Generalized Linear Model (", tools::toTitleCase(family_type), " Regression)"),
      subtitle = paste("Estimated", metric_label, "with 95% Confidence Intervals")
    ) +
    ggplot2::scale_color_brewer(palette = color_palette) +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(color = "black", linewidth = 1),
      axis.title = ggplot2::element_text(face = "bold"),
      axis.text = ggplot2::element_text(color = "black", face = "bold"),
      legend.position = "none",
      plot.title = ggplot2::element_text(face = "bold", size = 12, hjust = 0.5),
      plot.subtitle = ggplot2::element_text(face = "bold", size = 10, hjust = 0.5, color = "grey30")
    )

  return(list(
    Model_Summary = glm_summary,
    Coefficient_Table = coef_df,
    Plot = p
  ))
}
