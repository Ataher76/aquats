#' Flexible Regression Analysis with Shaded Confidence Intervals and Statistics
#'
#' This function fits linear, logarithmic, or polynomial regression models to data,
#' calculates coefficients, R-squared, and p-values, and generates a publication-ready
#' plot with confidence interval shading and equation annotations.
#'
#' @param data A data frame containing the variables.
#' @param x_var Character; the name of the independent variable (x-axis).
#' @param y_var Character; the name of the dependent response variable (y-axis).
#' @param fit_type Character; type of model: "linear", "logarithmic", or "polynomial". Default is "linear".
#' @param group_var Optional character; the name of a categorical variable to group/color regressions by.
#' @param color_palette Character; name of a color palette from RColorBrewer. Default is "Set1".
#'
#' @return A list containing model summaries, regression parameters, and the ggplot object.
#'
#' @import ggplot2 RColorBrewer dplyr
#' @importFrom stats lm summary.lm as.formula poly
#' @importFrom dplyr %>%
#' @export
regression_analysis <- function(data,
                                x_var,
                                y_var,
                                fit_type = "linear",
                                group_var = NULL,
                                color_palette = "Set1") {

  # 1. Input Validation
  required_cols <- c(x_var, y_var)
  if (!is.null(group_var)) required_cols <- c(required_cols, group_var)
  if (!all(required_cols %in% names(data))) stop("Specified columns not found in data.")
  if (!is.numeric(data[[x_var]])) stop(paste(x_var, "must be numeric."))
  if (!is.numeric(data[[y_var]])) stop(paste(y_var, "must be numeric."))

  orig_x <- x_var
  orig_y <- y_var

  # 2. Formula Construction for Model Fitting (lm)
  if (fit_type == "linear") {
    mod_formula <- if (is.null(group_var)) as.formula(paste(y_var, "~", x_var)) else as.formula(paste(y_var, "~", x_var, "*", group_var))
  } else if (fit_type == "logarithmic") {
    data$log_x <- log(data[[x_var]])
    data$log_y <- log(data[[y_var]])
    x_var <- "log_x"
    y_var <- "log_y"
    mod_formula <- if (is.null(group_var)) as.formula("log_y ~ log_x") else as.formula(paste(y_var, "~", x_var, "*", group_var))
  } else if (fit_type == "polynomial") {
    mod_formula <- if (is.null(group_var)) as.formula(paste(y_var, "~ poly(", x_var, ", 2)")) else as.formula(paste(y_var, "~ poly(", x_var, ", 2) *", group_var))
  } else {
    stop("Invalid fit_type. Choose 'linear', 'logarithmic', or 'polynomial'.")
  }

  # 3. Fit Model & Extract Stats
  model <- stats::lm(mod_formula, data = data)
  model_summary <- summary(model)

  r_squared <- model_summary$r.squared
  f_stat <- model_summary$fstatistic
  p_val <- if (!is.null(f_stat)) stats::pf(f_stat[1], f_stat[2], f_stat[3], lower.tail = FALSE) else NA

  # 4. Base Plot Initialization (Using standard aesthetic mapping references x and y)
  smooth_formula <- if (fit_type == "polynomial") y ~ poly(x, 2) else y ~ x

  if (is.null(group_var)) {
    p <- ggplot(data, aes(x = .data[[x_var]], y = .data[[y_var]])) +
      geom_point(alpha = 0.6, size = 2, color = "darkblue") +
      geom_smooth(method = "lm", formula = smooth_formula,
                  se = TRUE, color = "red", fill = "pink", alpha = 0.3, linewidth = 1)
  } else {
    data[[group_var]] <- as.factor(data[[group_var]])
    p <- ggplot(data, aes(x = .data[[x_var]], y = .data[[y_var]], color = .data[[group_var]], fill = .data[[group_var]])) +
      geom_point(alpha = 0.6, size = 2) +
      geom_smooth(method = "lm", formula = smooth_formula,
                  se = TRUE, alpha = 0.2, linewidth = 1) +
      scale_color_brewer(palette = color_palette) +
      scale_fill_brewer(palette = color_palette)
  }

  # 5. Clean Axis Labels
  clean_x <- gsub("_", " ", sub("_cm|_g", "", orig_x))
  clean_y <- gsub("_", " ", sub("_cm|_g", "", orig_y))

  axis_x_label <- if (fit_type == "logarithmic") paste0("ln(", clean_x, ")") else clean_x
  axis_y_label <- if (fit_type == "logarithmic") paste0("ln(", clean_y, ")") else clean_y

  p <- p + theme_bw() +
    labs(
      x = axis_x_label,
      y = axis_y_label,
      title = paste0("Regression Analysis (", stringr::str_to_title(fit_type), " Fit)"),
      subtitle = sprintf("R2 = %.3f | p-value = %.4g", r_squared, p_val),
      color = group_var,
      fill = group_var
    ) +
    theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))

  return(list(
    Model_Summary = model_summary,
    R_Squared = r_squared,
    P_Value = p_val,
    Plot = p
  ))
}

