#' Parametric and Non-Linear Regression Analysis with Fitted Equations
#'
#' Fits linear, allometric power, exponential, logarithmic, quadratic polynomial,
#' or cubic regression models. Automatically formats algebraic mathematical formulas,
#' R-squared values, and p-values, displays confidence or prediction intervals, and renders
#' publication-ready ggplot2 graphics.
#'
#' @param data A \code{data.frame}, \code{matrix}, or tibble.
#' @param x_var Character. Name of the independent predictor variable.
#' @param y_var Character. Name of the dependent response variable.
#' @param fit_type Character. Model type: \code{"linear"} (default), \code{"power"} (allometric),
#'   \code{"exponential"}, \code{"logarithmic"}, \code{"polynomial"} (quadratic), or \code{"cubic"}.
#' @param group_var Optional character. Grouping factor for multi-category regression.
#' @param facet Logical. If \code{TRUE} and \code{group_var} is provided, facets by group. Default is \code{FALSE}.
#' @param interval Character. Uncertainty band type: \code{"confidence"} (default),
#'   \code{"prediction"}, or \code{"none"}.
#' @param show_equation Logical. If \code{TRUE} (default), prints the regression equation,
#'   \eqn{R^2}, and p-value on the canvas.
#' @param color_palette Character. RColorBrewer palette name (default: \code{"Set1"}).
#' @param point_size Numeric. Size of scatter points (default: 2.2).
#' @param point_alpha Numeric. Opacity of scatter points (default: 0.75).
#' @param title Optional character. Plot title.
#'
#' @return A list containing:
#'   \item{Plot}{Publication-grade \code{ggplot} object.}
#'   \item{Models}{Named list of fitted \code{lm} model objects per group.}
#'   \item{Equation_Table}{Data frame of algebraic formulas, \eqn{R^2}, and p-values.}
#'   \item{Predictions}{Data frame of calculated predictions with confidence/prediction limits.}
#' @export
#'
#' @import ggplot2
#' @importFrom stats lm predict as.formula complete.cases pf poly
#' @importFrom tools toTitleCase
#' @importFrom rlang .data
regression_analysis <- function(data,
                                x_var,
                                y_var,
                                fit_type = c("linear", "power", "exponential", "logarithmic", "polynomial", "cubic"),
                                group_var = NULL,
                                facet = FALSE,
                                interval = c("confidence", "prediction", "none"),
                                show_equation = TRUE,
                                color_palette = "Set1",
                                point_size = 2.2,
                                point_alpha = 0.75,
                                title = NULL) {

  fit_type <- match.arg(fit_type)
  interval <- match.arg(interval)

  # 1. Defensive Ingestion Guard
  if (is.matrix(data)) {
    data <- as.data.frame(data)
  } else if (!is.data.frame(data)) {
    stop("Input 'data' must be a data.frame, tibble, or matrix.")
  }

  required_cols <- c(x_var, y_var, group_var)
  missing_cols <- required_cols[!required_cols %in% names(data)]
  if (length(missing_cols) > 0) {
    stop(sprintf("Specified column(s) not found in data: %s", paste(missing_cols, collapse = ", ")))
  }

  if (!is.numeric(data[[x_var]]) || !is.numeric(data[[y_var]])) {
    stop("Both 'x_var' and 'y_var' must be numeric continuous variables.")
  }

  clean_data <- stats::na.omit(data[, required_cols, drop = FALSE])

  # Positive value guards for logarithmic transformations
  if (fit_type %in% c("power", "logarithmic")) {
    invalid_rows <- clean_data[[x_var]] <= 0 | clean_data[[y_var]] <= 0
    if (any(invalid_rows)) {
      clean_data <- clean_data[!invalid_rows, , drop = FALSE]
      message(sprintf("Note: Removed %d non-positive observation(s) required for log-transformation.", sum(invalid_rows)))
    }
  } else if (fit_type == "exponential") {
    invalid_rows <- clean_data[[y_var]] <= 0
    if (any(invalid_rows)) {
      clean_data <- clean_data[!invalid_rows, , drop = FALSE]
      message(sprintf("Note: Removed %d non-positive y observation(s) required for exponential log-link.", sum(invalid_rows)))
    }
  }

  min_req <- if (fit_type == "cubic") 5 else 4
  if (nrow(clean_data) < min_req) {
    stop(sprintf("Insufficient complete observations to fit %s regression model.", fit_type))
  }

  if (!is.null(group_var)) {
    clean_data[[group_var]] <- as.factor(clean_data[[group_var]])
  }

  # 2. Fit Models & Compute Predictions per Group
  split_factor <- if (is.null(group_var)) factor(rep("All", nrow(clean_data))) else clean_data[[group_var]]
  split_data <- split(clean_data, split_factor)

  model_list <- list()
  eq_rows <- list()
  pred_list <- list()

  for (grp_lvl in names(split_data)) {
    sub_df <- split_data[[grp_lvl]]

    if (nrow(sub_df) < min_req) {
      warning(sprintf("Group '%s' contains fewer than %d samples. Skipping fit.", grp_lvl, min_req))
      next
    }

    df_fit <- data.frame(
      x = sub_df[[x_var]],
      y = sub_df[[y_var]]
    )

    if (fit_type == "linear") {
      mod <- stats::lm(y ~ x, data = df_fit)
      b0 <- stats::coef(mod)[1]
      b1 <- stats::coef(mod)[2]
      sign_b1 <- ifelse(b1 >= 0, "+", "-")
      eq_text <- sprintf("y == '%.2f' %s '%.2f'*x", b0, sign_b1, abs(b1))

    } else if (fit_type == "power") {
      df_fit$log_x <- log(df_fit$x)
      df_fit$log_y <- log(df_fit$y)
      mod <- stats::lm(log_y ~ log_x, data = df_fit)
      ln_a <- stats::coef(mod)[1]
      b <- stats::coef(mod)[2]
      a <- exp(ln_a)
      eq_text <- sprintf("y == '%.3f'*x^'%.2f'", a, b)

    } else if (fit_type == "exponential") {
      df_fit$log_y <- log(df_fit$y)
      mod <- stats::lm(log_y ~ x, data = df_fit)
      ln_a <- stats::coef(mod)[1]
      k <- stats::coef(mod)[2]
      a <- exp(ln_a)
      sign_k <- ifelse(k >= 0, "+", "-")
      eq_text <- sprintf("y == '%.3f'*e^(%s'%.3f'*x)", a, sign_k, abs(k))

    } else if (fit_type == "logarithmic") {
      df_fit$log_x <- log(df_fit$x)
      mod <- stats::lm(y ~ log_x, data = df_fit)
      b0 <- stats::coef(mod)[1]
      b1 <- stats::coef(mod)[2]
      sign_b1 <- ifelse(b1 >= 0, "+", "-")
      eq_text <- sprintf("y == '%.2f' %s '%.2f'*ln(x)", b0, sign_b1, abs(b1))

    } else if (fit_type == "polynomial") {
      df_fit$x2 <- df_fit$x^2
      mod <- stats::lm(y ~ x + x2, data = df_fit)
      b0 <- stats::coef(mod)[1]
      b1 <- stats::coef(mod)[2]
      b2 <- stats::coef(mod)[3]
      s1 <- ifelse(b1 >= 0, "+", "-")
      s2 <- ifelse(b2 >= 0, "+", "-")
      eq_text <- sprintf("y == '%.2f' %s '%.2f'*x %s '%.3f'*x^2", b0, s1, abs(b1), s2, abs(b2))

    } else if (fit_type == "cubic") {
      df_fit$x2 <- df_fit$x^2
      df_fit$x3 <- df_fit$x^3
      mod <- stats::lm(y ~ x + x2 + x3, data = df_fit)
      b0 <- stats::coef(mod)[1]
      b1 <- stats::coef(mod)[2]
      b2 <- stats::coef(mod)[3]
      b3 <- stats::coef(mod)[4]
      s1 <- ifelse(b1 >= 0, "+", "-")
      s2 <- ifelse(b2 >= 0, "+", "-")
      s3 <- ifelse(b3 >= 0, "+", "-")
      eq_text <- sprintf("y == '%.2f' %s '%.2f'*x %s '%.3f'*x^2 %s '%.4f'*x^3", b0, s1, abs(b1), s2, abs(b2), s3, abs(b3))
    }

    # Model metrics
    mod_sum <- summary(mod)
    r2 <- mod_sum$r.squared
    f_stat <- mod_sum$fstatistic
    p_val <- if (!is.null(f_stat)) {
      stats::pf(f_stat[1], f_stat[2], f_stat[3], lower.tail = FALSE)
    } else {
      NA_real_
    }

    p_formatted <- if (is.na(p_val)) {
      "p == NA"
    } else if (p_val < 0.001) {
      "p < 0.001"
    } else {
      sprintf("p == '%.3f'", p_val)
    }

    full_eq_label <- sprintf("atop(%s, paste(R^2 == '%.2f', ', ', %s))", eq_text, r2, p_formatted)

    eq_rows[[grp_lvl]] <- data.frame(
      Group = grp_lvl,
      Fit_Type = fit_type,
      Equation = eq_text,
      R_Squared = round(r2, 4),
      P_Value = round(p_val, 5),
      Label = full_eq_label,
      stringsAsFactors = FALSE
    )

    model_list[[grp_lvl]] <- mod

    # Build evaluation grid
    grid_x <- seq(min(df_fit$x), max(df_fit$x), length.out = 120)

    if (fit_type == "linear") {
      newdata <- data.frame(x = grid_x)
    } else if (fit_type == "power") {
      newdata <- data.frame(log_x = log(grid_x))
    } else if (fit_type == "exponential") {
      newdata <- data.frame(x = grid_x)
    } else if (fit_type == "logarithmic") {
      newdata <- data.frame(log_x = log(grid_x))
    } else if (fit_type == "polynomial") {
      newdata <- data.frame(x = grid_x, x2 = grid_x^2)
    } else if (fit_type == "cubic") {
      newdata <- data.frame(x = grid_x, x2 = grid_x^2, x3 = grid_x^3)
    }

    if (interval != "none") {
      pred_out <- as.data.frame(stats::predict(mod, newdata = newdata, interval = interval, level = 0.95))
      if (fit_type %in% c("power", "exponential")) {
        pred_out <- exp(pred_out)
      }
      p_df <- data.frame(
        x = grid_x,
        y_fit = pred_out$fit,
        y_lwr = pred_out$lwr,
        y_upr = pred_out$upr,
        Group = grp_lvl,
        stringsAsFactors = FALSE
      )
    } else {
      fit_vals <- stats::predict(mod, newdata = newdata)
      if (fit_type %in% c("power", "exponential")) fit_vals <- exp(fit_vals)
      p_df <- data.frame(
        x = grid_x,
        y_fit = fit_vals,
        y_lwr = NA_real_,
        y_upr = NA_real_,
        Group = grp_lvl,
        stringsAsFactors = FALSE
      )
    }

    pred_list[[grp_lvl]] <- p_df
  }

  eq_table <- do.call(rbind, eq_rows)
  pred_master <- do.call(rbind, pred_list)
  rownames(eq_table) <- NULL
  rownames(pred_master) <- NULL

  # 3. Base ggplot Assembly
  p <- ggplot2::ggplot()

  if (is.null(group_var)) {
    if (interval != "none") {
      p <- p + ggplot2::geom_ribbon(
        data = pred_master,
        ggplot2::aes(x = .data$x, ymin = .data$y_lwr, ymax = .data$y_upr),
        fill = "#4A90E2",
        alpha = 0.22
      )
    }
    p <- p + ggplot2::geom_line(
      data = pred_master,
      ggplot2::aes(x = .data$x, y = .data$y_fit),
      color = "#1A5276",
      linewidth = 1.1
    )
    p <- p + ggplot2::geom_point(
      data = clean_data,
      ggplot2::aes(x = .data[[x_var]], y = .data[[y_var]]),
      color = "#2E4053",
      size = point_size,
      alpha = point_alpha
    )

  } else {
    if (interval != "none") {
      p <- p + ggplot2::geom_ribbon(
        data = pred_master,
        ggplot2::aes(x = .data$x, ymin = .data$y_lwr, ymax = .data$y_upr, fill = .data$Group),
        alpha = 0.18,
        color = NA
      )
    }
    p <- p + ggplot2::geom_line(
      data = pred_master,
      ggplot2::aes(x = .data$x, y = .data$y_fit, color = .data$Group),
      linewidth = 1.1
    )
    p <- p + ggplot2::geom_point(
      data = clean_data,
      ggplot2::aes(x = .data[[x_var]], y = .data[[y_var]], color = .data[[group_var]], fill = .data[[group_var]]),
      size = point_size,
      alpha = point_alpha
    ) +
      ggplot2::scale_color_brewer(palette = color_palette, name = group_var) +
      ggplot2::scale_fill_brewer(palette = color_palette, name = group_var)
  }

  # 4. Inset Regression Equation
  if (show_equation) {
    x_min <- min(clean_data[[x_var]])
    x_max <- max(clean_data[[x_var]])
    y_min <- min(clean_data[[y_var]])
    y_max <- max(clean_data[[y_var]])

    ann_x <- x_min + 0.05 * (x_max - x_min)

    if (is.null(group_var) || facet) {
      ann_df <- eq_table
      ann_df$x <- ann_x
      ann_df$y <- y_max - 0.08 * (y_max - y_min)
      if (!is.null(group_var)) colnames(ann_df)[colnames(ann_df) == "Group"] <- group_var

      p <- p + ggplot2::geom_text(
        data = ann_df,
        ggplot2::aes(x = .data$x, y = .data$y, label = .data$Label),
        parse = TRUE,
        hjust = 0,
        vjust = 1,
        size = 3.6,
        color = "black",
        fontface = "bold",
        inherit.aes = FALSE
      )
    } else {
      n_grps <- nrow(eq_table)
      ann_df <- eq_table
      ann_df$x <- ann_x
      y_offsets <- seq(0.08, by = 0.12, length.out = n_grps)
      ann_df$y <- y_max - y_offsets * (y_max - y_min)
      colnames(ann_df)[colnames(ann_df) == "Group"] <- group_var

      p <- p + ggplot2::geom_text(
        data = ann_df,
        ggplot2::aes(x = .data$x, y = .data$y, label = .data$Label, color = .data[[group_var]]),
        parse = TRUE,
        hjust = 0,
        vjust = 1,
        size = 3.4,
        fontface = "bold",
        show.legend = FALSE,
        inherit.aes = FALSE
      )
    }
  }

  # 5. Faceting & Publication Styling
  if (!is.null(group_var) && facet) {
    p <- p + ggplot2::facet_wrap(stats::as.formula(paste("~", group_var)), scales = "free")
  }

  int_subtitle <- if (interval == "confidence") {
    "Fitted curve with 95% Confidence Interval (Mean Response)"
  } else if (interval == "prediction") {
    "Fitted curve with 95% Prediction Interval (Individual Observations)"
  } else {
    "Fitted curve without uncertainty intervals"
  }

  p <- p +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(
      axis.line = ggplot2::element_line(linewidth = 0.7, color = "black"),
      axis.text = ggplot2::element_text(color = "black", face = "bold"),
      axis.title = ggplot2::element_text(face = "bold"),
      legend.position = if (is.null(group_var)) "none" else "top",
      legend.title = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(face = "bold", size = 12, hjust = 0.5),
      plot.subtitle = ggplot2::element_text(size = 9.5, hjust = 0.5, color = "grey35"),
      strip.background = ggplot2::element_rect(fill = "grey92", color = "black"),
      strip.text = ggplot2::element_text(face = "bold")
    ) +
    ggplot2::labs(
      x = x_var,
      y = y_var,
      title = if (!is.null(title)) title else paste(tools::toTitleCase(fit_type), "Regression Analysis"),
      subtitle = int_subtitle
    )

  return(list(
    Plot = p,
    Models = model_list,
    Equation_Table = eq_table,
    Predictions = pred_master
  ))
}
