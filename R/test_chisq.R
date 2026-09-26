#' Chi-Square Test of Independence and Goodness-of-Fit with Visualization
#'
#' Performs Pearson's Chi-Square (\eqn{\chi^2}) test of independence or goodness-of-fit
#' on categorical survey, ecological, or environmental frequency data. Automatically
#' evaluates expected cell frequencies, applies Yates' continuity correction or Monte Carlo
#' simulation when sample sizes are small, computes effect sizes (Cramér's \eqn{V} and
#' Phi coefficient \eqn{\phi}), extracts standardized Pearson residuals, and generates a
#' publication-ready \code{ggplot2} stacked or dodged bar chart.
#'
#' @param data A \code{data.frame}, \code{matrix}, \code{table}, or character string file path
#'   to a \code{.csv} or \code{.xlsx} file containing the observations.
#' @param x Character string; column name in \code{data} representing the primary categorical
#'   factor (rows of the contingency table).
#' @param y Optional character string; column name in \code{data} representing the secondary
#'   categorical factor (columns of the contingency table). If \code{NULL}, a Chi-Square
#'   goodness-of-fit test is performed on \code{x}. Default is \code{NULL}.
#' @param p Optional numeric vector; theoretical probabilities for the goodness-of-fit test
#'   when \code{y = NULL}. Must sum to 1 and match the number of levels in \code{x}.
#'   Default is \code{NULL} (equal probabilities across categories).
#' @param correct Logical; if \code{TRUE}, applies Yates' continuity correction for
#'   \eqn{2 \times 2} contingency tables. Default is \code{TRUE}.
#' @param simulate_p Logical; if \code{TRUE}, computes \eqn{p}-values via Monte Carlo simulation
#'   when expected cell counts fall below 5. Default is \code{TRUE}.
#' @param B Numeric integer; number of Monte Carlo replicates used when \code{simulate_p = TRUE}.
#'   Default is \code{2000}.
#' @param plot_type Character string; visual geometry: \code{"stacked"} (100\% stacked proportion bar chart; default)
#'   or \code{"dodged"} (grouped absolute count bar chart). Default is \code{"stacked"}.
#' @param palette Character string; a valid \code{RColorBrewer} palette name. Default is \code{"Set2"}.
#' @param title Optional character string; custom plot title. Default is \code{NULL}.
#'
#' @return A list containing:
#'   \item{Chi_Square_Test}{Tidy data frame containing the \eqn{\chi^2} statistic, degrees of freedom, \eqn{p}-value, test method, and effect sizes (Cramér's \eqn{V} and Phi \eqn{\phi}).}
#'   \item{Observed_Counts}{Cross-tabulated data frame of observed cell frequencies with marginal totals.}
#'   \item{Expected_Counts}{Data frame of theoretical expected cell frequencies under the null hypothesis of independence.}
#'   \item{Standardized_Residuals}{Data frame of standardized Pearson residuals (\eqn{|r| > 2} indicates significant deviation).}
#'   \item{Plot}{A publication-ready \code{ggplot2} bar chart annotated with test statistics.}
#'   \item{Model_Object}{The underlying \code{htest} object returned by \code{stats::chisq.test}.}
#' @export
#'
#' @import ggplot2
#' @importFrom stats chisq.test na.omit
#' @importFrom tools file_ext
#' @importFrom utils read.csv
#' @importFrom rlang .data
test_chisq <- function(
    data,
    x,
    y = NULL,
    p = NULL,
    correct = TRUE,
    simulate_p = TRUE,
    B = 2000,
    plot_type = "stacked",
    palette = "Set2",
    title = NULL
) {

  # 1. Defensive Ingestion Guard
  if (inherits(data, "table")) {
    data <- as.data.frame(data)
  } else if (is.character(data) && length(data) == 1) {
    if (!file.exists(data)) {
      stop(sprintf("File does not exist: '%s'", data))
    }
    ext <- tolower(tools::file_ext(data))
    if (ext == "csv") {
      data <- utils::read.csv(data, stringsAsFactors = FALSE)
    } else if (ext %in% c("xlsx", "xls")) {
      if (!requireNamespace("readxl", quietly = TRUE)) {
        stop("Package 'readxl' is required to read Excel files. Install it or provide a data.frame.")
      }
      data <- as.data.frame(readxl::read_excel(data))
    } else {
      stop("Unsupported file type. Please provide a .csv, .xlsx file, or a data.frame.")
    }
  } else if (is.matrix(data) || inherits(data, "tbl_df") || inherits(data, "tbl")) {
    data <- as.data.frame(data)
  } else if (!is.data.frame(data)) {
    stop("Input 'data' must be a data.frame, tibble, matrix, table, or file path.")
  }

  data <- as.data.frame(data)

  # 2. Variable Validation
  if (!is.character(x) || length(x) != 1) {
    stop("'x' must be a single character string representing a column name.")
  }
  target_vars <- x
  if (!is.null(y)) {
    if (!is.character(y) || length(y) != 1) {
      stop("'y' must be a single character string representing a column name.")
    }
    target_vars <- c(target_vars, y)
  }

  missing_vars <- setdiff(target_vars, names(data))
  if (length(missing_vars) > 0) {
    stop(sprintf("Variable(s) not found in data: %s", paste(missing_vars, collapse = ", ")))
  }

  # 3. Clean and Isolate Observations
  sub_df <- data[, target_vars, drop = FALSE]
  n_initial <- nrow(sub_df)
  sub_df <- stats::na.omit(sub_df)
  n_dropped <- n_initial - nrow(sub_df)
  if (n_dropped > 0) {
    warning(sprintf("Removed %d observation(s) with missing values.", n_dropped), call. = FALSE)
  }

  if (nrow(sub_df) == 0) {
    stop("No complete observations remaining after removing missing values.")
  }

  sub_df[[x]] <- as.factor(sub_df[[x]])

  # 4. Statistical Execution (Independence vs. Goodness-of-Fit)
  is_gof <- is.null(y)

  if (is_gof) {
    # Goodness-of-Fit Test
    obs_tab <- table(sub_df[[x]])
    n_total <- sum(obs_tab)

    if (is.null(p)) {
      p_expected <- rep(1 / length(obs_tab), length(obs_tab))
    } else {
      if (length(p) != length(obs_tab)) {
        stop(sprintf("Length of 'p' (%d) must match the number of levels in 'x' (%d).", length(p), length(obs_tab)))
      }
      if (abs(sum(p) - 1) > 1e-5) {
        stop("Probabilities in 'p' must sum to 1.")
      }
      p_expected <- p
    }

    # Preliminary expected check
    exp_counts <- n_total * p_expected
    small_cells <- any(exp_counts < 5)

    if (small_cells && isTRUE(simulate_p)) {
      warning("Expected cell counts < 5 detected. Applying Monte Carlo simulation for p-value estimation.", call. = FALSE)
      chi_obj <- stats::chisq.test(obs_tab, p = p_expected, rescale.p = FALSE, simulate.p.value = TRUE, B = B)
    } else {
      chi_obj <- stats::chisq.test(obs_tab, p = p_expected, rescale.p = FALSE, correct = correct)
    }

    cramers_v <- NA_real_
    phi <- NA_real_

    # Format output tables
    obs_df <- as.data.frame(obs_tab)
    names(obs_df) <- c(x, "Observed")
    obs_df$Percentage <- round(100 * obs_df$Observed / n_total, 2)

    exp_df <- data.frame(Category = names(obs_tab), Expected = round(as.numeric(chi_obj$expected), 2), stringsAsFactors = FALSE)
    names(exp_df)[1] <- x

    res_df <- data.frame(Category = names(obs_tab), Std_Residual = round(as.numeric(chi_obj$residuals), 3), stringsAsFactors = FALSE)
    names(res_df)[1] <- x

  } else {
    # Test of Independence
    sub_df[[y]] <- as.factor(sub_df[[y]])
    obs_tab <- table(sub_df[[x]], sub_df[[y]])
    n_total <- sum(obs_tab)

    r <- nrow(obs_tab)
    c <- ncol(obs_tab)

    if (r < 2 || c < 2) {
      stop("Both 'x' and 'y' must contain at least 2 distinct factor levels.")
    }

    # Assess expected counts
    initial_exp <- outer(rowSums(obs_tab), colSums(obs_tab)) / n_total
    small_cells <- any(initial_exp < 5)

    if (small_cells && isTRUE(simulate_p)) {
      warning("Expected cell counts < 5 detected. Applying Monte Carlo simulation for p-value estimation.", call. = FALSE)
      chi_obj <- stats::chisq.test(obs_tab, simulate.p.value = TRUE, B = B)
    } else {
      chi_obj <- stats::chisq.test(obs_tab, correct = correct)
    }

    # Effect Sizes
    chi_val <- unname(chi_obj$statistic)
    min_dim <- min(r - 1, c - 1)
    cramers_v <- if (min_dim > 0 && n_total > 0) sqrt(chi_val / (n_total * min_dim)) else 0
    phi <- if (r == 2 && c == 2 && n_total > 0) sqrt(chi_val / n_total) else NA_real_

    # Format tables
    obs_df <- as.data.frame.matrix(obs_tab)
    obs_df$Total <- rowSums(obs_df)
    total_row <- as.data.frame(t(colSums(obs_df)))
    rownames(total_row) <- "Total"
    obs_df <- rbind(obs_df, total_row)
    obs_df <- cbind(Category = rownames(obs_df), obs_df)
    names(obs_df)[1] <- x
    rownames(obs_df) <- NULL

    exp_df <- as.data.frame(round(chi_obj$expected, 2))
    exp_df <- cbind(Category = rownames(exp_df), exp_df)
    names(exp_df)[1] <- x
    rownames(exp_df) <- NULL

    res_df <- as.data.frame(round(chi_obj$stdres, 3))
    res_df <- cbind(Category = rownames(res_df), res_df)
    names(res_df)[1] <- x
    rownames(res_df) <- NULL
  }

  # 5. Summary Test Table
  p_val <- chi_obj$p.value
  p_val_display <- if (p_val < 0.001) "< 0.001" else sprintf("%.4f", p_val)

  test_summary <- data.frame(
    Test_Type   = chi_obj$method,
    Chi_Square  = round(unname(chi_obj$statistic), 4),
    df          = if (!is.null(chi_obj$parameter)) unname(chi_obj$parameter) else NA_integer_,
    p_value     = round(p_val, 4),
    Cramers_V   = round(cramers_v, 4),
    Phi         = round(phi, 4),
    stringsAsFactors = FALSE
  )

  # 6. Publication-Ready ggplot2 Graphic
  sub_label <- sprintf(
    "Chi-Square: X2 = %.2f%s, p = %s%s",
    unname(chi_obj$statistic),
    if (!is.na(test_summary$df)) sprintf(", df = %d", test_summary$df) else "",
    p_val_display,
    if (!is.na(cramers_v)) sprintf(" | Cram\\u00e9r's V = %.2f", cramers_v) else ""
  )

  if (is_gof) {
    plot_df <- as.data.frame(obs_tab)
    names(plot_df) <- c("Category", "Count")
    plot_df$Prop <- plot_df$Count / sum(plot_df$Count)

    p_plot <- ggplot2::ggplot(plot_df, ggplot2::aes(x = .data$Category, y = .data$Count, fill = .data$Category)) +
      ggplot2::geom_col(width = 0.6, color = "black", linewidth = 0.5, alpha = 0.85) +
      ggplot2::geom_text(
        ggplot2::aes(label = sprintf("%d (%.1f%%)", .data$Count, .data$Prop * 100)),
        vjust = -0.4, fontface = "bold", size = 3.8
      ) +
      ggplot2::scale_fill_brewer(palette = palette) +
      ggplot2::labs(
        title = if (!is.null(title)) title else sprintf("Goodness-of-Fit Distribution: %s", x),
        subtitle = sub_label,
        x = x,
        y = "Observed Frequency"
      ) +
      ggplot2::theme_classic(base_size = 12) +
      ggplot2::theme(
        legend.position = "none",
        plot.title = ggplot2::element_text(face = "bold", size = 13, hjust = 0.5),
        plot.subtitle = ggplot2::element_text(color = "grey25", size = 10, hjust = 0.5),
        axis.title = ggplot2::element_text(face = "bold"),
        axis.text = ggplot2::element_text(color = "black")
      )
  } else {
    plot_df <- as.data.frame(table(sub_df[[x]], sub_df[[y]]))
    names(plot_df) <- c("X_Var", "Y_Var", "Count")

    if (plot_type == "stacked") {
      p_plot <- ggplot2::ggplot(plot_df, ggplot2::aes(x = .data$X_Var, y = .data$Count, fill = .data$Y_Var)) +
        ggplot2::geom_col(position = "fill", color = "black", linewidth = 0.5, alpha = 0.85, width = 0.65) +
        ggplot2::scale_y_continuous(labels = function(x) sprintf("%.0f%%", x * 100)) +
        ggplot2::scale_fill_brewer(palette = palette, name = y) +
        ggplot2::labs(
          title = if (!is.null(title)) title else sprintf("Contingency Proportion: %s across %s", y, x),
          subtitle = sub_label,
          x = x,
          y = "Relative Proportion"
        )
    } else {
      p_plot <- ggplot2::ggplot(plot_df, ggplot2::aes(x = .data$X_Var, y = .data$Count, fill = .data$Y_Var)) +
        ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.8), color = "black", linewidth = 0.5, alpha = 0.85, width = 0.7) +
        ggplot2::scale_fill_brewer(palette = palette, name = y) +
        ggplot2::labs(
          title = if (!is.null(title)) title else sprintf("Contingency Counts: %s across %s", y, x),
          subtitle = sub_label,
          x = x,
          y = "Observed Count"
        )
    }

    p_plot <- p_plot +
      ggplot2::theme_classic(base_size = 12) +
      ggplot2::theme(
        legend.position = "right",
        legend.title = ggplot2::element_text(face = "bold"),
        plot.title = ggplot2::element_text(face = "bold", size = 13, hjust = 0.5),
        plot.subtitle = ggplot2::element_text(color = "grey25", size = 10, hjust = 0.5),
        axis.title = ggplot2::element_text(face = "bold"),
        axis.text = ggplot2::element_text(color = "black")
      )
  }

  # 7. Return Structured Output List
  return(list(
    Chi_Square_Test      = test_summary,
    Observed_Counts      = obs_df,
    Expected_Counts      = exp_df,
    Standardized_Residuals = res_df,
    Plot                 = p_plot,
    Model_Object         = chi_obj
  ))
}
