#' Two-Sample Comparison with Automated F-Test, t-Test, and Visualization
#'
#' Evaluates equality of variances using Fisher's \eqn{F}-test, checks group-level
#' normality via the Shapiro-Wilk test, and performs Student's, Welch's, or paired
#' \eqn{t}-test to compare means between two groups. Automatically computes effect sizes
#' (Cohen's \eqn{d}), compiles group descriptive statistics, and generates a
#' publication-ready \code{ggplot2} graphic with significance annotations.
#'
#' @param data A \code{data.frame}, \code{matrix}, or character string file path
#'   to a \code{.csv} or \code{.xlsx} file containing the observations.
#' @param response_var Character string; column name in \code{data} representing
#'   the continuous numeric response variable.
#' @param group_var Character string; column name in \code{data} representing
#'   the categorical grouping factor (must contain exactly two distinct levels).
#' @param var_equal Character string or logical; variance assumption for independent \eqn{t}-tests:
#'   \code{"auto"} (default; runs Fisher's \eqn{F}-test and adopts Welch's \eqn{t}-test if
#'   \eqn{p < 0.05}), \code{TRUE} (forces Student's \eqn{t}-test), or \code{FALSE}
#'   (forces Welch's \eqn{t}-test). Ignored if \code{paired = TRUE}. Default is \code{"auto"}.
#' @param alternative Character string; alternative hypothesis: \code{"two.sided"} (default),
#'   \code{"greater"}, or \code{"less"}. Default is \code{"two.sided"}.
#' @param paired Logical; if \code{TRUE}, performs a paired \eqn{t}-test. When \code{TRUE},
#'   observations must be matched row-wise across groups. Default is \code{FALSE}.
#' @param conf_level Numeric; confidence level for difference intervals and summary bounds
#'   (e.g., \code{0.95} for 95\% confidence intervals). Default is \code{0.95}.
#' @param plot_type Character string; visual geometry: \code{"box"} (default), \code{"violin"},
#'   or \code{"jitter"}. Default is \code{"box"}.
#' @param palette Character string; a valid \code{RColorBrewer} palette name. Default is \code{"Dark2"}.
#' @param title Optional character string; custom plot title. Default is \code{NULL}.
#'
#' @return A list containing:
#'   \item{Summary_Table}{Tidy data frame containing sample sizes (\eqn{N}), Mean, SD, SE, 95\% CI per group, mean difference, and Cohen's \eqn{d}.}
#'   \item{T_Test}{Tidy summary data frame of the \eqn{t}-test results (method, test statistic \eqn{t}, degrees of freedom, \eqn{p}-value, and confidence interval).}
#'   \item{F_Test}{Tidy summary data frame of Fisher's \eqn{F}-test for variance equality (variance ratio, degrees of freedom, and \eqn{p}-value).}
#'   \item{Normality_Test}{Tidy data frame of Shapiro-Wilk normality test results for each group.}
#'   \item{Plot}{A publication-ready \code{ggplot2} object annotated with test results.}
#'   \item{Model_Objects}{A list containing raw base R test outputs (\code{t_test}, \code{f_test}, and \code{shapiro_tests}).}
#' @export
#'
#' @import ggplot2
#' @importFrom stats t.test var.test shapiro.test sd qt na.omit
#' @importFrom tools file_ext
#' @importFrom utils read.csv
#' @importFrom rlang .data
compare_two_groups <- function(
    data,
    response_var,
    group_var,
    var_equal = "auto",
    alternative = "two.sided",
    paired = FALSE,
    conf_level = 0.95,
    plot_type = "box",
    palette = "Dark2",
    title = NULL
) {

  # 1. Defensive Ingestion Guard
  if (is.character(data) && length(data) == 1) {
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
    stop("Input 'data' must be a data.frame, tibble, matrix, or file path.")
  }

  data <- as.data.frame(data)

  # 2. Variable Validation
  if (!is.character(response_var) || length(response_var) != 1) {
    stop("'response_var' must be a single character string.")
  }
  if (!is.character(group_var) || length(group_var) != 1) {
    stop("'group_var' must be a single character string.")
  }

  missing_vars <- setdiff(c(response_var, group_var), names(data))
  if (length(missing_vars) > 0) {
    stop(sprintf("Variable(s) not found in data: %s", paste(missing_vars, collapse = ", ")))
  }

  # 3. Data Cleaning and Isolation
  sub_df <- data[, c(response_var, group_var)]
  names(sub_df) <- c("response", "group")

  n_initial <- nrow(sub_df)
  sub_df <- stats::na.omit(sub_df)
  n_dropped <- n_initial - nrow(sub_df)
  if (n_dropped > 0) {
    warning(sprintf("Removed %d observation(s) with missing values.", n_dropped), call. = FALSE)
  }

  if (nrow(sub_df) == 0) {
    stop("No complete observations remaining after removing missing values.")
  }

  if (!is.numeric(sub_df$response)) {
    stop(sprintf("Response variable '%s' must be numeric.", response_var))
  }

  sub_df$group <- as.factor(sub_df$group)
  grp_levels <- levels(sub_df$group)
  if (length(grp_levels) != 2) {
    stop(sprintf(
      "Grouping variable '%s' must have exactly 2 levels for a two-sample test. Found %d levels: %s.",
      group_var, length(grp_levels), paste(grp_levels, collapse = ", ")
    ))
  }

  # Split data by group vectors
  y1 <- sub_df$response[sub_df$group == grp_levels[1]]
  y2 <- sub_df$response[sub_df$group == grp_levels[2]]
  n1 <- length(y1)
  n2 <- length(y2)

  if (n1 < 2 || n2 < 2) {
    stop("Each group must contain at least 2 observations to calculate variances and t-tests.")
  }

  if (paired && n1 != n2) {
    stop(sprintf("Paired t-test requires equal sample sizes per group. Found n1 = %d and n2 = %d.", n1, n2))
  }

  # 4. Diagnostic Tests: Normality (Shapiro-Wilk) & Variance Equality (F-test)
  shapiro_list <- list()
  shapiro_rows <- list()
  for (grp in grp_levels) {
    y_g <- sub_df$response[sub_df$group == grp]
    if (length(y_g) >= 3 && length(y_g) <= 5000) {
      if (length(unique(y_g)) == 1) {
        shapiro_rows[[grp]] <- data.frame(
          Group = grp, Statistic_W = NA_real_, p_value = NA_real_,
          Normal_05 = "Identical Values", stringsAsFactors = FALSE
        )
      } else {
        sw <- stats::shapiro.test(y_g)
        shapiro_list[[grp]] <- sw
        shapiro_rows[[grp]] <- data.frame(
          Group = grp,
          Statistic_W = round(sw$statistic, 4),
          p_value = round(sw$p.value, 4),
          Normal_05 = if (sw$p.value >= 0.05) "Yes" else "No",
          stringsAsFactors = FALSE
        )
      }
    } else {
      shapiro_rows[[grp]] <- data.frame(
        Group = grp, Statistic_W = NA_real_, p_value = NA_real_,
        Normal_05 = "N out of bounds (3-5000)", stringsAsFactors = FALSE
      )
    }
  }
  shapiro_df <- do.call(rbind, shapiro_rows)
  rownames(shapiro_df) <- NULL

  # Fisher's F-test on the two numeric vectors (bypasses formula limitations)
  f_test_res <- stats::var.test(x = y1, y = y2, conf.level = conf_level)
  f_df <- data.frame(
    Ratio_Var = round(unname(f_test_res$estimate), 4),
    F_Statistic = round(unname(f_test_res$statistic), 4),
    df_num = unname(f_test_res$parameter[1]),
    df_denom = unname(f_test_res$parameter[2]),
    p_value = round(f_test_res$p.value, 4),
    Equal_Var_05 = if (f_test_res$p.value >= 0.05) "Yes" else "No",
    stringsAsFactors = FALSE
  )
  rownames(f_df) <- NULL

  # 5. Determine Variance Assumption & Execute t-Test via Vectors
  if (paired) {
    use_var_equal <- FALSE
  } else if (isTRUE(var_equal)) {
    use_var_equal <- TRUE
  } else if (isFALSE(var_equal)) {
    use_var_equal <- FALSE
  } else if (is.character(var_equal) && tolower(var_equal) == "auto") {
    use_var_equal <- f_test_res$p.value >= 0.05
  } else {
    stop("'var_equal' must be 'auto', TRUE, or FALSE.")
  }

  # Vector-based call avoids t.test.formula rejection of paired argument
  if (paired) {
    t_test_res <- stats::t.test(
      x = y1,
      y = y2,
      paired = TRUE,
      alternative = alternative,
      conf.level = conf_level
    )
  } else {
    t_test_res <- stats::t.test(
      x = y1,
      y = y2,
      paired = FALSE,
      var.equal = use_var_equal,
      alternative = alternative,
      conf.level = conf_level
    )
  }

  t_df <- data.frame(
    Test_Type = t_test_res$method,
    Alternative = alternative,
    t_Statistic = round(unname(t_test_res$statistic), 4),
    df = round(unname(t_test_res$parameter), 2),
    p_value = round(t_test_res$p.value, 4),
    CI_Lower = round(t_test_res$conf.int[1], 4),
    CI_Upper = round(t_test_res$conf.int[2], 4),
    stringsAsFactors = FALSE
  )
  rownames(t_df) <- NULL

  # 6. Group Summaries and Effect Size (Cohen's d)
  m1 <- mean(y1); s1 <- stats::sd(y1); se1 <- s1 / sqrt(n1)
  m2 <- mean(y2); s2 <- stats::sd(y2); se2 <- s2 / sqrt(n2)
  crit_t1 <- stats::qt((1 + conf_level) / 2, df = n1 - 1)
  crit_t2 <- stats::qt((1 + conf_level) / 2, df = n2 - 1)

  if (paired) {
    diffs <- y1 - y2
    sd_diff <- stats::sd(diffs)
    cohen_d <- if (sd_diff > 0) abs(mean(diffs)) / sd_diff else 0
  } else {
    if (use_var_equal) {
      pooled_sd <- sqrt(((n1 - 1) * s1^2 + (n2 - 1) * s2^2) / (n1 + n2 - 2))
    } else {
      pooled_sd <- sqrt((s1^2 + s2^2) / 2)
    }
    cohen_d <- if (pooled_sd > 0) abs(m1 - m2) / pooled_sd else 0
  }

  summary_table <- data.frame(
    Group = grp_levels,
    N = c(n1, n2),
    Mean = round(c(m1, m2), 3),
    SD = round(c(s1, s2), 3),
    SE = round(c(se1, se2), 3),
    CI_Lower = round(c(m1 - crit_t1 * se1, m2 - crit_t2 * se2), 3),
    CI_Upper = round(c(m1 + crit_t1 * se1, m2 + crit_t2 * se2), 3),
    Mean_Diff = round(c(m1 - m2, m2 - m1), 3),
    Cohen_d = round(cohen_d, 3),
    stringsAsFactors = FALSE
  )

  # 7. Publication-Ready ggplot2 Graphic
  p_val_display <- if (t_test_res$p.value < 0.001) {
    "p < 0.001"
  } else {
    sprintf("p = %.4f", t_test_res$p.value)
  }

  if (paired) {
    sub_label <- sprintf("Paired t-test: t(%.1f) = %.2f, %s", t_test_res$parameter, t_test_res$statistic, p_val_display)
  } else {
    sub_label <- sprintf(
      "%s: t(%.1f) = %.2f, %s | F-test: %s (var.equal = %s)",
      if (use_var_equal) "Student's t" else "Welch's t",
      t_test_res$parameter,
      t_test_res$statistic,
      p_val_display,
      if (f_test_res$p.value >= 0.05) "p >= 0.05" else "p < 0.05",
      as.character(use_var_equal)
    )
  }

  max_val <- max(sub_df$response)
  min_val <- min(sub_df$response)
  y_span <- if (max_val == min_val) 1 else (max_val - min_val)
  bracket_y <- max_val + 0.08 * y_span
  tip_y <- bracket_y - 0.02 * y_span
  text_y <- bracket_y + 0.03 * y_span

  p <- ggplot2::ggplot(sub_df, ggplot2::aes(x = .data$group, y = .data$response, fill = .data$group))

  if (plot_type == "box") {
    p <- p +
      ggplot2::geom_boxplot(width = 0.5, alpha = 0.8, outlier.shape = NA, color = "black", linewidth = 0.6) +
      ggplot2::geom_jitter(width = 0.15, alpha = 0.5, size = 2, shape = 21, color = "black")
  } else if (plot_type == "violin") {
    p <- p +
      ggplot2::geom_violin(trim = FALSE, alpha = 0.7, color = "black", linewidth = 0.6) +
      ggplot2::geom_boxplot(width = 0.15, fill = "white", color = "black", outlier.shape = NA, linewidth = 0.5)
  } else {
    p <- p +
      ggplot2::geom_jitter(width = 0.2, alpha = 0.7, size = 2.5, shape = 21, color = "black") +
      ggplot2::stat_summary(fun.data = "mean_cl_normal", geom = "errorbar", width = 0.2, linewidth = 0.8, color = "black") +
      ggplot2::stat_summary(fun = "mean", geom = "point", size = 3.5, color = "red")
  }

  p <- p +
    ggplot2::scale_fill_brewer(palette = palette) +
    ggplot2::annotate("segment", x = 1, xend = 2, y = bracket_y, yend = bracket_y, linewidth = 0.6) +
    ggplot2::annotate("segment", x = 1, xend = 1, y = bracket_y, yend = tip_y, linewidth = 0.6) +
    ggplot2::annotate("segment", x = 2, xend = 2, y = bracket_y, yend = tip_y, linewidth = 0.6) +
    ggplot2::annotate("text", x = 1.5, y = text_y, label = p_val_display, fontface = "bold", size = 4) +
    ggplot2::coord_cartesian(ylim = c(min_val - 0.05 * y_span, bracket_y + 0.08 * y_span)) +
    ggplot2::labs(
      title = if (!is.null(title)) title else sprintf("Group Comparison: %s by %s", response_var, group_var),
      subtitle = sub_label,
      x = group_var,
      y = response_var
    ) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(
      legend.position = "none",
      plot.title = ggplot2::element_text(face = "bold", size = 13, hjust = 0.5),
      plot.subtitle = ggplot2::element_text(color = "grey25", size = 10, hjust = 0.5),
      axis.title = ggplot2::element_text(face = "bold"),
      axis.text = ggplot2::element_text(color = "black")
    )

  # 8. Return Structured Output List
  return(list(
    Summary_Table = summary_table,
    T_Test = t_df,
    F_Test = f_df,
    Normality_Test = shapiro_df,
    Plot = p,
    Model_Objects = list(
      t_test = t_test_res,
      f_test = f_test_res,
      shapiro_tests = shapiro_list
    )
  ))
}
