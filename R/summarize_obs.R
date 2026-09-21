#' Comprehensive Observation Summary with Excel Export
#'
#' Computes descriptive statistics (N, Mean, SD, SE, Median, Min, Max, and Mean \eqn{\pm} SD)
#' for numeric variables across optional grouping factors. Supports direct file
#' ingestion and exporting publication-ready tables directly to Excel (.xlsx) or CSV.
#'
#' @param data A \code{data.frame}, \code{matrix}, or file path to a \code{.csv} or \code{.xlsx} file.
#' @param vars Character vector. Names of numeric continuous variables to summarize.
#'   If \code{NULL}, all numeric variables in the dataset are evaluated.
#' @param group Optional character vector. Grouping variable(s) (e.g., \code{"Habitat"} or \code{c("Habitat", "Season")}).
#' @param digits Integer. Number of decimal places for rounding (default: 2).
#' @param file Optional character. File path to export the resulting summary table
#'   (must end in \code{.xlsx} or \code{.csv}).
#'
#' @return A tidy \code{data.frame} containing computed descriptive metrics.
#' @export
#'
#' @importFrom stats sd median na.omit
#' @importFrom utils write.csv read.csv
#' @importFrom tools file_ext
#'
#' @examples
#' \dontrun{
#' summary_tbl <- summarize_obs(
#'   data = fish_survey,
#'   vars = c("Length", "Weight"),
#'   group = "River",
#'   digits = 2
#' )
#' head(summary_tbl)
#' }
summarize_obs <- function(data,
                          vars = NULL,
                          group = NULL,
                          digits = 2,
                          file = NULL) {

  # 1. Ingest Data Formats (Base R, zero tidyverse overhead)
  if (is.character(data) && length(data) == 1) {
    if (!file.exists(data)) stop(sprintf("File '%s' not found.", data))
    ext <- tolower(tools::file_ext(data))
    if (ext == "csv") {
      df <- utils::read.csv(data, stringsAsFactors = FALSE)
    } else if (ext %in% c("xlsx", "xls")) {
      if (!requireNamespace("readxl", quietly = TRUE)) {
        stop("Please install 'readxl' to load Excel files.")
      }
      df <- as.data.frame(readxl::read_excel(data))
    } else {
      stop("Unsupported file type. Use data.frame, matrix, .csv, or .xlsx.")
    }
  } else if (is.matrix(data)) {
    df <- as.data.frame(data)
  } else if (is.data.frame(data)) {
    df <- as.data.frame(data)
  } else {
    stop("'data' must be a data.frame, matrix, or valid file path.")
  }

  # 2. Variable Selection & Validation
  if (is.null(vars)) {
    vars <- names(df)[vapply(df, is.numeric, logical(1))]
    if (length(vars) == 0) stop("No numeric variables detected in the dataset.")
  } else {
    missing_vars <- vars[!vars %in% names(df)]
    if (length(missing_vars) > 0) {
      stop(sprintf("Variable(s) not found: %s", paste(missing_vars, collapse = ", ")))
    }
    non_num <- vars[!vapply(df[vars], is.numeric, logical(1))]
    if (length(non_num) > 0) {
      stop(sprintf("Variable(s) must be numeric: %s", paste(non_num, collapse = ", ")))
    }
  }

  if (!is.null(group)) {
    missing_grp <- group[!group %in% names(df)]
    if (length(missing_grp) > 0) {
      stop(sprintf("Grouping variable(s) not found: %s", paste(missing_grp, collapse = ", ")))
    }
  }

  # 3. Compute Descriptive Statistics
  compute_stats <- function(v) {
    v_clean <- stats::na.omit(v)
    n <- length(v_clean)
    if (n == 0) {
      return(c(N = 0, Mean = NA, SD = NA, SE = NA, Median = NA, Min = NA, Max = NA))
    }
    m <- mean(v_clean)
    s <- if (n > 1) stats::sd(v_clean) else 0
    se <- if (n > 1) s / sqrt(n) else 0
    med <- stats::median(v_clean)
    mn <- min(v_clean)
    mx <- max(v_clean)
    c(N = n, Mean = m, SD = s, SE = se, Median = med, Min = mn, Max = mx)
  }

  out_list <- list()

  for (v_name in vars) {
    if (is.null(group)) {
      st <- compute_stats(df[[v_name]])
      row_df <- data.frame(
        Variable = v_name,
        N = as.integer(st["N"]),
        Mean = round(st["Mean"], digits),
        SD = round(st["SD"], digits),
        SE = round(st["SE"], digits),
        Median = round(st["Median"], digits),
        Min = round(st["Min"], digits),
        Max = round(st["Max"], digits),
        Mean_SD = sprintf(paste0("%.", digits, "f \u00b1 %.", digits, "f"), st["Mean"], st["SD"]),
        stringsAsFactors = FALSE
      )
      out_list[[v_name]] <- row_df
    } else {
      grp_list <- df[, group, drop = FALSE]
      split_vals <- split(df[[v_name]], grp_list, drop = TRUE, sep = "___")

      sub_rows <- list()
      for (grp_key in names(split_vals)) {
        st <- compute_stats(split_vals[[grp_key]])
        grp_parts <- unlist(strsplit(grp_key, "___", fixed = TRUE))

        meta <- as.data.frame(as.list(grp_parts), stringsAsFactors = FALSE)
        colnames(meta) <- group

        meta$Variable <- v_name
        meta$N <- as.integer(st["N"])
        meta$Mean <- round(st["Mean"], digits)
        meta$SD <- round(st["SD"], digits)
        meta$SE <- round(st["SE"], digits)
        meta$Median <- round(st["Median"], digits)
        meta$Min <- round(st["Min"], digits)
        meta$Max <- round(st["Max"], digits)
        meta$Mean_SD <- sprintf(paste0("%.", digits, "f \u00b1 %.", digits, "f"), st["Mean"], st["SD"])

        sub_rows[[length(sub_rows) + 1]] <- meta
      }
      out_list[[v_name]] <- do.call(rbind, sub_rows)
    }
  }

  final_summary <- do.call(rbind, out_list)
  rownames(final_summary) <- NULL

  # 4. Optional Export to Excel (.xlsx) or CSV
  if (!is.null(file)) {
    out_ext <- tolower(tools::file_ext(file))
    if (out_ext == "xlsx") {
      if (!requireNamespace("writexl", quietly = TRUE)) {
        stop("Please install 'writexl' using install.packages('writexl') to export to Excel.")
      }
      writexl::write_xlsx(final_summary, path = file)
      message(sprintf("Summary successfully saved to Excel file: %s", file))
    } else if (out_ext == "csv") {
      utils::write.csv(final_summary, file = file, row.names = FALSE)
      message(sprintf("Summary successfully saved to CSV file: %s", file))
    } else {
      warning("Unsupported export format extension. Table was computed but not saved to disk.")
    }
  }

  return(final_summary)
}
