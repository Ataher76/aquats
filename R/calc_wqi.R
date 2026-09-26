#' Weighted Arithmetic Water Quality Index (WAWQI) with Confidence Intervals
#'
#' Computes the Weighted Arithmetic Water Quality Index (WAWQI) per observation
#' or replicate. When replicated station sampling is detected, it calculates
#' station-level summary statistics including Mean, SD, SE, 95\% confidence intervals,
#' and overall ecological water quality rating.
#'
#' @param data A \code{data.frame}, \code{matrix}, or character string file path
#'   to a \code{.csv} or \code{.xlsx} file containing physicochemical water quality parameters.
#' @param id_col Character string; column name in \code{data} identifying sampling stations
#'   or monitoring sites. Default is \code{"Station"}.
#' @param standards Optional named numeric vector of standard regulatory limits (\eqn{S_i})
#'   in matching concentration units. If \code{NULL}, default benchmark guidelines (e.g., WHO/BIS)
#'   are applied. Default is \code{NULL}.
#' @param ideal_values Optional named numeric vector of ideal zero-effect reference values (\eqn{V_0})
#'   (e.g., 7.0 for pH, 14.6 mg/L for dissolved oxygen). Default is \code{NULL}.
#' @param ci Numeric; confidence level for interval estimation (e.g., \code{0.95} for 95\% confidence intervals).
#'   Default is \code{0.95}.
#' @param file Optional character string; export file path for tabular results
#'   (\code{.xlsx} or \code{.csv}). Default is \code{NULL}.
#'
#' @return A list containing:
#'   \item{Replicate_Data}{Data frame containing computed sub-indices, relative weights, and overall WAWQI for each individual sample observation.}
#'   \item{Summary_Data}{Data frame of station-level aggregated statistics (Mean, SD, SE, 95\% CI, and water quality classification rating).}
#' @export
#'
#' @importFrom stats qt sd na.omit
#' @importFrom tools file_ext
#' @importFrom utils read.csv write.csv
calc_wqi <- function(
    data,
    id_col = "Station",
    standards = c(
      pH = 8.5,
      DO = 5.0,
      BOD = 5.0,
      TDS = 500,
      Turbidity = 5.0
    ),
    ideal_values = NULL,
    ci = 0.95,
    file = NULL
) {

  # 1. Standard Defensive Ingestion & File Support
  if (is.character(data) && length(data) == 1) {
    if (!file.exists(data)) stop(sprintf("File '%s' not found.", data))
    ext <- tolower(tools::file_ext(data))
    if (ext == "csv") {
      df <- utils::read.csv(data, stringsAsFactors = FALSE)
    } else if (ext %in% c("xlsx", "xls")) {
      if (!requireNamespace("readxl", quietly = TRUE)) {
        stop("Please install 'readxl' to read Excel files.")
      }
      df <- as.data.frame(readxl::read_excel(data))
    } else {
      stop("Unsupported file format. Provide a .csv or .xlsx file.")
    }
  } else if (is.matrix(data)) {
    df <- as.data.frame(data)
  } else if (is.data.frame(data)) {
    df <- as.data.frame(data)
  } else {
    stop("Input 'data' must be a data.frame, matrix, or valid file path.")
  }

  if (!id_col %in% names(df)) {
    df[[id_col]] <- paste0("Sample_", seq_len(nrow(df)))
  }

  # 2. Check Matching Parameters
  param_names <- names(standards)
  matched_params <- param_names[param_names %in% names(df)]

  if (length(matched_params) < 2) {
    stop(sprintf(
      "At least 2 matching parameters are required. Found: %s",
      paste(matched_params, collapse = ", ")
    ))
  }

  std_active <- standards[matched_params]

  # 3. Setup Ideal Values (V0)
  v0 <- rep(0, length(matched_params))
  names(v0) <- matched_params
  if ("pH" %in% matched_params) v0["pH"] <- 7.0
  if ("DO" %in% matched_params) v0["DO"] <- 14.6

  if (!is.null(ideal_values)) {
    for (nm in names(ideal_values)) {
      if (nm %in% matched_params) v0[nm] <- ideal_values[nm]
    }
  }

  # 4. Compute Unit Weights (Wi)
  k_val <- 1 / sum(1 / std_active)
  unit_weights <- (k_val / std_active)

  # 5. Compute Sub-Indices (Qi) & WQI per Replicate
  n_rows <- nrow(df)
  qi_matrix <- matrix(NA_real_, nrow = n_rows, ncol = length(matched_params))
  colnames(qi_matrix) <- paste0("Q_", matched_params)

  for (j in seq_along(matched_params)) {
    p_name <- matched_params[j]
    obs_val <- as.numeric(df[[p_name]])
    s_val <- std_active[p_name]
    ideal <- v0[p_name]

    if (p_name == "DO") {
      q_val <- ((ideal - obs_val) / (ideal - s_val)) * 100
    } else {
      q_val <- ((obs_val - ideal) / (s_val - ideal)) * 100
    }
    qi_matrix[, j] <- pmax(0, q_val)
  }

  wqi_values <- apply(qi_matrix, 1, function(row_q) {
    if (any(is.na(row_q))) return(NA_real_)
    sum(row_q * unit_weights) / sum(unit_weights)
  })

  raw_out <- data.frame(
    Station = df[[id_col]],
    df[, matched_params, drop = FALSE],
    round(qi_matrix, 2),
    WQI = round(wqi_values, 2),
    stringsAsFactors = FALSE
  )

  # 6. Station-Level Aggregation & 95% Confidence Intervals
  compute_ci <- function(vec) {
    v <- stats::na.omit(vec)
    n <- length(v)
    m <- mean(v)
    if (n < 2) {
      return(c(N = n, Mean = m, SD = 0, SE = 0, LCI = m, UCI = m))
    }
    s <- stats::sd(v)
    se <- s / sqrt(n)
    alpha <- 1 - ci
    t_crit <- stats::qt(1 - alpha / 2, df = n - 1)
    c(
      N = n,
      Mean = m,
      SD = s,
      SE = se,
      LCI = m - (t_crit * se),
      UCI = m + (t_crit * se)
    )
  }

  split_wqi <- split(raw_out$WQI, raw_out$Station)
  summary_rows <- lapply(names(split_wqi), function(stn) {
    stats_m <- compute_ci(split_wqi[[stn]])
    data.frame(
      Station  = stn,
      N        = as.integer(stats_m["N"]),
      WQI_Mean = round(stats_m["Mean"], 2),
      WQI_SD   = round(stats_m["SD"], 2),
      WQI_SE   = round(stats_m["SE"], 2),
      WQI_LCI  = round(stats_m["LCI"], 2),
      WQI_UCI  = round(stats_m["UCI"], 2),
      stringsAsFactors = FALSE
    )
  })

  summary_df <- do.call(rbind, summary_rows)

  summary_df$Status <- cut(
    summary_df$WQI_Mean,
    breaks = c(-Inf, 25, 50, 75, 100, Inf),
    labels = c("Excellent", "Good", "Poor", "Very Poor", "Unsuitable"),
    right = TRUE
  )

  # 7. Two-Sheet Excel / CSV Export
  if (!is.null(file)) {
    ext <- tolower(tools::file_ext(file))
    if (ext == "xlsx") {
      if (!requireNamespace("writexl", quietly = TRUE)) {
        stop("Please install 'writexl' to export to Excel.")
      }
      writexl::write_xlsx(
        list(
          Replicate_Scores = raw_out,
          Station_Summary  = summary_df
        ),
        path = file
      )
      message(sprintf("WQI replicate and summary data exported to: %s", file))
    } else if (ext == "csv") {
      utils::write.csv(summary_df, file = file, row.names = FALSE)
      message(sprintf("Station summary exported to CSV: %s", file))
    }
  }

  return(list(
    Replicate_Data = raw_out,
    Summary_Data = summary_df
  ))
}
