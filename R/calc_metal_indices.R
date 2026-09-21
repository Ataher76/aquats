#' Heavy Metal Pollution Indices with Confidence Intervals
#'
#' Computes five heavy metal pollution indices (HPI, HEI, Cd, PLI, and Nemerow PN)
#' per observation or replicate. If replicated sampling is detected, it calculates
#' station-level summary statistics including Mean, SD, SE, and 95% Confidence Intervals.
#'
#' @param data A \code{data.frame}, \code{matrix}, or file path (.csv, .xlsx).
#' @param id_col Character. Column identifying stations/sites (default: \code{"Station"}).
#' @param standards Named numeric vector of permissible guideline limits (Si) (ug/L or mg/L).
#' @param ci Numeric. Confidence level for interval estimation (default: 0.95).
#' @param file Optional character. File path to export results (.xlsx or .csv).
#'
#' @return A list containing:
#'   \item{Replicate_Data}{Data frame containing computed indices for each individual sample.}
#'   \item{Summary_Data}{Station-level summary statistics (Mean, SD, SE, 95\% CI, Status).}
#' @export
#'
#' @importFrom stats qt sd na.omit
#' @importFrom tools file_ext
#' @importFrom utils read.csv write.csv
calc_metal_indices <- function(
    data,
    id_col = "Station",
    standards = c(
      As = 10,
      Cd = 3,
      Cr = 50,
      Cu = 2000,
      Pb = 10,
      Zn = 3000
    ),
    ci = 0.95,
    file = NULL
) {

  # 1. Standard Defensive Ingestion Guard & File Support
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

  # 2. Match Target Metals
  metal_names <- names(standards)
  matched_metals <- metal_names[metal_names %in% names(df)]

  if (length(matched_metals) < 2) {
    stop(sprintf(
      "At least 2 matching metals are required. Detected: %s",
      paste(matched_metals, collapse = ", ")
    ))
  }

  std_active <- standards[matched_metals]
  n_metals <- length(matched_metals)
  n_rows <- nrow(df)

  conc_mat <- as.matrix(df[, matched_metals, drop = FALSE])
  mode(conc_mat) <- "numeric"

  # 3. Compute Individual Contamination Factors (CF_i = C_i / S_i)
  cf_mat <- sweep(conc_mat, MARGIN = 2, STATS = std_active, FUN = "/")
  colnames(cf_mat) <- paste0("CF_", matched_metals)

  # 4. Compute Core Indices per Replicate
  hei_vals <- rowSums(cf_mat, na.rm = TRUE)
  cd_vals <- rowSums(cf_mat - 1, na.rm = TRUE)

  pli_vals <- apply(cf_mat, 1, function(r) {
    r_clean <- r[!is.na(r) & r > 0]
    if (length(r_clean) == 0) return(NA_real_)
    prod(r_clean)^(1 / length(r_clean))
  })

  cf_mean <- rowMeans(cf_mat, na.rm = TRUE)
  cf_max <- apply(cf_mat, 1, max, na.rm = TRUE)
  pn_vals <- sqrt((cf_mean^2 + cf_max^2) / 2)

  # HPI calculation
  k_factor <- 1 / sum(1 / std_active)
  unit_weights <- k_factor / std_active
  qi_mat <- matrix(NA_real_, nrow = n_rows, ncol = n_metals)
  for (j in seq_len(n_metals)) {
    qi_mat[, j] <- (abs(conc_mat[, j]) / std_active[j]) * 100
  }

  hpi_vals <- apply(qi_mat, 1, function(row_q) {
    if (any(is.na(row_q))) return(NA_real_)
    sum(row_q * unit_weights) / sum(unit_weights)
  })

  # Assemble Replicate Master Frame
  raw_out <- data.frame(
    Station = df[[id_col]],
    round(conc_mat, 2),
    round(cf_mat, 2),
    HPI = round(hpi_vals, 2),
    HEI = round(hei_vals, 2),
    Cd = round(cd_vals, 2),
    PLI = round(pli_vals, 2),
    Nemerow_PN = round(pn_vals, 2),
    stringsAsFactors = FALSE
  )

  # 5. Station-Level Summary with 95% Confidence Intervals
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

  metrics <- c("HPI", "HEI", "Cd", "PLI", "Nemerow_PN")
  split_indices <- split(raw_out[, metrics], raw_out$Station)

  summary_rows <- lapply(names(split_indices), function(stn) {
    sub_df <- split_indices[[stn]]
    res_list <- list(Station = stn, N = nrow(sub_df))

    for (m in metrics) {
      stats_m <- compute_ci(sub_df[[m]])
      res_list[[paste0(m, "_Mean")]] <- round(stats_m["Mean"], 2)
      res_list[[paste0(m, "_SD")]]   <- round(stats_m["SD"], 2)
      res_list[[paste0(m, "_SE")]]   <- round(stats_m["SE"], 2)
      res_list[[paste0(m, "_LCI")]]  <- round(stats_m["LCI"], 2)
      res_list[[paste0(m, "_UCI")]]  <- round(stats_m["UCI"], 2)
    }

    as.data.frame(res_list, stringsAsFactors = FALSE)
  })

  summary_df <- do.call(rbind, summary_rows)

  # Assign Standard Status Tiers based on Station Mean Values
  summary_df$HPI_Status <- cut(
    summary_df$HPI_Mean,
    breaks = c(-Inf, 100, Inf),
    labels = c("Suitable", "Polluted/Critical"),
    right = TRUE
  )
  summary_df$HEI_Status <- cut(
    summary_df$HEI_Mean,
    breaks = c(-Inf, 10, 20, Inf),
    labels = c("Low", "Medium", "High"),
    right = FALSE
  )
  summary_df$Cd_Status <- cut(
    summary_df$Cd_Mean,
    breaks = c(-Inf, 6, 12, Inf),
    labels = c("Low Contamination", "Moderate Contamination", "High Contamination"),
    right = FALSE
  )
  summary_df$PLI_Status <- cut(
    summary_df$PLI_Mean,
    breaks = c(-Inf, 1, Inf),
    labels = c("Unpolluted", "Polluted"),
    right = TRUE
  )
  summary_df$Nemerow_Status <- cut(
    summary_df$Nemerow_PN_Mean,
    breaks = c(-Inf, 0.7, 1.0, 2.0, 3.0, Inf),
    labels = c("Clean", "Warning Limit", "Slightly Polluted", "Moderately Polluted", "Heavily Polluted"),
    right = TRUE
  )

  # 6. Two-Sheet Excel / CSV Export
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
      message(sprintf("Two-sheet pollution indices exported to: %s", file))
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
