#' Publication-Ready Species Accumulation Curves
#'
#' Computes and plots species accumulation curves using \code{vegan::specaccum}.
#' Supports pooled curves or comparative grouped curves across habitats,
#' handles methods without standard deviations defensively, and formats
#' outputs with publication-ready styling.
#'
#' @param data A \code{data.frame}, \code{matrix}, or \code{list} containing
#'   community abundance counts with optional grouping metadata column(s).
#' @param group_col Optional character or integer vector specifying the grouping
#'   column(s) to exclude from the community matrix or to group by. Default is 1.
#' @param by_group Logical. If \code{TRUE}, generates comparative accumulation
#'   curves for each level of the first grouping column (default: \code{FALSE}).
#' @param method Character. Accumulation method passed to \code{vegan::specaccum}
#'   (default: \code{"exact"}).
#' @param palette Character. RColorBrewer palette name for grouped curves (default: \code{"Set2"}).
#' @param line_color Character. Hex code or color name for single curves (default: \code{"#2A9D8F"}).
#' @param ci Numeric. Multiplier for confidence intervals (default: 1.96).
#' @param show_ci Logical. If \code{TRUE}, displays confidence ribbons when standard
#'   deviations are available (default: \code{TRUE}).
#' @param title Optional character. Plot title.
#'
#' @return A list containing:
#'   \item{Specaccum_Object}{The raw \code{specaccum} object (or list of objects if grouped).}
#'   \item{Curve_Data}{A data frame with site counts, richness, and bounds.}
#'   \item{Plot}{A publication-ready \code{ggplot} object.}
#' @export
#'
#' @import ggplot2
#' @importFrom vegan specaccum
#' @importFrom rlang .data
specaccum_plot <- function(
    data,
    group_col = 1,
    by_group = FALSE,
    method = "exact",
    palette = "Set2",
    line_color = "#2A9D8F",
    ci = 1.96,
    show_ci = TRUE,
    title = NULL
) {

  # 1. Standard Defensive Ingestion Guard
  if (is.list(data) && !is.data.frame(data)) {
    grp_name <- names(data)[1]
    message(sprintf("Note: 'data' is a list. Analyzing first element: '%s'.", grp_name))
    df <- as.data.frame(data[[1]])
  } else if (is.matrix(data)) {
    df <- as.data.frame(data)
  } else if (is.data.frame(data)) {
    df <- as.data.frame(data)
  } else {
    stop("'data' must be a data.frame, matrix, or list.")
  }

  # 2. Extract Metadata vs Community Matrix (with auto-detection of non-numeric columns)
  non_num_cols <- names(df)[!vapply(df, is.numeric, logical(1))]

  if (!is.null(group_col)) {
    if (is.character(group_col)) {
      if (!all(group_col %in% names(df))) {
        missing_cols <- group_col[!group_col %in% names(df)]
        stop(sprintf("Columns not found: %s", paste(missing_cols, collapse = ", ")))
      }
      grp_names <- group_col
    } else {
      grp_names <- names(df)[as.integer(group_col)]
    }
    all_meta_cols <- unique(c(grp_names, non_num_cols))
    meta_df <- df[, all_meta_cols, drop = FALSE]
    comm_mat <- df[, !names(df) %in% all_meta_cols, drop = FALSE]
  } else {
    if (length(non_num_cols) > 0) {
      meta_df <- df[, non_num_cols, drop = FALSE]
      comm_mat <- df[, !names(df) %in% non_num_cols, drop = FALSE]
    } else {
      meta_df <- NULL
      comm_mat <- df
    }
    if (isTRUE(by_group)) {
      warning("'by_group = TRUE' requested without 'group_col'. Defaulting to pooled curve.")
      by_group <- FALSE
    }
  }

  # 3. Numeric Validation
  if (ncol(comm_mat) < 2) {
    stop("Community data must contain at least 2 species columns.")
  }

  if (any(is.na(comm_mat))) {
    stop("NA values detected. Replace missing values with 0 before running accumulation.")
  }

  # 4. Remove Empty Samples
  n_abund <- rowSums(comm_mat)
  valid_rows <- n_abund > 0

  if (!all(valid_rows)) {
    n_removed <- sum(!valid_rows)
    comm_mat <- comm_mat[valid_rows, , drop = FALSE]
    if (!is.null(meta_df)) meta_df <- meta_df[valid_rows, , drop = FALSE]
    message(sprintf("Note: Removed %d empty sample(s) with zero total abundance.", n_removed))
  }

  if (nrow(comm_mat) < 2) {
    stop("Species accumulation requires at least 2 valid samples.")
  }

  # 5. Accumulation Modeling
  if (!isTRUE(by_group)) {
    acc <- vegan::specaccum(comm_mat, method = method)
    has_sd <- !is.null(acc$sd) && !all(is.na(acc$sd))

    acc_df <- data.frame(
      Sites = acc$sites,
      Richness = acc$richness,
      SD = if (has_sd) acc$sd else NA_real_,
      Lower = if (has_sd) pmax(0, acc$richness - ci * acc$sd) else NA_real_,
      Upper = if (has_sd) acc$richness + ci * acc$sd else NA_real_,
      stringsAsFactors = FALSE
    )

    p <- ggplot2::ggplot(
      acc_df,
      ggplot2::aes(
        x = .data$Sites,
        y = .data$Richness
      )
    )

    if (show_ci && has_sd) {
      p <- p + ggplot2::geom_ribbon(
        ggplot2::aes(
          ymin = .data$Lower,
          ymax = .data$Upper
        ),
        fill = line_color,
        alpha = 0.2
      )
    }

    p <- p + ggplot2::geom_line(
      color = line_color,
      linewidth = 1.2
    )

    plot_obj <- acc

  } else {
    grp_var_name <- if (!is.null(meta_df)) names(meta_df)[1] else "Group"
    groups <- as.factor(meta_df[[grp_var_name]])
    split_mats <- split(comm_mat, groups)

    acc_list <- list()
    df_list <- list()

    for (grp_lvl in names(split_mats)) {
      sub_mat <- split_mats[[grp_lvl]]
      sub_valid <- rowSums(sub_mat) > 0
      sub_mat <- sub_mat[sub_valid, , drop = FALSE]

      if (nrow(sub_mat) < 2) {
        warning(sprintf("Group '%s' has fewer than 2 valid samples. Excluded from curve.", grp_lvl))
        next
      }

      sub_acc <- vegan::specaccum(sub_mat, method = method)
      sub_has_sd <- !is.null(sub_acc$sd) && !all(is.na(sub_acc$sd))

      acc_list[[grp_lvl]] <- sub_acc
      df_list[[grp_lvl]] <- data.frame(
        Group = grp_lvl,
        Sites = sub_acc$sites,
        Richness = sub_acc$richness,
        SD = if (sub_has_sd) sub_acc$sd else NA_real_,
        Lower = if (sub_has_sd) pmax(0, sub_acc$richness - ci * sub_acc$sd) else NA_real_,
        Upper = if (sub_has_sd) sub_acc$richness + ci * sub_acc$sd else NA_real_,
        stringsAsFactors = FALSE
      )
    }

    if (length(df_list) == 0) {
      stop("No groups contained at least 2 samples for accumulation analysis.")
    }

    acc_df <- do.call(rbind, df_list)
    rownames(acc_df) <- NULL

    p <- ggplot2::ggplot(
      acc_df,
      ggplot2::aes(
        x = .data$Sites,
        y = .data$Richness,
        color = .data$Group,
        fill = .data$Group
      )
    )

    has_any_sd <- any(!is.na(acc_df$SD))
    if (show_ci && has_any_sd) {
      p <- p + ggplot2::geom_ribbon(
        ggplot2::aes(
          ymin = .data$Lower,
          ymax = .data$Upper
        ),
        color = NA,
        alpha = 0.15
      )
    }

    p <- p +
      ggplot2::geom_line(linewidth = 1.1) +
      ggplot2::scale_color_brewer(palette = palette) +
      ggplot2::scale_fill_brewer(palette = palette)

    plot_obj <- acc_list
  }

  # 6. Styling & Layout
  p <- p +
    ggplot2::scale_y_continuous(
      limits = c(0, NA),
      expand = ggplot2::expansion(mult = c(0, 0.08))
    ) +
    ggplot2::scale_x_continuous(
      expand = ggplot2::expansion(mult = c(0.01, 0.04))
    ) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(
      axis.line = ggplot2::element_line(linewidth = 0.7, color = "black"),
      axis.title = ggplot2::element_text(face = "bold"),
      axis.text = ggplot2::element_text(color = "black", face = "bold"),
      legend.position = if (isTRUE(by_group)) "top" else "none",
      legend.title = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(face = "bold", size = 12, hjust = 0.5)
    ) +
    ggplot2::labs(
      x = "Number of Samples / Sampling Effort",
      y = "Cumulative Species Richness",
      color = if (isTRUE(by_group) && !is.null(meta_df)) names(meta_df)[1] else NULL,
      fill = if (isTRUE(by_group) && !is.null(meta_df)) names(meta_df)[1] else NULL,
      title = if (!is.null(title)) title else "Species Accumulation Curve"
    )

  return(list(
    Specaccum_Object = plot_obj,
    Curve_Data = acc_df,
    Plot = p
  ))
}
