#' Integrated Environmental Correlation and Multi-Community Mantel Network
#'
#' Computes pairwise environmental correlations (with significance testing) while
#' performing simple or partial Mantel tests across single or multiple biological
#' communities (e.g., Phytoplankton, Zooplankton, Benthos, Fish). Generates
#' publication-ready lower or upper triangular matrices coupled with smooth
#' Bezier or straight multi-community network linkages.
#'
#' @param comm_data A \code{data.frame}, \code{matrix}, or a \code{named list} of
#'   data frames representing distinct biological communities.
#' @param env_data A \code{data.frame} or \code{matrix} of numeric environmental variables.
#' @param plot_type Character string; visualization layout: \code{"network"} (default)
#'   or \code{"bars"} (synchronized side-by-side barplot). Default is \code{"network"}.
#' @param matrix_type Character string; triangle orientation for network mode:
#'   \code{"lower"} (default, matrix at bottom-left, communities at top-right) or
#'   \code{"upper"} (matrix at top-right, communities at bottom-left).
#' @param cell_shape Character string; geometry of correlation cells:
#'   \code{"tile"} (default full square), \code{"circle"} (proportional circle scaled by |r|),
#'   or \code{"square"} (proportional square scaled by |r|).
#' @param line_style Character string; linkage geometry: \code{"curve"} (default) or \code{"straight"}.
#' @param mantel_type Character string; \code{"simple"} (default standard Mantel test via
#'   \code{vegan::mantel}) or \code{"partial"} (partial Mantel test via \code{vegan::mantel.partial}
#'   controlling for all remaining environmental variables).
#' @param method Character string; correlation method for environmental variables:
#'   \code{"pearson"} (default), \code{"spearman"}, or \code{"kendall"}.
#' @param mantel_method Optional character string; correlation method for the Mantel test
#'   (\code{"pearson"}, \code{"spearman"}, or \code{"kendall"}). If \code{NULL} (default),
#'   inherits the value of \code{method}.
#' @param spec_dist Character string; dissimilarity metric for community abundance data
#'   via \code{vegan::vegdist} (e.g., \code{"bray"}, \code{"jaccard"}). Default is \code{"bray"}.
#' @param env_dist Character string; distance metric for environmental parameters. Default is \code{"euclidean"}.
#' @param transform Character string; transformation applied to community counts via
#'   \code{vegan::decostand}: \code{"none"} (default), \code{"hellinger"}, \code{"log"},
#'   \code{"pa"}, or \code{"sqrt"}.
#' @param p_adjust Character string; multiple-testing correction method for Mantel p-values
#'   passed to \code{stats::p.adjust} (e.g., \code{"none"}, \code{"BH"}, \code{"holm"},
#'   \code{"bonferroni"}). Default is \code{"none"}.
#' @param r_breaks Numeric vector of length 2 defining the three Mantel r tiers.
#'   Default is \code{c(0.20, 0.40)}.
#' @param p_breaks Numeric vector of length 2 defining the three significance tiers.
#'   Default is \code{c(0.01, 0.05)}.
#' @param only_significant Logical; if \code{TRUE}, displays only statistically significant
#'   Mantel linkages (\eqn{p < p\_breaks[2]}). Default is \code{FALSE}.
#' @param hide_ns_lines Logical; if \code{TRUE}, renders non-significant linkages with low
#'   opacity so significant drivers stand out clearly. Default is \code{TRUE}.
#' @param line_color_by Character string; variable mapping for linkage line colors:
#'   \code{"significance"} (default) or \code{"community"}.
#' @param curviness Numeric; base curvature multiplier for network arcs. Default is \code{0.20}.
#' @param permutations Numeric integer; number of Monte Carlo permutations for Mantel testing.
#'   Default is \code{999}.
#' @param seed Optional integer; random seed for reproducible permutations. Default is \code{42}.
#' @param show_cor_text Logical; if \code{TRUE}, displays correlation coefficients inside cells.
#'   Default is \code{TRUE}.
#' @param show_cor_sig Logical; if \code{TRUE}, appends significance stars (\code{*}, \code{**},
#'   \code{***}) to correlation values in the heatmap. Default is \code{FALSE}.
#' @param diag_label Character string; placement of variable names along the diagonal:
#'   \code{"outside"} (default, collision-free beside anchor pins), \code{"inside"} (inside
#'   diagonal tiles), or \code{"axis"} (only on outer x/y axes).
#' @param color_palette Character string of a diverging \code{RColorBrewer} palette (default \code{"RdBu"})
#'   or a character vector of 3 custom colors (low, mid, high) for the correlation heatmap.
#' @param sig_colors Optional character vector of 3 colors for the Mantel significance
#'   tiers (strong, moderate, non-significant).
#' @param title Optional character string; master figure title. Default is \code{NULL}.
#'
#' @return A list containing:
#'   \item{Plot}{The primary publication-ready \code{ggplot2} object.}
#'   \item{Correlation_Matrix}{Symmetric matrix of environmental correlation coefficients.}
#'   \item{Correlation_P_Matrix}{Symmetric matrix of pairwise environmental correlation p-values.}
#'   \item{Mantel_Results}{Tidy data frame containing community-specific Mantel \eqn{r} values, raw and adjusted \eqn{p}-values, and significance tiers.}
#' @export
#'
#' @import ggplot2
#' @importFrom vegan vegdist mantel mantel.partial decostand
#' @importFrom stats cor cor.test complete.cases sd p.adjust setNames
#' @importFrom tools toTitleCase
#' @importFrom rlang .data
mantel_heatmap_analysis <- function(
    comm_data,
    env_data,
    plot_type = c("network", "bars"),
    matrix_type = c("lower", "upper"),
    cell_shape = c("tile", "circle", "square"),
    line_style = c("curve", "straight"),
    mantel_type = c("simple", "partial"),
    method = c("pearson", "spearman", "kendall"),
    mantel_method = NULL,
    spec_dist = "bray",
    env_dist = "euclidean",
    transform = c("none", "hellinger", "log", "pa", "sqrt"),
    p_adjust = c("none", "BH", "holm", "bonferroni"),
    r_breaks = c(0.20, 0.40),
    p_breaks = c(0.01, 0.05),
    only_significant = FALSE,
    hide_ns_lines = TRUE,
    line_color_by = c("significance", "community"),
    curviness = 0.20,
    permutations = 999,
    seed = 42,
    show_cor_text = TRUE,
    show_cor_sig = FALSE,
    diag_label = c("outside", "inside", "axis"),
    color_palette = "RdBu",
    sig_colors = NULL,
    title = NULL
) {

  # 1. DEFENSIVE INGESTION: ENVIRONMENTAL DATA
  if (is.matrix(env_data)) {
    env_data <- as.data.frame(env_data)
  } else if (!is.data.frame(env_data)) {
    stop("Input 'env_data' must be a data.frame, tibble, or matrix.")
  }

  # DEFENSIVE INGESTION: COMMUNITY DATA
  if (is.matrix(comm_data) || is.data.frame(comm_data)) {
    comm_list <- list("Aquatic Community" = as.data.frame(comm_data))
  } else if (is.list(comm_data)) {
    if (is.null(names(comm_data)) || any(names(comm_data) == "")) {
      names(comm_data) <- paste0("Community_", seq_along(comm_data))
    }
    comm_list <- lapply(comm_data, as.data.frame)
  } else {
    stop("Input 'comm_data' must be a data.frame, matrix, or a named list of data frames/matrices.")
  }

  plot_type     <- match.arg(plot_type)
  matrix_type   <- match.arg(matrix_type)
  cell_shape    <- match.arg(cell_shape)
  line_style    <- match.arg(line_style)
  mantel_type   <- match.arg(mantel_type)
  method        <- match.arg(method)
  transform     <- match.arg(transform)
  p_adjust      <- match.arg(p_adjust)
  line_color_by <- match.arg(line_color_by)
  diag_label    <- match.arg(diag_label)

  if (is.null(mantel_method)) {
    mantel_method <- method
  } else {
    mantel_method <- match.arg(mantel_method, c("pearson", "spearman", "kendall"))
  }

  if (length(r_breaks) != 2 || r_breaks[1] >= r_breaks[2]) {
    stop("'r_breaks' must be a numeric vector of length 2 in ascending order (e.g., c(0.20, 0.40)).")
  }
  if (length(p_breaks) != 2 || p_breaks[1] >= p_breaks[2]) {
    stop("'p_breaks' must be a numeric vector of length 2 in ascending order (e.g., c(0.01, 0.05)).")
  }

  if (!is.null(seed)) set.seed(seed)

  env_num <- env_data[, vapply(env_data, is.numeric, logical(1)), drop = FALSE]
  if (ncol(env_num) < 2) {
    stop("Environmental matrix must contain at least 2 numeric parameter columns.")
  }

  n_rows_env <- nrow(env_num)
  for (c_name in names(comm_list)) {
    c_mat <- comm_list[[c_name]]
    c_num <- c_mat[, vapply(c_mat, is.numeric, logical(1)), drop = FALSE]
    if (nrow(c_num) != n_rows_env) {
      stop(sprintf(
        "Row count mismatch: Community '%s' has %d rows, but 'env_data' has %d rows.",
        c_name, nrow(c_num), n_rows_env
      ))
    }
    comm_list[[c_name]] <- c_num
  }

  zero_var <- names(env_num)[vapply(env_num, function(x) stats::sd(x, na.rm = TRUE) == 0, logical(1))]
  if (length(zero_var) > 0) {
    stop(sprintf("Zero variance detected in environmental parameter(s): %s.", paste(zero_var, collapse = ", ")))
  }

  # 2. ENVIRONMENTAL CORRELATION & P-VALUE MATRICES
  env_vars <- colnames(env_num)
  n_env    <- length(env_vars)
  cor_res  <- stats::cor(env_num, method = method, use = "pairwise.complete.obs")

  cor_p <- matrix(0, nrow = n_env, ncol = n_env, dimnames = list(env_vars, env_vars))
  for (i in seq_len(n_env)) {
    for (j in seq_len(i)) {
      if (i != j) {
        ct <- suppressWarnings(stats::cor.test(env_num[[i]], env_num[[j]], method = method))
        cor_p[i, j] <- ct$p.value
        cor_p[j, i] <- ct$p.value
      }
    }
  }

  # 3. ITERATIVE COMMUNITY MANTEL / PARTIAL MANTEL TESTS
  mantel_rows <- list()

  for (c_name in names(comm_list)) {
    c_mat <- comm_list[[c_name]]

    valid_idx <- stats::complete.cases(c_mat) &
      stats::complete.cases(env_num) &
      (rowSums(c_mat, na.rm = TRUE) > 0)
    sub_comm  <- c_mat[valid_idx, , drop = FALSE]
    sub_env   <- env_num[valid_idx, , drop = FALSE]

    if (nrow(sub_comm) < 3) {
      warning(sprintf("Community '%s' has fewer than 3 valid complete samples. Skipping.", c_name))
      next
    }

    if (transform == "hellinger") {
      sub_comm <- vegan::decostand(sub_comm, method = "hellinger")
    } else if (transform == "log") {
      sub_comm <- vegan::decostand(sub_comm, method = "log")
    } else if (transform == "pa") {
      sub_comm <- vegan::decostand(sub_comm, method = "pa")
    } else if (transform == "sqrt") {
      sub_comm <- sqrt(sub_comm)
    }

    comm_d   <- vegan::vegdist(sub_comm, method = spec_dist, binary = (transform == "pa"))
    scaled_e <- as.data.frame(scale(sub_env))

    for (var_name in env_vars) {
      single_var <- scaled_e[, var_name, drop = FALSE]
      env_d      <- vegan::vegdist(single_var, method = env_dist)

      if (mantel_type == "partial" && n_env > 2) {
        other_vars <- scaled_e[, setdiff(env_vars, var_name), drop = FALSE]
        covar_d    <- vegan::vegdist(other_vars, method = env_dist)
        m_test     <- vegan::mantel.partial(
          comm_d, env_d, covar_d,
          method = mantel_method,
          permutations = permutations
        )
      } else {
        m_test <- vegan::mantel(
          comm_d, env_d,
          method = mantel_method,
          permutations = permutations
        )
      }

      mantel_rows[[length(mantel_rows) + 1]] <- data.frame(
        Community = c_name,
        Variable  = var_name,
        Mantel_r  = round(m_test$statistic, 3),
        p_raw     = round(m_test$signif, 4),
        stringsAsFactors = FALSE
      )
    }
  }

  mantel_results <- do.call(rbind, mantel_rows)
  mantel_results$p_value <- round(stats::p.adjust(mantel_results$p_raw, method = p_adjust), 4)

  p_labels <- c(
    sprintf("p < %g", p_breaks[1]),
    sprintf("%g \u2264 p < %g", p_breaks[1], p_breaks[2]),
    sprintf("p \u2265 %g", p_breaks[2])
  )
  mantel_results$Significance <- cut(
    mantel_results$p_value,
    breaks = c(-Inf, p_breaks[1], p_breaks[2], Inf),
    labels = p_labels,
    right = FALSE
  )
  mantel_results$Significance <- factor(mantel_results$Significance, levels = p_labels)

  r_labels <- c(
    sprintf("< %.2f", r_breaks[1]),
    sprintf("%.2f - %.2f", r_breaks[1], r_breaks[2]),
    sprintf("\u2265 %.2f", r_breaks[2])
  )
  mantel_results$Mantel_R_Tier <- cut(
    abs(mantel_results$Mantel_r),
    breaks = c(-Inf, r_breaks[1], r_breaks[2], Inf),
    labels = r_labels,
    right = FALSE
  )
  mantel_results$Mantel_R_Tier <- factor(mantel_results$Mantel_R_Tier, levels = r_labels)
  mantel_results$Line_Type <- ifelse(mantel_results$p_value < p_breaks[2], "solid", "dashed")

  if (is.null(sig_colors) || length(sig_colors) < 3) {
    sig_pal <- c("#D9534F", "#2B8CBE", "grey78")
  } else {
    sig_pal <- unname(sig_colors[1:3])
  }
  names(sig_pal) <- p_labels

  add_corr_fill_scale <- function(label_name) {
    if (length(color_palette) >= 3) {
      ggplot2::scale_fill_gradient2(
        low = color_palette[1],
        mid = color_palette[2],
        high = color_palette[3],
        midpoint = 0,
        limit = c(-1, 1),
        name = label_name,
        guide = ggplot2::guide_colorbar(barwidth = ggplot2::unit(0.45, "cm"), barheight = ggplot2::unit(2.6, "cm"))
      )
    } else {
      ggplot2::scale_fill_distiller(
        palette = color_palette[1],
        limit = c(-1, 1),
        direction = -1,
        name = label_name,
        guide = ggplot2::guide_colorbar(barwidth = ggplot2::unit(0.45, "cm"), barheight = ggplot2::unit(2.6, "cm"))
      )
    }
  }

  # =========================================================================
  # PLOT STYLE A: NETWORK
  # =========================================================================
  if (plot_type == "network") {

    tile_list <- list()
    for (i in seq_len(n_env)) {
      for (j in seq_len(n_env)) {
        keep_cell <- if (matrix_type == "lower") (j <= i) else (j >= i)
        if (keep_cell) {
          r_val <- cor_res[env_vars[i], env_vars[j]]
          p_val <- cor_p[env_vars[i], env_vars[j]]
          stars <- if (p_val < 0.001) "***" else if (p_val < 0.01) "**" else if (p_val < 0.05) "*" else ""
          lbl   <- if (show_cor_sig) sprintf("%.2f%s", r_val, stars) else sprintf("%.2f", r_val)

          tile_list[[length(tile_list) + 1]] <- data.frame(
            x = j,
            y = n_env - i + 1,
            Var1 = env_vars[j],
            Var2 = env_vars[i],
            Correlation = r_val,
            Abs_Corr = abs(r_val),
            Label = lbl,
            is_diag = (i == j),
            stringsAsFactors = FALSE
          )
        }
      }
    }
    tiles_df <- do.call(rbind, tile_list)
    off_diag <- tiles_df[!tiles_df$is_diag, , drop = FALSE]
    on_diag  <- tiles_df[tiles_df$is_diag, , drop = FALSE]

    diag_pts <- data.frame(
      x = seq_len(n_env),
      y = n_env - seq_len(n_env) + 1,
      Variable = env_vars,
      stringsAsFactors = FALSE
    )

    pin_offset <- if (matrix_type == "lower") 0.42 else -0.42
    diag_pts$pin_x <- diag_pts$x + pin_offset
    diag_pts$pin_y <- diag_pts$y + pin_offset

    u_comms <- unique(mantel_results$Community)
    n_comms <- length(u_comms)

    if (matrix_type == "lower") {
      comm_x_pos <- n_env + max(1.75, n_env * 0.30)
      comm_y_pos <- if (n_comms == 1) n_env * 0.65 else seq(n_env * 0.86, n_env * 0.20, length.out = n_comms)
    } else {
      comm_x_pos <- 1 - max(1.85, n_env * 0.32)
      comm_y_pos <- if (n_comms == 1) n_env * 0.35 else seq(n_env * 0.70, n_env * 0.12, length.out = n_comms)
    }

    node_palette <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#3A86FF")
    assigned_colors <- node_palette[((seq_len(n_comms) - 1) %% length(node_palette)) + 1]

    comm_nodes <- data.frame(
      Community  = u_comms,
      comm_x     = comm_x_pos,
      comm_y     = comm_y_pos,
      node_color = assigned_colors,
      stringsAsFactors = FALSE
    )

    arc_df <- merge(mantel_results, diag_pts, by = "Variable")
    arc_df <- merge(arc_df, comm_nodes, by = "Community")

    arc_df$x_start <- arc_df$comm_x
    arc_df$y_start <- arc_df$comm_y
    arc_df$x_end   <- arc_df$pin_x
    arc_df$y_end   <- arc_df$pin_y

    # Fan-graduated negative curvature bows arcs away from both hubs and diagonal pins
    if (line_style == "curve") {
      if (matrix_type == "lower") {
        arc_df$curvature <- -curviness * (0.85 + 0.04 * (n_env - arc_df$x))
      } else {
        arc_df$curvature <- -curviness * (0.72 + 0.06 * arc_df$x)
      }
    } else {
      arc_df$curvature <- 0
    }

    if (only_significant) {
      arc_df <- arc_df[arc_df$p_value < p_breaks[2], , drop = FALSE]
    }

    arc_df <- arc_df[order(arc_df$p_value, decreasing = TRUE), , drop = FALSE]
    arc_df$Alpha <- ifelse(
      arc_df$p_value < p_breaks[2],
      0.88,
      if (hide_ns_lines) 0.32 else 0.75
    )

    build_bezier_paths <- function(df, n_pts = 40) {
      t_seq <- seq(0, 1, length.out = n_pts)
      path_list <- vector("list", nrow(df))
      for (k in seq_len(nrow(df))) {
        r <- df[k, ]
        dx <- r$x_end - r$x_start
        dy <- r$y_end - r$y_start
        mid_x <- (r$x_start + r$x_end) / 2
        mid_y <- (r$y_start + r$y_end) / 2
        ctrl_x <- mid_x - r$curvature * dy
        ctrl_y <- mid_y + r$curvature * dx

        bx <- (1 - t_seq)^2 * r$x_start + 2 * (1 - t_seq) * t_seq * ctrl_x + t_seq^2 * r$x_end
        by <- (1 - t_seq)^2 * r$y_start + 2 * (1 - t_seq) * t_seq * ctrl_y + t_seq^2 * r$y_end

        path_list[[k]] <- data.frame(
          edge_id       = k,
          x             = bx,
          y             = by,
          Community     = r$Community,
          Significance  = r$Significance,
          Mantel_R_Tier = r$Mantel_R_Tier,
          Line_Type     = r$Line_Type,
          Alpha         = r$Alpha,
          stringsAsFactors = FALSE
        )
      }
      out <- do.call(rbind, path_list)
      out$Significance  <- factor(out$Significance, levels = p_labels)
      out$Mantel_R_Tier <- factor(out$Mantel_R_Tier, levels = r_labels)
      out
    }

    p <- ggplot2::ggplot()

    # Base grid cells
    if (cell_shape == "tile") {
      p <- p + ggplot2::geom_tile(
        data = off_diag,
        ggplot2::aes(x = .data$x, y = .data$y, fill = .data$Correlation),
        color = "white", linewidth = 0.75
      )
    } else {
      min_pt <- if (show_cor_text) 5.8 else 3.2
      p <- p +
        ggplot2::geom_tile(
          data = off_diag,
          ggplot2::aes(x = .data$x, y = .data$y),
          fill = "grey98", color = "grey88", linewidth = 0.45
        ) +
        ggplot2::geom_point(
          data = off_diag,
          ggplot2::aes(x = .data$x, y = .data$y, fill = .data$Correlation, size = .data$Abs_Corr),
          shape = if (cell_shape == "circle") 21 else 22,
          color = if (show_cor_text) "grey55" else "grey30",
          stroke = 0.35
        ) +
        ggplot2::scale_size_continuous(range = c(min_pt, 11.2), guide = "none")
    }

    # Correlation text inside off-diagonal cells
    if (show_cor_text && nrow(off_diag) > 0) {
      dark_tiles  <- off_diag[off_diag$Abs_Corr > 0.62, , drop = FALSE]
      light_tiles <- off_diag[off_diag$Abs_Corr <= 0.62, , drop = FALSE]

      txt_size <- if (show_cor_sig) 2.75 else 3.10
      if (nrow(dark_tiles) > 0) {
        p <- p + ggplot2::geom_text(
          data = dark_tiles,
          ggplot2::aes(x = .data$x, y = .data$y, label = .data$Label),
          color = "white", fontface = "bold", size = txt_size
        )
      }
      if (nrow(light_tiles) > 0) {
        p <- p + ggplot2::geom_text(
          data = light_tiles,
          ggplot2::aes(x = .data$x, y = .data$y, label = .data$Label),
          color = "grey15", fontface = "bold", size = txt_size
        )
      }
    }

    # Diagonal & outer axis labels
    if (diag_label == "inside") {
      p <- p +
        ggplot2::geom_tile(
          data = on_diag,
          ggplot2::aes(x = .data$x, y = .data$y),
          fill = "grey93", color = "white", linewidth = 0.75
        ) +
        ggplot2::geom_text(
          data = on_diag,
          ggplot2::aes(x = .data$x, y = .data$y, label = .data$Var1),
          fontface = "bold", color = "grey15", size = 3.4
        )
    } else if (diag_label == "outside") {
      if (matrix_type == "lower") {
        x_ax <- data.frame(x = seq_len(n_env - 1), y = 0.38, label = env_vars[seq_len(n_env - 1)])
        y_ax <- data.frame(x = 0.38, y = n_env - (2:n_env) + 1, label = env_vars[2:n_env])
        p <- p +
          ggplot2::geom_text(data = y_ax, ggplot2::aes(x = .data$x, y = .data$y, label = .data$label), hjust = 1, fontface = "bold", size = 3.4, color = "grey15") +
          ggplot2::geom_text(data = x_ax, ggplot2::aes(x = .data$x, y = .data$y, label = .data$label), angle = 45, hjust = 1, vjust = 1, fontface = "bold", size = 3.4, color = "grey15") +
          ggplot2::geom_label(
            data = diag_pts,
            ggplot2::aes(x = .data$x + 0.02, y = .data$y + 0.02, label = .data$Variable),
            fontface = "bold", size = 3.1, color = "grey20", fill = "white", alpha = 0.90, label.size = 0, label.padding = ggplot2::unit(0.12, "lines")
          )
      } else {
        x_ax <- data.frame(x = 2:n_env, y = n_env + 0.62, label = env_vars[2:n_env])
        y_ax <- data.frame(x = n_env + 0.62, y = n_env - seq_len(n_env - 1) + 1, label = env_vars[seq_len(n_env - 1)])
        p <- p +
          ggplot2::geom_text(data = y_ax, ggplot2::aes(x = .data$x, y = .data$y, label = .data$label), hjust = 0, fontface = "bold", size = 3.4, color = "grey15") +
          ggplot2::geom_text(data = x_ax, ggplot2::aes(x = .data$x, y = .data$y, label = .data$label), angle = 45, hjust = 0, vjust = 0, fontface = "bold", size = 3.4, color = "grey15") +
          ggplot2::geom_label(
            data = diag_pts,
            ggplot2::aes(x = .data$x - 0.02, y = .data$y - 0.02, label = .data$Variable),
            fontface = "bold", size = 3.1, color = "grey20", fill = "white", alpha = 0.90, label.size = 0, label.padding = ggplot2::unit(0.12, "lines")
          )
      }
    } else {
      if (matrix_type == "lower") {
        x_ax <- data.frame(x = seq_len(n_env - 1), y = 0.38, label = env_vars[seq_len(n_env - 1)])
        y_ax <- data.frame(x = 0.38, y = n_env - (2:n_env) + 1, label = env_vars[2:n_env])
      } else {
        x_ax <- data.frame(x = 2:n_env, y = n_env + 0.62, label = env_vars[2:n_env])
        y_ax <- data.frame(x = n_env + 0.62, y = n_env - seq_len(n_env - 1) + 1, label = env_vars[seq_len(n_env - 1)])
      }
      p <- p +
        ggplot2::geom_text(data = y_ax, ggplot2::aes(x = .data$x, y = .data$y, label = .data$label), hjust = if (matrix_type == "lower") 1 else 0, fontface = "bold", size = 3.5, color = "grey15") +
        ggplot2::geom_text(data = x_ax, ggplot2::aes(x = .data$x, y = .data$y, label = .data$label), angle = 45, hjust = if (matrix_type == "lower") 1 else 0, fontface = "bold", size = 3.5, color = "grey15")
    }

    # Single vectorized geom_path layer
    if (nrow(arc_df) > 0) {
      bezier_df <- build_bezier_paths(arc_df, n_pts = if (line_style == "curve") 40 else 2)

      if (line_color_by == "significance") {
        p <- p +
          ggplot2::geom_path(
            data = bezier_df,
            ggplot2::aes(
              x = .data$x, y = .data$y,
              group = .data$edge_id,
              color = .data$Significance,
              linewidth = .data$Mantel_R_Tier,
              linetype = .data$Line_Type,
              alpha = .data$Alpha
            ),
            lineend = "round",
            show.legend = TRUE
          ) +
          ggplot2::scale_color_manual(
            values = sig_pal,
            limits = p_labels,
            name = if (p_adjust == "none") "Mantel p-value" else sprintf("Mantel p (%s)", p_adjust),
            drop = FALSE,
            guide = ggplot2::guide_legend(override.aes = list(linewidth = 1.6, alpha = 1, linetype = "solid"))
          )
      } else {
        comm_pal <- stats::setNames(comm_nodes$node_color, comm_nodes$Community)
        p <- p +
          ggplot2::geom_path(
            data = bezier_df,
            ggplot2::aes(
              x = .data$x, y = .data$y,
              group = .data$edge_id,
              color = .data$Community,
              linewidth = .data$Mantel_R_Tier,
              linetype = .data$Line_Type,
              alpha = .data$Alpha
            ),
            lineend = "round",
            show.legend = TRUE
          ) +
          ggplot2::scale_color_manual(
            values = comm_pal,
            limits = comm_nodes$Community,
            name = "Community",
            drop = FALSE,
            guide = ggplot2::guide_legend(override.aes = list(linewidth = 1.6, alpha = 1, linetype = "solid"))
          )
      }
      p <- p + ggplot2::scale_alpha_identity()
    }

    # Diagonal anchor pins (drawn ON TOP of curve endpoints)
    p <- p +
      ggplot2::geom_point(
        data = diag_pts,
        ggplot2::aes(x = .data$pin_x, y = .data$pin_y),
        shape = 21, size = 3.4, fill = "#2A9D8F", color = "white", stroke = 0.9
      )

    # Community hubs (drawn ON TOP of curve start points)
    lbl_nudge <- if (matrix_type == "lower") 0.28 else -0.28
    lbl_hjust <- if (matrix_type == "lower") 0 else 1

    p <- p +
      ggplot2::geom_point(
        data = comm_nodes,
        ggplot2::aes(x = .data$comm_x, y = .data$comm_y),
        shape = 21, size = 6.5, fill = comm_nodes$node_color, color = "white", stroke = 1.2
      ) +
      ggplot2::geom_point(
        data = comm_nodes,
        ggplot2::aes(x = .data$comm_x, y = .data$comm_y),
        size = 2.0, color = "white"
      ) +
      ggplot2::geom_text(
        data = comm_nodes,
        ggplot2::aes(x = .data$comm_x + lbl_nudge, y = .data$comm_y, label = .data$Community),
        fontface = "bold", size = 3.8, hjust = lbl_hjust, color = "grey10"
      )

    r_widths <- stats::setNames(c(0.55, 1.45, 2.75), r_labels)

    max_comm_chars <- max(nchar(as.character(u_comms)), na.rm = TRUE)
    label_pad      <- max(2.2, max_comm_chars * 0.17)
    max_var_chars  <- max(nchar(as.character(env_vars)), na.rm = TRUE)
    axis_pad       <- max(0.5, max_var_chars * 0.12)

    x_limits <- if (matrix_type == "lower") {
      c(1 - axis_pad - 0.5, comm_x_pos + label_pad)
    } else {
      c(comm_x_pos - label_pad, n_env + axis_pad + 0.8)
    }
    y_limits <- if (matrix_type == "lower") {
      c(1 - axis_pad - 0.5, n_env + 0.8)
    } else {
      c(-0.3, n_env + axis_pad + 0.9)
    }

    r_legend_title <- if (mantel_type == "partial") "Partial Mantel's |r|" else "Mantel's |r|"

    p <- p +
      add_corr_fill_scale(paste0("Corr (", tools::toTitleCase(method), ")")) +
      ggplot2::scale_linewidth_manual(
        values = r_widths,
        limits = r_labels,
        name = r_legend_title,
        drop = FALSE,
        guide = ggplot2::guide_legend(override.aes = list(color = "grey25", alpha = 1, linetype = "solid"))
      ) +
      ggplot2::scale_linetype_manual(
        values = c("solid" = "solid", "dashed" = "dashed"),
        guide = "none"
      ) +
      ggplot2::coord_fixed(
        xlim = x_limits,
        ylim = y_limits,
        clip = "off"
      ) +
      ggplot2::theme_void(base_size = 11) +
      ggplot2::theme(
        legend.position = "right",
        legend.box = "vertical",
        legend.spacing.y = ggplot2::unit(6, "pt"),
        legend.key.width = ggplot2::unit(1.1, "cm"),
        legend.title = ggplot2::element_text(face = "bold", size = 9.5, color = "grey15"),
        legend.text = ggplot2::element_text(size = 8.8, color = "grey20"),
        plot.title = ggplot2::element_text(face = "bold", size = 12.5, hjust = 0.4, margin = ggplot2::margin(b = 8)),
        plot.margin = ggplot2::margin(t = 12, r = 24, b = 16, l = 16)
      ) +
      ggplot2::labs(
        title = if (!is.null(title)) title else "Multi-Community Environmental Mantel Network"
      )

    final_plot <- p

    # =========================================================================
    # PLOT STYLE B: BARS
    # =========================================================================
  } else {

    cor_long <- as.data.frame(as.table(cor_res), stringsAsFactors = FALSE)
    colnames(cor_long) <- c("Var1", "Var2", "Correlation")
    cor_long$Var1 <- factor(cor_long$Var1, levels = env_vars)
    cor_long$Var2 <- factor(cor_long$Var2, levels = env_vars)
    cor_long$TextColor <- ifelse(abs(cor_long$Correlation) > 0.62, "white", "grey15")

    mantel_results$Variable <- factor(mantel_results$Variable, levels = env_vars)
    mantel_results$Stars <- ifelse(
      mantel_results$p_value < p_breaks[1], "**",
      ifelse(mantel_results$p_value < p_breaks[2], "*", "")
    )
    mantel_results$Label <- sprintf("%.2f%s", mantel_results$Mantel_r, mantel_results$Stars)

    p_heat <- ggplot2::ggplot(cor_long, ggplot2::aes(x = .data$Var1, y = .data$Var2, fill = .data$Correlation)) +
      ggplot2::geom_tile(color = "white", linewidth = 0.5) +
      add_corr_fill_scale("r") +
      ggplot2::coord_fixed() +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::labs(title = "Environmental Correlation", x = NULL, y = NULL) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, face = "bold", color = "grey15"),
        axis.text.y = ggplot2::element_text(face = "bold", color = "grey15"),
        panel.grid = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
      )

    if (show_cor_text) {
      p_heat <- p_heat +
        ggplot2::geom_text(
          ggplot2::aes(label = sprintf("%.2f", .data$Correlation), color = .data$TextColor),
          fontface = "bold", size = 3.1
        ) +
        ggplot2::scale_color_identity()
    }

    max_r <- max(c(0.25, max(mantel_results$Mantel_r, na.rm = TRUE) * 1.32))
    min_r <- min(c(0, min(mantel_results$Mantel_r, na.rm = TRUE)))
    lower_lim <- if (min_r < 0) min(c(min_r * 1.45, min_r - 0.12)) else -0.02

    p_bar <- ggplot2::ggplot(mantel_results, ggplot2::aes(x = .data$Variable, y = .data$Mantel_r, fill = .data$Significance)) +
      ggplot2::geom_col(width = 0.6, color = "grey20", linewidth = 0.35, show.legend = TRUE) +
      ggplot2::geom_hline(yintercept = 0, color = "grey20", linewidth = 0.5) +
      ggplot2::geom_text(
        ggplot2::aes(label = .data$Label, hjust = ifelse(.data$Mantel_r >= 0, -0.15, 1.15)),
        fontface = "bold", size = 3.0
      ) +
      ggplot2::coord_flip(clip = "off") +
      ggplot2::facet_wrap(~Community, ncol = 1) +
      ggplot2::scale_fill_manual(values = sig_pal, limits = p_labels, drop = FALSE) +
      ggplot2::scale_y_continuous(
        limits = c(lower_lim, max_r),
        expand = ggplot2::expansion(mult = c(0.08, 0.15))
      ) +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::labs(
        title = if (mantel_type == "partial") "Partial Mantel Linkages" else "Mantel Linkages by Community",
        x = NULL,
        y = if (mantel_type == "partial") "Partial Mantel r" else "Mantel r",
        fill = "Significance"
      ) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(face = "bold"),
        axis.text.y = ggplot2::element_text(face = "bold"),
        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
        strip.text = ggplot2::element_text(face = "bold")
      )

    if (requireNamespace("patchwork", quietly = TRUE)) {
      final_plot <- patchwork::wrap_plots(p_heat, p_bar, ncol = 2, widths = c(1.35, 1.05))
    } else if (requireNamespace("cowplot", quietly = TRUE)) {
      final_plot <- cowplot::plot_grid(p_heat, p_bar, ncol = 2, rel_widths = c(1.35, 1.05))
    } else {
      message("Note: Install 'patchwork' or 'cowplot' for side-by-side composite bar layout. Returning heatmap.")
      final_plot <- p_heat
    }
  }

  return(list(
    Plot = final_plot,
    Correlation_Matrix = cor_res,
    Correlation_P_Matrix = cor_p,
    Mantel_Results = mantel_results
  ))
}
