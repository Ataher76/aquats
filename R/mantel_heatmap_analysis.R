#' Integrated Environmental Correlation and Multi-Community Mantel Network
#'
#' Computes pairwise environmental correlations while performing Mantel tests
#' across single or multiple biological communities (e.g., Phytoplankton,
#' Zooplankton, Benthos, Fish). Generates publication-ready network diagrams
#' with customizable straight or curved linkages and anchor pins.
#'
#' @param comm_data A \code{data.frame}, \code{matrix}, or a \code{named list} of
#'   data frames representing distinct biological communities.
#' @param env_data A \code{data.frame} or \code{matrix} of numeric environmental variables.
#' @param plot_type Character. Visualization layout: \code{"network"} (default curved/straight diagram)
#'   or \code{"bars"} (synchronized side-by-side barplot).
#' @param line_style Character. Linkage geometry: \code{"curve"} (default) or \code{"straight"}.
#' @param method Character. Correlation method: \code{"pearson"} (default) or \code{"spearman"}.
#' @param spec_dist Character. Dissimilarity metric for community data (default: \code{"bray"}).
#' @param env_dist Character. Distance metric for environmental parameters (default: \code{"euclidean"}).
#' @param transform Character. Transformation applied to community counts via \code{vegan::decostand}:
#'   \code{"none"} (default), \code{"hellinger"} (recommended for abundance), \code{"log"}, \code{"pa"}, or \code{"sqrt"}.
#' @param only_significant Logical. If \code{TRUE}, displays only significant Mantel linkages (p < 0.05).
#' @param line_color_by Character. Color linkage lines by \code{"significance"} (default) or \code{"community"}.
#' @param permutations Integer. Number of Monte Carlo permutations for Mantel tests (default: 999).
#' @param seed Optional integer. Random seed for reproducible permutations (default: 42).
#' @param show_cor_text Logical. If \code{TRUE}, displays correlation values inside tiles (default: \code{TRUE}).
#' @param color_palette Character. Diverging palette name for the correlation heatmap (default: \code{"RdBu"}).
#' @param title Optional character. Master figure title.
#'
#' @return A list containing:
#'   \item{Plot}{The primary publication-ready \code{ggplot} object.}
#'   \item{Correlation_Matrix}{Symmetric matrix of environmental correlation coefficients.}
#'   \item{Mantel_Results}{Tidy data frame containing community-specific Mantel r, p-values, and significance.}
#' @export
#'
#' @import ggplot2
#' @importFrom vegan vegdist mantel decostand
#' @importFrom stats cor complete.cases sd
#' @importFrom tools toTitleCase
#' @importFrom rlang .data
mantel_heatmap_analysis <- function(
    comm_data,
    env_data,
    plot_type = c("network", "bars"),
    line_style = c("curve", "straight"),
    method = c("pearson", "spearman"),
    spec_dist = "bray",
    env_dist = "euclidean",
    transform = c("none", "hellinger", "log", "pa", "sqrt"),
    only_significant = FALSE,
    line_color_by = c("significance", "community"),
    permutations = 999,
    seed = 42,
    show_cor_text = TRUE,
    color_palette = "RdBu",
    title = NULL
) {

  # DEFENSIVE INGESTION: ENVIRONMENTAL DATA
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
    comm_list <- lapply(comm_data, function(item) {
      if (is.matrix(item)) as.data.frame(item) else as.data.frame(item)
    })
  } else {
    stop("Input 'comm_data' must be a data.frame, matrix, or a named list of data frames/matrices.")
  }

  plot_type     <- match.arg(plot_type)
  line_style    <- match.arg(line_style)
  method        <- match.arg(method)
  transform     <- match.arg(transform)
  line_color_by <- match.arg(line_color_by)

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

  # Zero-variance check
  zero_var <- names(env_num)[vapply(env_num, function(x) stats::sd(x, na.rm = TRUE) == 0, logical(1))]
  if (length(zero_var) > 0) {
    stop(sprintf("Zero variance detected in environmental parameter(s): %s.", paste(zero_var, collapse = ", ")))
  }

  # 2. Environmental Correlation Matrix
  env_vars <- colnames(env_num)
  n_env    <- length(env_vars)
  cor_res  <- stats::cor(env_num, method = method, use = "complete.obs")

  # 3. Iterative Community Mantel Tests
  mantel_rows <- list()

  for (c_name in names(comm_list)) {
    c_mat <- comm_list[[c_name]]

    valid_idx <- stats::complete.cases(c_mat) & (rowSums(c_mat) > 0)
    sub_comm  <- c_mat[valid_idx, , drop = FALSE]
    sub_env   <- env_num[valid_idx, , drop = FALSE]

    if (nrow(sub_comm) < 3) {
      warning(sprintf("Community '%s' has fewer than 3 valid samples. Skipping.", c_name))
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

    comm_d <- vegan::vegdist(sub_comm, method = spec_dist, binary = (transform == "pa"))

    for (var_name in env_vars) {
      single_var <- sub_env[, var_name, drop = FALSE]
      if (env_dist == "euclidean") {
        single_var <- scale(single_var)
      }
      env_d  <- vegan::vegdist(single_var, method = env_dist)
      m_test <- vegan::mantel(comm_d, env_d, method = method, permutations = permutations)

      mantel_rows[[length(mantel_rows) + 1]] <- data.frame(
        Community = c_name,
        Variable  = var_name,
        Mantel_r  = round(m_test$statistic, 3),
        p_value   = round(m_test$signif, 4),
        stringsAsFactors = FALSE
      )
    }
  }

  mantel_results <- do.call(rbind, mantel_rows)

  mantel_results$Significance <- cut(
    mantel_results$p_value,
    breaks = c(-Inf, 0.01, 0.05, Inf),
    labels = c("p < 0.01", "p < 0.05", "p \u2265 0.05"),
    right = TRUE
  )
  mantel_results$Significance <- factor(mantel_results$Significance, levels = c("p < 0.01", "p < 0.05", "p \u2265 0.05"))

  mantel_results$Mantel_R_Tier <- cut(
    mantel_results$Mantel_r,
    breaks = c(-Inf, 0.20, 0.40, Inf),
    labels = c("< 0.20", "0.20 - 0.40", "\u2265 0.40"),
    right = FALSE
  )
  mantel_results$Mantel_R_Tier <- factor(mantel_results$Mantel_R_Tier, levels = c("< 0.20", "0.20 - 0.40", "\u2265 0.40"))
  mantel_results$Line_Type <- ifelse(mantel_results$p_value < 0.05, "solid", "dashed")

  # =========================================================================
  # PLOT STYLE A: NETWORK
  # =========================================================================
  if (plot_type == "network") {

    tile_list <- list()
    for (i in seq_len(n_env)) {
      for (j in seq_len(i)) {
        tile_list[[length(tile_list) + 1]] <- data.frame(
          x = j,
          y = n_env - i + 1,
          Var1 = env_vars[j],
          Var2 = env_vars[i],
          Correlation = cor_res[env_vars[i], env_vars[j]],
          is_diag = (i == j),
          stringsAsFactors = FALSE
        )
      }
    }
    tiles_df <- do.call(rbind, tile_list)

    diag_pts <- data.frame(
      x = seq_len(n_env),
      y = n_env - seq_len(n_env) + 1,
      Variable = env_vars,
      stringsAsFactors = FALSE
    )

    u_comms <- unique(mantel_results$Community)
    n_comms <- length(u_comms)
    comm_x_pos <- n_env + 1.35

    if (n_comms == 1) {
      comm_y_pos <- n_env * 0.70
    } else {
      comm_y_pos <- seq(n_env * 0.88, n_env * 0.20, length.out = n_comms)
    }

    node_palette <- c("#2A9D8F", "#E76F51", "#4575B4", "#9B5DE5", "#F4A261", "#E63946", "#3A86FF", "#8338EC")
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
    arc_df$x_end   <- arc_df$x + 0.38
    arc_df$y_end   <- arc_df$y + 0.38

    if (line_style == "curve") {
      arc_df$curvature <- ifelse(
        arc_df$y_end >= arc_df$y_start,
        -0.20 - (arc_df$y_end - arc_df$y_start) * 0.025,
        0.16 + (arc_df$y_start - arc_df$y_end) * 0.025
      )
    } else {
      arc_df$curvature <- 0
    }

    if (only_significant) {
      arc_df <- arc_df[arc_df$p_value < 0.05, , drop = FALSE]
    }

    p <- ggplot2::ggplot() +
      ggplot2::geom_tile(
        data = tiles_df[!tiles_df$is_diag, ],
        ggplot2::aes(x = .data$x, y = .data$y, fill = .data$Correlation),
        color = "white", linewidth = 0.65
      )

    if (show_cor_text) {
      dark_tiles  <- tiles_df[!tiles_df$is_diag & abs(tiles_df$Correlation) > 0.65, ]
      light_tiles <- tiles_df[!tiles_df$is_diag & abs(tiles_df$Correlation) <= 0.65, ]

      if (nrow(dark_tiles) > 0) {
        p <- p + ggplot2::geom_text(
          data = dark_tiles,
          ggplot2::aes(x = .data$x, y = .data$y, label = sprintf("%.2f", .data$Correlation)),
          color = "white", fontface = "bold", size = 3.2
        )
      }
      if (nrow(light_tiles) > 0) {
        p <- p + ggplot2::geom_text(
          data = light_tiles,
          ggplot2::aes(x = .data$x, y = .data$y, label = sprintf("%.2f", .data$Correlation)),
          color = "black", fontface = "bold", size = 3.2
        )
      }
    }

    p <- p +
      ggplot2::geom_tile(
        data = tiles_df[tiles_df$is_diag, ],
        ggplot2::aes(x = .data$x, y = .data$y),
        fill = "grey93", color = "white", linewidth = 0.65
      ) +
      ggplot2::geom_text(
        data = tiles_df[tiles_df$is_diag, ],
        ggplot2::aes(x = .data$x, y = .data$y, label = .data$Var1),
        fontface = "bold", color = "black", size = 3.6
      )

    for (d_i in seq_len(nrow(diag_pts))) {
      p <- p +
        ggplot2::annotate("point", x = diag_pts$x[d_i] + 0.38, y = diag_pts$y[d_i] + 0.38, size = 3.0, color = "#2A9D8F", alpha = 0.85) +
        ggplot2::annotate("point", x = diag_pts$x[d_i] + 0.38, y = diag_pts$y[d_i] + 0.38, size = 1.2, color = "white")
    }

    if (nrow(arc_df) > 0) {
      if (line_color_by == "significance") {
        sig_colors <- c("p < 0.01" = "#D9534F", "p < 0.05" = "#4575B4", "p \u2265 0.05" = "grey75")
        for (k in seq_len(nrow(arc_df))) {
          row_k <- arc_df[k, ]
          p <- p + ggplot2::geom_curve(
            data = row_k,
            ggplot2::aes(
              x = .data$x_start, y = .data$y_start,
              xend = .data$x_end, yend = .data$y_end,
              color = .data$Significance,
              linewidth = .data$Mantel_R_Tier,
              linetype = .data$Line_Type
            ),
            curvature = row_k$curvature,
            alpha = 0.85
          )
        }
        p <- p + ggplot2::scale_color_manual(values = sig_colors, name = "Mantel p-value", drop = FALSE)
      } else {
        for (k in seq_len(nrow(arc_df))) {
          row_k <- arc_df[k, ]
          p <- p + ggplot2::geom_curve(
            data = row_k,
            ggplot2::aes(
              x = .data$x_start, y = .data$y_start,
              xend = .data$x_end, yend = .data$y_end,
              color = .data$Community,
              linewidth = .data$Mantel_R_Tier,
              linetype = .data$Line_Type
            ),
            curvature = row_k$curvature,
            alpha = 0.85
          )
        }
        p <- p + ggplot2::scale_color_brewer(palette = "Set1", name = "Community")
      }
    }

    for (m in seq_len(nrow(comm_nodes))) {
      p <- p +
        ggplot2::annotate("point", x = comm_nodes$comm_x[m], y = comm_nodes$comm_y[m], size = 6.5, color = comm_nodes$node_color[m]) +
        ggplot2::annotate("point", x = comm_nodes$comm_x[m], y = comm_nodes$comm_y[m], size = 2.5, color = "white") +
        ggplot2::annotate("text", x = comm_nodes$comm_x[m] + 0.32, y = comm_nodes$comm_y[m], label = comm_nodes$Community[m], fontface = "bold", size = 3.8, hjust = 0)
    }

    r_widths <- c("< 0.20" = 0.8, "0.20 - 0.40" = 1.8, "\u2265 0.40" = 3.2)

    p <- p +
      ggplot2::scale_fill_distiller(
        palette = color_palette,
        limit = c(-1, 1),
        direction = -1,
        name = paste0("Corr (", tools::toTitleCase(method), ")")
      ) +
      ggplot2::scale_linewidth_manual(
        values = r_widths,
        name = "Mantel's r",
        drop = FALSE
      ) +
      ggplot2::scale_linetype_manual(
        values = c("solid" = "solid", "dashed" = "dashed"),
        guide = "none"
      ) +
      ggplot2::coord_fixed(
        xlim = c(0.5, n_env + 3.2),
        ylim = c(0.4, n_env + 0.8),
        clip = "off"
      ) +
      ggplot2::theme_void(base_size = 11) +
      ggplot2::theme(
        legend.position = "right",
        legend.box = "vertical",
        legend.spacing.y = ggplot2::unit(4, "pt"),
        legend.title = ggplot2::element_text(face = "bold", size = 9),
        legend.text = ggplot2::element_text(size = 8.5),
        plot.title = ggplot2::element_text(face = "bold", size = 12, hjust = 0.35, margin = ggplot2::margin(b = 6)),
        plot.margin = ggplot2::margin(t = 10, r = 20, b = 10, l = 10)
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
    cor_long$TextColor <- ifelse(abs(cor_long$Correlation) > 0.65, "white", "black")

    mantel_results$Variable <- factor(mantel_results$Variable, levels = env_vars)
    mantel_results$Stars <- ifelse(mantel_results$p_value < 0.01, "**", ifelse(mantel_results$p_value < 0.05, "*", ""))
    mantel_results$Label <- sprintf("%.2f%s", mantel_results$Mantel_r, mantel_results$Stars)

    p_heat <- ggplot2::ggplot(cor_long, ggplot2::aes(x = .data$Var1, y = .data$Var2, fill = .data$Correlation)) +
      ggplot2::geom_tile(color = "white", linewidth = 0.5) +
      ggplot2::scale_fill_distiller(palette = color_palette, limit = c(-1, 1), direction = -1) +
      ggplot2::coord_fixed() +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::labs(title = "Environmental Correlation", x = NULL, y = NULL, fill = "r") +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, face = "bold", color = "black"),
        axis.text.y = ggplot2::element_text(face = "bold", color = "black"),
        panel.grid = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
      )

    if (show_cor_text) {
      p_heat <- p_heat +
        ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", .data$Correlation), color = .data$TextColor), fontface = "bold", size = 3.2) +
        ggplot2::scale_color_identity()
    }

    max_r <- max(c(0.2, max(mantel_results$Mantel_r, na.rm = TRUE) * 1.30))
    p_bar <- ggplot2::ggplot(mantel_results, ggplot2::aes(x = .data$Variable, y = .data$Mantel_r, fill = .data$Significance)) +
      ggplot2::geom_col(width = 0.6, color = "black", linewidth = 0.35) +
      ggplot2::geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
      ggplot2::geom_text(ggplot2::aes(label = .data$Label, hjust = ifelse(.data$Mantel_r >= 0, -0.2, 1.2)), fontface = "bold", size = 3.0) +
      ggplot2::coord_flip() +
      ggplot2::facet_wrap(~Community, ncol = 1) +
      ggplot2::scale_fill_manual(
        values = c("p < 0.01" = "#D9534F", "p < 0.05" = "#4575B4", "p \u2265 0.05" = "grey80"),
        drop = FALSE
      ) +
      ggplot2::scale_y_continuous(limits = c(min(c(0, min(mantel_results$Mantel_r) * 1.2)), max_r), expand = ggplot2::expansion(mult = c(0.05, 0.15))) +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::labs(title = "Mantel Linkages by Community", x = NULL, y = "Mantel r", fill = "Significance") +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(face = "bold"),
        axis.text.y = ggplot2::element_text(face = "bold"),
        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
        strip.text = ggplot2::element_text(face = "bold")
      )

    if (requireNamespace("patchwork", quietly = TRUE)) {
      final_plot <- patchwork::wrap_plots(p_heat, p_bar, ncol = 2, widths = c(1.4, 1.0))
    } else if (requireNamespace("cowplot", quietly = TRUE)) {
      final_plot <- cowplot::plot_grid(p_heat, p_bar, ncol = 2, rel_widths = c(1.4, 1.0))
    } else {
      message("Note: Install 'patchwork' or 'cowplot' for side-by-side composite bar layout. Returning heatmap.")
      final_plot <- p_heat
    }
  }

  return(list(
    Plot = final_plot,
    Correlation_Matrix = cor_res,
    Mantel_Results = mantel_results
  ))
}
