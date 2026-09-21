# R/globals.R

utils::globalVariables(c(
  # General & Heatmaps
  "Value", "Variable", "Species", "Comparison", "Contribution",
  "X", "Y", "x", "y", "cos2", "Correlation", "Var1", "Var2",
  # ANOVA & Group Summaries
  "Group", "mean_val", "sd_val", "se_val", "error_margin",
  "upper_limit", ".group",
  # GLM & Regression
  "Term", "Estimate", "Effect_Size", "Lower_CI", "Upper_CI",
  # MANOVA & LDA
  "Measurement", "LD1", "LD2",
  # Mantel & Network
  "Mantel_r", "Significance", "Env_X", "Env_Y", "Comm_X", "Comm_Y",
  "p_cat", "r_cat", "Community",
  # Ordination (NMDS, PCA, RDA)
  "NMDS1", "NMDS2", "PC1", "PC2", "RDA1", "RDA2",
  # Accumulation & Misc
  "Sites", "Richness", "Lower", "Upper", ":="
))
