usethis::create_package("D:/R/MyRpackage/aanova")


# 1. Generate the dummy data
set.seed(42)
hilsa_weight <- data.frame(
  Habitat = rep(c("Marine", "Estuary", "River"), each = 30),
  Weight_g = c(
    rnorm(30, mean = 980, sd = 110),  # Ocean/Marine weights
    rnorm(30, mean = 850, sd = 95),   # Estuary weights
    rnorm(30, mean = 720, sd = 80)    # River weights
  )
)

# Ensure the categorical variable is a factor
hilsa_weight$Habitat <- as.factor(hilsa_weight$Habitat)

# 2. Save it directly to the package as an .rda file
usethis::use_data(hilsa_weight, overwrite = TRUE)


usethis::use_r("data")



one_way_anova(
  data = hilsa_weight,
  factor_var = "Habitat",
  numeric_var = "Weight_g",
  factor_levels = c("River", "Estuary", "Marine"),
  y_limits = c(400, 1500),
  plot_type = "boxplot",
  sig_display = "stars"
)







# 1. Generate and save the dataset into the package
set.seed(42)
hilsa_two_way <- data.frame(
  Habitat = rep(c("Marine", "Estuary", "River"), each = 30),
  Season = rep(c("Monsoon", "Dry"), times = 45),
  Weight_g = c(
    rnorm(30, mean = 1050, sd = 90),
    rnorm(30, mean = 920, sd = 80),
    rnorm(30, mean = 780, sd = 70)
  )
)
hilsa_two_way$Weight_g <- hilsa_two_way$Weight_g + ifelse(hilsa_two_way$Season == "Monsoon", 60, -60)
hilsa_two_way$Habitat <- as.factor(hilsa_two_way$Habitat)
hilsa_two_way$Season <- as.factor(hilsa_two_way$Season)


library(ggplot2)
library(emmeans)
library(multcomp)
library(dplyr)
library(ggpubr)

source("R/two_way_anova.R")





usethis::use_data(hilsa_two_way, overwrite = TRUE)

# 2. Reload the package functions and data
devtools::load_all()

# Test the refactored function with the new hilsa_two_way dataset
res_two_way <- two_way_anova(
  data = hilsa_two_way,
  factor1_var = "Habitat",
  factor2_var = "Season",
  numeric_var = "Weight_g",
  factor1_levels = c("River", "Estuary", "Marine"),
  factor2_levels = c("Dry", "Monsoon"),
  y_limits = c(700, 1300),
  plot_type = "barplot"
)

print(res_two_way$Plot)
res_two_way$Plot


# three way
# 1. Load required libraries
library(ggplot2)
library(emmeans)
library(multcomp)
library(dplyr)

# 2. Generate the three-way dummy dataset
set.seed(42)
hilsa_three_way <- data.frame(
  Habitat = rep(c("Marine", "Estuary", "River"), each = 40),
  Season = rep(c("Monsoon", "Dry"), each = 20, times = 3),
  Size_Class = rep(c("Juvenile", "Adult"), times = 60),
  Weight_g = rnorm(120, mean = 850, sd = 80)
)

hilsa_three_way$Weight_g <- hilsa_three_way$Weight_g +
  ifelse(hilsa_three_way$Habitat == "Marine", 150, ifelse(hilsa_three_way$Habitat == "Estuary", 50, -100)) +
  ifelse(hilsa_three_way$Season == "Monsoon", 50, -50) +
  ifelse(hilsa_three_way$Size_Class == "Adult", 200, -200)

hilsa_three_way$Habitat <- as.factor(hilsa_three_way$Habitat)
hilsa_three_way$Season <- as.factor(hilsa_three_way$Season)
hilsa_three_way$Size_Class <- as.factor(hilsa_three_way$Size_Class)

# 3. Source the function directly from your R folder
source("R/three_way_anova.R")

# 4. Execute the function and print the plot
res_three_way <- three_way_anova(
  data = hilsa_three_way,
  factor1_var = "Habitat",
  factor2_var = "Season",
  factor3_var = "Size_Class",
  numeric_var = "Weight_g",
  factor1_levels = c("River", "Estuary", "Marine"),
  factor2_levels = c("Dry", "Monsoon"),
  factor3_levels = c("Juvenile", "Adult"),
  y_limits = c(0, 1600),
  plot_type = "barplot"
)

print(res_three_way$Plot)



set.seed(42)
hilsa_morphology <- data.frame(
  Habitat = rep(c("Marine", "Estuary", "River"), each = 30),
  Body_Depth_cm = c(rnorm(30, mean = 12.5, sd = 1.1), rnorm(30, mean = 10.2, sd = 0.9), rnorm(30, mean = 8.1, sd = 0.7)),
  Head_Length_cm = c(rnorm(30, mean = 15.2, sd = 1.3), rnorm(30, mean = 13.0, sd = 1.1), rnorm(30, mean = 10.8, sd = 0.8)),
  Fin_Length_cm = c(rnorm(30, mean = 9.4, sd = 0.8), rnorm(30, mean = 8.2, sd = 0.7), rnorm(30, mean = 6.9, sd = 0.6))
)

hilsa_morphology$Habitat <- as.factor(hilsa_morphology$Habitat)
usethis::use_data(hilsa_morphology, overwrite = TRUE)

library(ggplot2)
library(dplyr)
library(tidyr)
library(MASS)
library(RColorBrewer)

source("R/manova_analysis.R")

res_lda <- manova_analysis(
  data = hilsa_morphology,
  response_vars = c("Body_Depth_cm", "Head_Length_cm", "Fin_Length_cm"),
  factor_var = "Habitat",
  factor_levels = c("River", "Estuary", "Marine"),
  plot_type = "lda",
  color_palette = "Set1"
)

print(res_lda$MANOVA_Summary)
print(res_lda$Plot)


# ancova
set.seed(42)
hilsa_ancova <- data.frame(
  Habitat = rep(c("Marine", "Estuary", "River"), each = 30),
  Total_Length_cm = runif(90, min = 20, max = 45)
)

# Generate Weight as a function of Length and Habitat with an ANCOVA slope structure
hilsa_ancova$Weight_g <- 15 + 35 * hilsa_ancova$Total_Length_cm +
  ifelse(hilsa_ancova$Habitat == "Marine", 120,
         ifelse(hilsa_ancova$Habitat == "Estuary", 50, 0)) +
  rnorm(90, mean = 0, sd = 25)

hilsa_ancova$Habitat <- as.factor(hilsa_ancova$Habitat)
usethis::use_data(hilsa_ancova, overwrite = TRUE)



library(ggplot2)
library(dplyr)
library(emmeans)
library(RColorBrewer)

source("R/ancova_analysis.R")

res_ancova <- ancova_analysis(
  data = hilsa_ancova,
  response_var = "Weight_g",
  factor_var = "Habitat",
  covariate_var = "Total_Length_cm",
  factor_levels = c("River", "Estuary", "Marine"),
  color_palette = "Set1"
)

print(res_ancova$ANCOVA_Table)
print(res_ancova$Adjusted_Means)
print(res_ancova$Plot)


# regression

set.seed(42)
hilsa_regression <- data.frame(
  Total_Length_cm = runif(100, min = 15, max = 50)
)
# Simulate a non-linear power relationship (Length-Weight: W = 0.01 * L^3.05 + noise)
hilsa_regression$Weight_g <- 0.012 * (hilsa_regression$Total_Length_cm^3.04) * exp(rnorm(100, mean = 0, sd = 0.08))
hilsa_regression$Habitat <- sample(c("River", "Estuary", "Marine"), 100, replace = TRUE)

usethis::use_data(hilsa_regression, overwrite = TRUE)


library(ggplot2)
library(dplyr)
library(RColorBrewer)

source("R/regression_analysis.R")

res_reg <- regression_analysis(
  data = hilsa_regression,
  x_var = "Total_Length_cm",
  y_var = "Weight_g",
  fit_type = "polynomial",
  group_var = "Habitat",
  color_palette = "Set1"
)

print(res_reg$Plot)
