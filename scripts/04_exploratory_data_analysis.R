# ============================================================
# 04_exploratory_data_analysis.R
# Exploratory data analysis of cholera model covariates
#
# Purpose:
#   - Compare covariate distributions by cholera occurrence
#   - Conduct Wilcoxon rank-sum tests
#   - Generate descriptive statistics
#   - Examine Spearman correlations
#   - Assess multicollinearity using VIF
#
# Input:
#   cholera_inla: integrated district-week dataset created by
#   the preceding data preparation and covariate integration scripts.
#
# Note:
#   The source dataset is not included in this public repository.
# ============================================================


# 1. Load packages ---------------------------------------------------------

library(tidyverse)
library(DT)
library(corrplot)
library(car)


# 2. Set output directory --------------------------------------------------

output_dir <- "outputs/eda"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}


# 3. Check input dataset ---------------------------------------------------

if (!exists("cholera_inla")) {
  stop(
    "The object 'cholera_inla' was not found. ",
    "Load the integrated dataset before running this script."
  )
}


# 4. Boxplots by cholera occurrence ---------------------------------------

box_data <- cholera_inla %>%
  select(
    y_cases1s0s,
    Rain_avg,
    Temp,
    flood_intensity,
    Pop_density,
    Poor_sanitation,
    Slope,
    building_density,
    poverty_index,
    water_exposure,
    water_pop_prop,
    conflict_binary
  )

box_long <- box_data %>%
  pivot_longer(
    cols = -y_cases1s0s,
    names_to = "Variable",
    values_to = "Value"
  )

boxplot_by_occurrence <- ggplot(
  box_long,
  aes(
    x = factor(y_cases1s0s),
    y = Value
  )
) +
  geom_boxplot(fill = "lightgrey") +
  facet_wrap(~ Variable, scales = "free") +
  labs(
    x = "Cholera occurrence (y)",
    y = "Value"
  ) +
  theme_bw()

print(boxplot_by_occurrence)

ggsave(
  filename = file.path(output_dir, "boxplots_by_cholera_occurrence.png"),
  plot = boxplot_by_occurrence,
  width = 12,
  height = 9,
  dpi = 300
)


# 5. Wilcoxon rank-sum tests ----------------------------------------------

# The outcome must contain exactly two groups for the rank-sum test.
if (dplyr::n_distinct(stats::na.omit(cholera_inla$y_cases1s0s)) != 2) {
  stop(
    "'y_cases1s0s' must contain exactly two non-missing groups ",
    "for the Wilcoxon rank-sum tests."
  )
}

covariates <- c(
  "Rain_avg",
  "Temp",
  "flood_intensity",
  "building_density",
  "Pop_density",
  "Poor_sanitation",
  "poverty_index",
  "Slope",
  "water_exposure",
  "water_pop_prop",
  "conflict_binary"
)

wilcox_results <- purrr::map_dfr(covariates, function(variable) {
  
  formula_test <- stats::reformulate(
    "y_cases1s0s",
    response = variable
  )
  
  test <- stats::wilcox.test(
    formula_test,
    data = cholera_inla,
    exact = FALSE
  )
  
  tibble(
    Variable = variable,
    W = unname(test$statistic),
    p_value = test$p.value
  )
})

print(wilcox_results)

write_csv(
  wilcox_results,
  file.path(output_dir, "wilcoxon_test_results.csv")
)


# 6. Descriptive statistics ------------------------------------------------

summary_variables <- c(
  "incidence_rate",
  "SMR",
  "Rain_avg",
  "flood_occurrence",
  "Temp",
  "Pop_adj",
  "Pop_density",
  "Poor_sanitation",
  "slope_mean",
  "building_density",
  "poverty_index",
  "water_prop",
  "pop_within1km_water",
  "TWI",
  "health_facilities_per1000",
  "conflict_events_weekly",
  "hospital_rate"
)

summary_stats <- cholera_inla %>%
  select(all_of(summary_variables)) %>%
  summarise(
    across(
      everything(),
      list(
        Min = ~ min(.x, na.rm = TRUE),
        Q1 = ~ quantile(.x, 0.25, na.rm = TRUE),
        Median = ~ median(.x, na.rm = TRUE),
        Mean = ~ mean(.x, na.rm = TRUE),
        SD = ~ sd(.x, na.rm = TRUE),
        Q3 = ~ quantile(.x, 0.75, na.rm = TRUE),
        Max = ~ max(.x, na.rm = TRUE)
      )
    )
  ) %>%
  pivot_longer(
    cols = everything(),
    names_to = c("Variable", ".value"),
    names_pattern = "^(.*)_(Min|Q1|Median|Mean|SD|Q3|Max)$"
  )

print(summary_stats)

DT::datatable(
  summary_stats,
  options = list(pageLength = 10),
  rownames = FALSE
)

write_csv(
  summary_stats,
  file.path(output_dir, "descriptive_statistics.csv")
)

# Additional descriptive check for elevation.
summary(cholera_inla$Elevation)


# 7. Spearman correlation analysis -----------------------------------------

correlation_variables <- c(
  "incidence_rate",
  "SMR",
  "Rain_avg",
  "flood_occurrence",
  "Temp",
  "Pop_adj",
  "Pop_density",
  "Poor_sanitation",
  "slope_mean",
  "building_density",
  "poverty_index",
  "water_prop",
  "pop_within1km_water",
  "TWI",
  "health_facilities_per1000",
  "conflict_events_weekly",
  "hospital_rate"
)

corr_df <- cholera_inla %>%
  select(all_of(correlation_variables))

cor_matrix <- cor(
  corr_df,
  method = "spearman",
  use = "complete.obs"
)

# Labels are kept in the same order as correlation_variables.
correlation_labels <- c(
  "Incidence rate",
  "SMR",
  "Precipitation",
  "Flood occurrence",
  "Temperature",
  "Population",
  "Population density",
  "Unimproved sanitation",
  "Slope",
  "Building density",
  "Poverty index",
  "Proportion water",
  "Population within 1 km of water",
  "TWI",
  "Healthcare facilities/1,000",
  "Conflict events",
  "Hospitalisation rate"
)

colnames(cor_matrix) <- correlation_labels
rownames(cor_matrix) <- correlation_labels

print(round(cor_matrix, 2))

png(
  filename = file.path(output_dir, "spearman_correlation_matrix.png"),
  width = 2400,
  height = 2000,
  res = 300
)

corrplot(
  cor_matrix,
  method = "color",
  type = "upper",
  addCoef.col = "black",
  tl.col = "black",
  tl.srt = 45,
  tl.cex = 0.7,
  number.cex = 0.6
)

dev.off()

write.csv(
  cor_matrix,
  file.path(output_dir, "spearman_correlation_matrix.csv"),
  row.names = TRUE
)


# 8. Variance inflation factor (VIF) ---------------------------------------

# Model 1: VIF assessment using the broader covariate set.
vif_model <- lm(
  y_casecount ~
    Rain_avg +
    Temp +
    flood_occurrence +
    Pop_adj +
    Pop_density +
    Poor_sanitation +
    slope_mean +
    Elevation +
    building_density +
    poverty_index +
    water_prop +
    pop_within1km_water +
    TWI +
    health_facilities_per1000 +
    conflict_events_weekly +
    hospital_rate,
  data = cholera_inla
)

vif_results_1 <- car::vif(vif_model)
print(vif_results_1)

write.csv(
  as.data.frame(vif_results_1),
  file.path(output_dir, "vif_model_1.csv"),
  row.names = TRUE
)


# Model 2: VIF assessment using the alternative covariate set.
vif_model2 <- lm(
  y_casecount ~
    Rain_avg +
    Temp +
    flood_intensity +
    Pop_density +
    Poor_sanitation +
    Slope +
    poverty_index +
    water_exposure +
    water_pop_prop +
    conflict_events,
  data = cholera_inla
)

vif_results_2 <- car::vif(vif_model2)
print(vif_results_2)

write.csv(
  as.data.frame(vif_results_2),
  file.path(output_dir, "vif_model_2.csv"),
  row.names = TRUE
)


# 9. End of analysis -------------------------------------------------------
# Outputs are saved in: outputs/eda/
