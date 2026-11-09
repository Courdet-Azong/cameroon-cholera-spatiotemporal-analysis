# ============================================================
# 07_fixed_and_random_effects.R
# ============================================================
# Purpose:
# Extract and summarise fixed effects, annual effects,
# hyperparameters and regional random effects from the final
# spatiotemporal BYM model.
#
# Model:
# model_st (BYM + temporal effects + covariates)
#
# Prerequisites:
# - cholera_inla: prepared modelling dataset
# - model_st: fitted INLA model
#
# Note:
# The underlying surveillance data are restricted and are not
# included in this public repository.
# ============================================================


# ------------------------------------------------------------
# 1. Load packages
# ------------------------------------------------------------

library(dplyr)
library(ggplot2)
library(DT)


# ------------------------------------------------------------
# 2. Extract fixed effects
# ------------------------------------------------------------

fixed <- model_st$summary.fixed

fixed_df <- data.frame(
  Covariate  = rownames(fixed),
  Mean       = round(fixed$mean, 4),
  SD         = round(fixed$sd, 4),
  Lower95    = round(fixed$`0.025quant`, 4),
  Upper95    = round(fixed$`0.975quant`, 4),
  RR         = round(exp(fixed$mean), 4),
  RR_Lower95 = round(exp(fixed$`0.025quant`), 4),
  RR_Upper95 = round(exp(fixed$`0.975quant`), 4)
)

rownames(fixed_df) <- NULL

datatable(
  fixed_df,
  caption = "Fixed Effects and Relative Risks",
  rownames = FALSE,
  options = list(
    pageLength = 20,
    scrollX = TRUE
  )
)


# ------------------------------------------------------------
# 3. Extract covariate effects
# ------------------------------------------------------------

covariates_df <- fixed_df %>%
  filter(Covariate != "(Intercept)")

datatable(
  covariates_df,
  caption = "Fixed Effect Coefficients and Relative Risks",
  rownames = FALSE,
  options = list(
    pageLength = 40,
    scrollX = TRUE
  )
)


# ------------------------------------------------------------
# 4. Extract annual effects
# ------------------------------------------------------------

year_effects <- model_st$summary.random$yrid_alpha

datatable(
  year_effects,
  caption = "Annual Random Effects"
)

national_rr <- data.frame(
  Year = 2018:2023,
  Mean = year_effects$mean,
  Lower_CrI = year_effects$`0.025quant`,
  Upper_CrI = year_effects$`0.975quant`,
  RR = exp(year_effects$mean),
  RR_Lower = exp(year_effects$`0.025quant`),
  RR_Upper = exp(year_effects$`0.975quant`)
)

datatable(
  national_rr,
  caption = "Annual Relative Risks",
  rownames = FALSE
)


# ------------------------------------------------------------
# 5. Visualise annual effects
# ------------------------------------------------------------

ggplot(
  national_rr,
  aes(x = Year, y = Mean)
) +
  
  geom_line(linewidth = 1) +
  
  geom_point(size = 3) +
  
  geom_line(
    aes(y = Lower_CrI),
    linetype = "dashed",
    linewidth = 0.8
  ) +
  
  geom_line(
    aes(y = Upper_CrI),
    linetype = "dashed",
    linewidth = 0.8
  ) +
  
  geom_hline(
    yintercept = 0,
    linetype = "dotted"
  ) +
  
  labs(
    x = "Year",
    y = "Posterior mean",
    title = "National Annual Cholera Transmission",
    subtitle = "Posterior mean and 95% credible intervals"
  ) +
  
  theme_bw()


# ------------------------------------------------------------
# 6. Extract model hyperparameters
# ------------------------------------------------------------

hyper_df <- model_st$summary.hyperpar

datatable(
  hyper_df,
  caption = "Model Hyperparameters",
  rownames = TRUE
)


# ------------------------------------------------------------
# 7. Calculate spatial variance components
# ------------------------------------------------------------

tau_u <- model_st$summary.hyperpar[
  "Precision for id_U (spatial component)",
  "mean"
]

tau_v <- model_st$summary.hyperpar[
  "Precision for id_U (iid component)",
  "mean"
]

spatial_variance <- data.frame(
  Structured_Variance = 1 / tau_u,
  Unstructured_Variance = 1 / tau_v,
  Spatial_Fraction =
    (1 / tau_u) / ((1 / tau_u) + (1 / tau_v))
)

spatial_variance


# ------------------------------------------------------------
# 8. Extract regional rainfall effects
# ------------------------------------------------------------

rain_eff <- model_st$summary.random$Rain_rw1

rain_eff$rid_Prec <- rep(1:10, each = 20)

region_lookup <- cholera_inla %>%
  distinct(rid_Prec, NAM_REG) %>%
  arrange(rid_Prec)

rain_eff <- rain_eff %>%
  left_join(
    region_lookup,
    by = "rid_Prec"
  )

rain_eff$RR <- exp(rain_eff$mean)
rain_eff$RR_low <- exp(rain_eff$`0.025quant`)
rain_eff$RR_high <- exp(rain_eff$`0.975quant`)


# ------------------------------------------------------------
# 9. Plot regional rainfall effects
# ------------------------------------------------------------

dir.create(
  "Rain_Nonlinear_mean",
  showWarnings = FALSE
)

for (reg in unique(rain_eff$NAM_REG)) {
  
  p <- ggplot(
    subset(rain_eff, NAM_REG == reg),
    aes(x = ID, y = mean)
  ) +
    
    geom_line(linewidth = 1) +
    
    geom_line(
      aes(y = `0.025quant`),
      linetype = "dashed"
    ) +
    
    geom_line(
      aes(y = `0.975quant`),
      linetype = "dashed"
    ) +
    
    geom_hline(
      yintercept = 0,
      linetype = "dotted"
    ) +
    
    labs(
      title = reg,
      x = "Lag 2 Precipitation (mm)",
      y = "Posterior means"
    ) +
    
    theme_bw() +
    
    theme(
      plot.title = element_text(
        face = "bold",
        hjust = 0.5
      ),
      axis.title.x = element_text(face = "bold"),
      axis.title.y = element_text(face = "bold"),
      axis.text.x = element_text(face = "bold"),
      axis.text.y = element_text(face = "bold")
    )
  
  ggsave(
    filename = file.path(
      "Rain_Nonlinear_mean",
      paste0("Rain_", reg, ".png")
    ),
    plot = p,
    width = 7,
    height = 5,
    dpi = 300
  )
}


# ------------------------------------------------------------
# 10. Extract regional temperature effects
# ------------------------------------------------------------

temp_eff <- model_st$summary.random$Temp_rw1

temp_eff$rid_Temp <- rep(1:10, each = 20)

region_lookup <- cholera_inla %>%
  distinct(rid_Temp, NAM_REG) %>%
  arrange(rid_Temp)

temp_eff <- temp_eff %>%
  left_join(
    region_lookup,
    by = "rid_Temp"
  )

temp_eff$RR <- exp(temp_eff$mean)
temp_eff$RR_low <- exp(temp_eff$`0.025quant`)
temp_eff$RR_high <- exp(temp_eff$`0.975quant`)


# ------------------------------------------------------------
# 11. Plot regional temperature effects
# ------------------------------------------------------------

dir.create(
  "Temperature_Nonlinear_mean",
  showWarnings = FALSE
)

for (reg in unique(temp_eff$NAM_REG)) {
  
  p <- ggplot(
    subset(temp_eff, NAM_REG == reg),
    aes(x = ID, y = mean)
  ) +
    
    geom_line(linewidth = 1) +
    
    geom_line(
      aes(y = `0.025quant`),
      linetype = "dashed"
    ) +
    
    geom_line(
      aes(y = `0.975quant`),
      linetype = "dashed"
    ) +
    
    geom_hline(
      yintercept = 0,
      linetype = "dotted"
    ) +
    
    labs(
      title = reg,
      x = "Lag 2 Temperature (°C)",
      y = "Posterior means"
    ) +
    
    theme_bw() +
    
    theme(
      plot.title = element_text(
        face = "bold",
        hjust = 0.5
      ),
      axis.title.x = element_text(face = "bold"),
      axis.title.y = element_text(face = "bold"),
      axis.text.x = element_text(face = "bold"),
      axis.text.y = element_text(face = "bold")
    )
  
  ggsave(
    filename = file.path(
      "Temperature_Nonlinear_mean",
      paste0("Temperature_", reg, ".png")
    ),
    plot = p,
    width = 7,
    height = 5,
    dpi = 300
  )
}


# ------------------------------------------------------------
# 12. Extract regional flood effects
# ------------------------------------------------------------

flood_eff <- model_st$summary.random$flood
