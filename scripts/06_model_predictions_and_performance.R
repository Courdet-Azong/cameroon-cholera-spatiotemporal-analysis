# ============================================================
# 06_model_predictions_and_performance.R
# ============================================================
# Purpose:
# Extract fitted values from the final spatiotemporal BYM model,
# evaluate model performance, and visualise observed versus
# predicted cholera cases at national and regional levels.
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
# 2. Extract fitted values and 95% credible intervals
# ------------------------------------------------------------

cholera_inla <- cholera_inla %>%
  mutate(
    pred_mean = model_st$summary.fitted.values$mean,
    pred_lo   = model_st$summary.fitted.values$`0.025quant`,
    pred_hi   = model_st$summary.fitted.values$`0.975quant`
  )

summary(cholera_inla$pred_mean)


# ------------------------------------------------------------
# 3. Model fit statistics
# ------------------------------------------------------------

DIC <- round(model_st$dic$dic, 2)

WAIC <- round(model_st$waic$waic, 2)

LCPO <- round(
  -sum(log(model_st$cpo$cpo), na.rm = TRUE),
  2
)

cat(
  "DIC :", DIC, "\n",
  "WAIC:", WAIC, "\n",
  "LCPO:", LCPO, "\n"
)


# ------------------------------------------------------------
# 4. Calculate overall model performance metrics
# ------------------------------------------------------------

RMSE <- sqrt(
  mean(
    (cholera_inla$y_casecount -
       cholera_inla$pred_mean)^2,
    na.rm = TRUE
  )
)

MAE <- mean(
  abs(
    cholera_inla$y_casecount -
      cholera_inla$pred_mean
  ),
  na.rm = TRUE
)

correlation <- cor(
  cholera_inla$y_casecount,
  cholera_inla$pred_mean,
  use = "complete.obs"
)

R2 <- correlation^2

fit_stats <- data.frame(
  RMSE = RMSE,
  MAE = MAE,
  Correlation = correlation,
  R2 = R2
)

datatable(
  fit_stats,
  caption = "Overall Model Performance",
  rownames = FALSE
)


# ------------------------------------------------------------
# 5. Observed versus predicted cases by region and year
# ------------------------------------------------------------

dir.create(
  "Observed_vs_Predicted_Regional",
  showWarnings = FALSE
)

for (reg in sort(unique(cholera_inla$NAM_REG))) {
  
  for (y in sort(unique(cholera_inla$Year))) {
    
    weekly_compare <- cholera_inla %>%
      filter(
        NAM_REG == reg,
        Year == y
      ) %>%
      group_by(Week) %>%
      summarise(
        Observed  = sum(y_casecount, na.rm = TRUE),
        Predicted = sum(pred_mean, na.rm = TRUE),
        Lower95   = sum(pred_lo, na.rm = TRUE),
        Upper95   = sum(pred_hi, na.rm = TRUE),
        .groups = "drop"
      )
    
    p <- ggplot(
      weekly_compare,
      aes(x = Week)
    ) +
      
      geom_ribbon(
        aes(
          ymin = Lower95,
          ymax = Upper95
        ),
        alpha = 0.25
      ) +
      
      geom_line(
        aes(
          y = Observed,
          colour = "Observed"
        ),
        linewidth = 1.2
      ) +
      
      geom_line(
        aes(
          y = Predicted,
          colour = "Predicted"
        ),
        linewidth = 1.2
      ) +
      
      scale_colour_manual(
        values = c(
          "Observed" = "black",
          "Predicted" = "red"
        )
      ) +
      
      labs(
        title = paste(
          "Observed vs Predicted Cholera Cases",
          reg,
          y
        ),
        x = "Week",
        y = "Cases",
        colour = ""
      ) +
      
      theme_minimal()
    
    ggsave(
      filename = file.path(
        "Observed_vs_Predicted_Regional",
        paste0(
          gsub("[^A-Za-z0-9]", "_", reg),
          "_",
          y,
          ".png"
        )
      ),
      plot = p,
      width = 8,
      height = 5,
      dpi = 300
    )
  }
}


# ------------------------------------------------------------
# 6. Observed versus predicted scatterplots by year
# ------------------------------------------------------------

dir.create(
  "Scatterplots_Year",
  showWarnings = FALSE
)

years <- sort(unique(cholera_inla$Year))

metrics_table <- data.frame(
  Year = years,
  Pearson_r = NA,
  Spearman_rho = NA,
  R2 = NA
)

for (i in seq_along(years)) {
  
  y <- years[i]
  
  plot_df <- cholera_inla %>%
    filter(Year == y)
  
  pearson_r <- suppressWarnings(
    cor(
      plot_df$y_casecount,
      plot_df$pred_mean,
      use = "complete.obs",
      method = "pearson"
    )
  )
  
  spearman_r <- suppressWarnings(
    cor(
      plot_df$y_casecount,
      plot_df$pred_mean,
      use = "complete.obs",
      method = "spearman"
    )
  )
  
  r2 <- pearson_r^2
  
  metrics_table$Pearson_r[i] <- pearson_r
  metrics_table$Spearman_rho[i] <- spearman_r
  metrics_table$R2[i] <- r2
  
  p <- ggplot(
    plot_df,
    aes(
      x = pred_mean,
      y = y_casecount
    )
  ) +
    
    geom_point(alpha = 0.5) +
    
    geom_abline(
      slope = 1,
      intercept = 0,
      linetype = "dashed"
    ) +
    
    labs(
      title = paste(
        "Observed vs Predicted Cases",
        y
      ),
      subtitle = paste0(
        "Pearson r = ", round(pearson_r, 3),
        " | Spearman rho = ", round(spearman_r, 3),
        " | R² = ", round(r2, 3)
      ),
      x = "Predicted Cases",
      y = "Observed Cases"
    ) +
    
    theme_bw()
  
  ggsave(
    filename = file.path(
      "Scatterplots_Year",
      paste0("Scatter_", y, ".png")
    ),
    plot = p,
    width = 7,
    height = 5,
    dpi = 300
  )
}

write.csv(
  metrics_table,
  file.path(
    "Scatterplots_Year",
    "Model_Performance_By_Year.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 7. Observed versus predicted scatterplots by region and year
# ------------------------------------------------------------

dir.create(
  "Scatterplots_Region_Year",
  showWarnings = FALSE
)

regions <- sort(unique(cholera_inla$NAM_REG))

metrics_region_year <- expand.grid(
  Year = years,
  NAM_REG = regions
) %>%
  arrange(Year, NAM_REG)

metrics_region_year$Pearson_r <- NA
metrics_region_year$Spearman_rho <- NA
metrics_region_year$R2 <- NA

for (i in seq_len(nrow(metrics_region_year))) {
  
  y <- metrics_region_year$Year[i]
  reg <- metrics_region_year$NAM_REG[i]
  
  plot_df <- cholera_inla %>%
    filter(
      Year == y,
      NAM_REG == reg
    )
  
  if (nrow(plot_df) > 1) {
    
    # Correlations
    
    pearson_r <- suppressWarnings(
      cor(
        plot_df$y_casecount,
        plot_df$pred_mean,
        use = "complete.obs",
        method = "pearson"
      )
    )
    
    spearman_r <- suppressWarnings(
      cor(
        plot_df$y_casecount,
        plot_df$pred_mean,
        use = "complete.obs",
        method = "spearman"
      )
    )
    
    r2 <- pearson_r^2
    
    # Store performance metrics
    
    metrics_region_year$Pearson_r[i] <- pearson_r
    metrics_region_year$Spearman_rho[i] <- spearman_r
    metrics_region_year$R2[i] <- r2
    
    # Generate scatterplot
    
    p <- ggplot(
      plot_df,
      aes(
        x = pred_mean,
        y = y_casecount
      )
    ) +
      
      geom_point(alpha = 0.5) +
      
      geom_abline(
        slope = 1,
        intercept = 0,
        linetype = "dashed"
      ) +
      
      labs(
        title = paste(
          "Observed vs Predicted:",
          reg,
          y
        ),
        subtitle = paste0(
          "Pearson r = ",
          round(pearson_r, 3),
          " | Spearman rho = ",
          round(spearman_r, 3),
          " | R² = ",
          round(r2, 3)
        ),
        x = "Predicted Cases",
        y = "Observed Cases"
      ) +
      
      theme_bw()
    
    ggsave(
      filename = file.path(
        "Scatterplots_Region_Year",
        paste0(
          gsub("[/ ]", "_", reg),
          "_",
          y,
          ".png"
        )
      ),
      plot = p,
      width = 7,
      height = 5,
      dpi = 300
    )
  }
}

write.csv(
  metrics_region_year,
  file.path(
    "Scatterplots_Region_Year",
    "Model_Performance_By_Region_Year.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 8. Display regional and annual performance summaries
# ------------------------------------------------------------

datatable(
  metrics_table,
  caption = "Model Performance by Year",
  rownames = FALSE
)

datatable(
  metrics_region_year,
  caption = "Model Performance by Region and Year",
  rownames = FALSE
)


# ============================================================
# END OF SCRIPT
# ============================================================
