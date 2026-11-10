# 08_temporal_effects_analysis.R
# Regional temporal effects from the spatiotemporal BYM model
#
# This script extracts the RW1 temporal effects from model_st,
# links them to the corresponding regions and dates, and produces
# regional temporal trend plots.
#
# Required objects from previous scripts:
#   - model_st
#   - cholera_inla
#
# The analysis uses 10 regions and 313 weekly time points.

##################################################
# Load packages
##################################################

library(dplyr)
library(ggplot2)

##################################################
# Extract temporal effects
##################################################

gamma_eff <- model_st$summary.random$id_gamma

# Assign regional and weekly identifiers
gamma_eff <- gamma_eff %>%
  mutate(
    rid_gamma = rep(1:10, each = 313),
    id_gamma = rep(1:313, times = 10)
  )

##################################################
# Link temporal effects to regional information
##################################################

gamma_df <- gamma_eff %>%
  left_join(
    cholera_inla %>%
      distinct(
        id_gamma,
        rid_gamma,
        NAM_REG
      ),
    by = c("id_gamma", "rid_gamma")
  )

##################################################
# Create weekly dates
##################################################

gamma_df <- gamma_df %>%
  mutate(
    Date = seq(
      from = as.Date("2018-01-01"),
      by = "week",
      length.out = 313
    )[id_gamma]
  )

# Exclude regions from the regional plots cause of zero cases
gamma_df <- gamma_df %>%
  filter(
    !NAM_REG %in% c("Adamawa", "North West")
  )

##################################################
# Calculate relative risks
##################################################

gamma_df <- gamma_df %>%
  mutate(
    RR = exp(mean),
    RR_low = exp(`0.025quant`),
    RR_high = exp(`0.975quant`)
  )

##################################################
# Create output directory
##################################################

dir.create(
  "Regional_Temporal_Effects",
  showWarnings = FALSE
)

##################################################
# Plot regional temporal effects
##################################################

for(reg in unique(gamma_df$NAM_REG)) {
  
  p <- ggplot(
    subset(gamma_df, NAM_REG == reg),
    aes(
      x = Date,
      y = mean
    )
  ) +
    
    geom_ribbon(
      aes(
        ymin = `0.025quant`,
        ymax = `0.975quant`
      ),
      fill = "grey80",
      alpha = 0.4
    ) +
    
    geom_line(
      linewidth = 1,
      colour = "black"
    ) +
    
    geom_hline(
      yintercept = 0,
      linetype = "dashed"
    ) +
    
    scale_x_date(
      breaks = seq(
        as.Date("2018-01-01"),
        as.Date("2023-01-01"),
        by = "1 year"
      ),
      date_labels = "%Y",
      limits = c(
        as.Date("2018-01-01"),
        as.Date("2023-12-31")
      )
    ) +
    
    labs(
      title = reg,
      x = "Year",
      y = "Posterior Mean"
    ) +
    
    theme_bw() +
    
    theme(
      plot.title = element_text(
        face = "bold",
        size = 18,
        hjust = 0.5
      ),
      axis.title.x = element_text(
        face = "bold",
        size = 14
      ),
      axis.title.y = element_text(
        face = "bold",
        size = 14
      ),
      axis.text.x = element_text(
        face = "bold",
        size = 12,
        angle = 45,
        hjust = 1
      ),
      axis.text.y = element_text(
        face = "bold",
        size = 12
      )
    )
  
  ggsave(
    filename = file.path(
      "Regional_Temporal_Effects",
      paste0(
        gsub(" ", "_", reg),
        "_Temporal_Effect.png"
      )
    ),
    plot = p,
    width = 8,
    height = 5,
    dpi = 300
  )
}
