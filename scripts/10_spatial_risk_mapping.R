# 10_spatial_risk_mapping.R
#
# Purpose:
#   Visualise annual district-level relative risks derived from
#   the structured and combined spatial effects of the BYM model.
#
# Required objects from previous scripts:
#   - districts_sf
#   - structured_risk_df
#   - combined_risk_df
#
# Required packages:
#   - dplyr
#   - ggplot2
#   - sf
#   - ggspatial
#
# Outputs:
#   - bym_RR_Maps/
#   - Combined_RR_Maps/

############################################################
# 1. Load required packages
############################################################

library(dplyr)
library(ggplot2)
library(sf)
library(ggspatial)

############################################################
# 2. Define relative-risk categories
############################################################

rr_breaks <- c(0, 0.9, 1.0, 2.0, 5.0, Inf)

rr_labels <- c(
  "[0,0.9]",
  "(0.9,1.0]",
  "(1.0,2.0]",
  "(2.0,5.0]",
  "(5.0,+]"
)

rr_colours <- c(
  "[0,0.9]"   = "#fff5f0",
  "(0.9,1.0]" = "#fcbba1",
  "(1.0,2.0]" = "#fb6a4a",
  "(2.0,5.0]" = "#de2d26",
  "(5.0,+]"   = "#7f0000"
)

############################################################
# 3. Create output directories
############################################################

dir.create(
  "bym_RR_Maps",
  showWarnings = FALSE
)

dir.create(
  "Combined_RR_Maps",
  showWarnings = FALSE
)

############################################################
# 4. Map structured spatial relative risks
############################################################

years <- sort(unique(structured_risk_df$Year))

for (y in years) {
  
  map_sf <- districts_sf %>%
    left_join(
      structured_risk_df %>%
        filter(Year == y),
      by = "ID_HD"
    ) %>%
    mutate(
      RR_cat = cut(
        RR_model,
        breaks = rr_breaks,
        labels = rr_labels,
        include.lowest = TRUE
      )
    )
  
  rr_counts <- table(
    factor(
      map_sf$RR_cat,
      levels = rr_labels
    )
  )
  
  legend_labels <- paste0(
    names(rr_counts),
    ": ",
    as.numeric(rr_counts)
  )
  
  p_rr <- ggplot(map_sf) +
    
    geom_sf(
      aes(fill = RR_cat),
      colour = "grey40",
      linewidth = 0.1
    ) +
    
    scale_fill_manual(
      values = rr_colours,
      labels = legend_labels,
      name = "Relative Risk",
      na.value = "grey80",
      drop = FALSE
    ) +
    
    annotation_scale(
      location = "bl",
      width_hint = 0.3
    ) +
    
    annotation_north_arrow(
      location = "tl",
      which_north = "true",
      style = north_arrow_fancy_orienteering
    ) +
    
    labs(
      title = paste(
        "Relative Risk Associated with Spatial Random Effects",
        y
      )
    ) +
    
    theme_minimal() +
    
    theme(
      axis.text = element_blank(),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      
      plot.title = element_text(
        face = "bold",
        size = 16,
        hjust = 0.5
      ),
      
      legend.title = element_text(
        face = "bold",
        size = 14
      ),
      
      legend.text = element_text(
        face = "bold",
        size = 12
      )
    )
  
  ggsave(
    filename = file.path(
      "bym_RR_Maps",
      paste0("bym_RR_", y, ".png")
    ),
    plot = p_rr,
    width = 8,
    height = 6,
    dpi = 300
  )
}

############################################################
# 5. Map combined spatial relative risks
############################################################

years <- sort(unique(combined_risk_df$Year))

for (y in years) {
  
  map_sf <- districts_sf %>%
    left_join(
      combined_risk_df %>%
        filter(Year == y),
      by = "ID_HD"
    ) %>%
    mutate(
      RR_cat = cut(
        RR_model,
        breaks = rr_breaks,
        labels = rr_labels,
        include.lowest = TRUE
      )
    )
  
  rr_counts <- table(
    factor(
      map_sf$RR_cat,
      levels = rr_labels
    )
  )
  
  legend_labels <- paste0(
    names(rr_counts),
    ": ",
    as.numeric(rr_counts)
  )
  
  p_rr <- ggplot(map_sf) +
    
    geom_sf(
      aes(fill = RR_cat),
      colour = "grey40",
      linewidth = 0.1
    ) +
    
    scale_fill_manual(
      values = rr_colours,
      labels = legend_labels,
      name = "Relative Risk",
      na.value = "grey80",
      drop = FALSE
    ) +
    
    labs(
      title = paste(
        "Combined Relative Risk",
        y
      )
    ) +
    
    theme_minimal() +
    
    theme(
      axis.text = element_blank(),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      
      plot.title = element_text(
        face = "bold",
        size = 14
      ),
      
      legend.title = element_text(
        face = "bold"
      )
    )
  
  ggsave(
    filename = file.path(
      "Combined_RR_Maps",
      paste0("Combined_RR_", y, ".png")
    ),
    plot = p_rr,
    width = 8,
    height = 6,
    dpi = 300
  )
}
