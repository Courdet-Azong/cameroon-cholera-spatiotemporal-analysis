# 12_hotspot_mapping.R
#
# Purpose:
#   Map district-level probabilities of elevated spatial risk
#   for each study year.
#
# Required objects:
#   - districts_sf
#   - hotspot_probabilities.csv
#
# Output:
#   Hotspot_Maps/

############################################################
# 1. Load required packages
############################################################

library(dplyr)
library(ggplot2)
library(sf)

############################################################
# 2. Load hotspot probabilities
############################################################

hotspot_df <- read.csv(
  "hotspot_probabilities.csv"
)

############################################################
# 3. Create output directory
############################################################

output_dir <- "Hotspot_Maps"

dir.create(
  output_dir,
  showWarnings = FALSE
)

############################################################
# 4. Produce annual hotspot maps
############################################################

years <- sort(unique(hotspot_df$Year))

for (y in years) {
  
  map_sf <- districts_sf %>%
    left_join(
      hotspot_df %>%
        filter(Year == y),
      by = "ID_HD"
    )
  
  p_hotspot <- ggplot(map_sf) +
    
    geom_sf(
      aes(fill = hotspot_probability),
      colour = "grey40",
      linewidth = 0.1
    ) +
    
    scale_fill_viridis_c(
      name = "Hotspot probability",
      limits = c(0, 1),
      na.value = "grey80"
    ) +
    
    labs(
      title = paste(
        "Probability of Elevated Spatial Risk",
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
        size = 12
      ),
      
      legend.text = element_text(
        size = 10
      )
    )
  
  ggsave(
    filename = file.path(
      output_dir,
      paste0("Hotspot_Probability_", y, ".png")
    ),
    plot = p_hotspot,
    width = 8,
    height = 6,
    dpi = 300
  )
}
