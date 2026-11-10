# 13_spatial_effect_visualisation.R
#
# Purpose:
#   Visualise the distributions of structured and unstructured
#   spatial random effects across the study years.
#
# Required objects:
#   - model_st
#   - cholera_inla
#
# Outputs:
#   Spatial_Effect_Visualisations/
#     Structured_Spatial_Effects.png
#     Unstructured_Spatial_Effects.png

############################################################
# 1. Load required packages
############################################################

library(dplyr)
library(ggplot2)

############################################################
# 2. Extract spatial random effects
############################################################

bym_eff <- model_st$summary.random$id_U

years <- sort(unique(cholera_inla$Year))

bym_eff$Year <- rep(years, each = 400)

############################################################
# 3. Separate structured and unstructured effects
############################################################

structured_eff <- bym_eff %>%
  filter(ID <= 200) %>%
  mutate(
    id_U = rep(1:200, times = 6),
    Effect = "Structured"
  )

unstructured_eff <- bym_eff %>%
  filter(ID > 200) %>%
  mutate(
    id_U = rep(1:200, times = 6),
    Effect = "Unstructured"
  )

############################################################
# 4. Create output directory
############################################################

output_dir <- "Spatial_Effect_Visualisations"

dir.create(
  output_dir,
  showWarnings = FALSE
)

############################################################
# 5. Visualise structured spatial effects
############################################################

p_structured <- ggplot(
  structured_eff,
  aes(x = mean)
) +
  
  geom_histogram(
    bins = 30,
    colour = "black"
  ) +
  
  geom_vline(
    xintercept = 0,
    linetype = "dashed"
  ) +
  
  facet_wrap(~ Year) +
  
  labs(
    title = "Distribution of Structured Spatial Effects",
    x = "Posterior mean spatial effect",
    y = "Number of districts"
  ) +
  
  theme_bw() +
  
  theme(
    plot.title = element_text(
      face = "bold",
      size = 14,
      hjust = 0.5
    )
  )

ggsave(
  filename = file.path(
    output_dir,
    "Structured_Spatial_Effects.png"
  ),
  plot = p_structured,
  width = 10,
  height = 7,
  dpi = 300
)

############################################################
# 6. Visualise unstructured spatial effects
############################################################

p_unstructured <- ggplot(
  unstructured_eff,
  aes(x = mean)
) +
  
  geom_histogram(
    bins = 30,
    colour = "black"
  ) +
  
  geom_vline(
    xintercept = 0,
    linetype = "dashed"
  ) +
  
  facet_wrap(~ Year) +
  
  labs(
    title = "Distribution of Unstructured Spatial Effects",
    x = "Posterior mean spatial effect",
    y = "Number of districts"
  ) +
  
  theme_bw() +
  
  theme(
    plot.title = element_text(
      face = "bold",
      size = 14,
      hjust = 0.5
    )
  )

ggsave(
  filename = file.path(
    output_dir,
    "Unstructured_Spatial_Effects.png"
  ),
  plot = p_unstructured,
  width = 10,
  height = 7,
  dpi = 300
)
