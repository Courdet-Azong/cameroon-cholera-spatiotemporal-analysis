# 09_spatial_effects_analysis.R
#
# Purpose:
#   Extract structured and unstructured spatial random effects
#   from the BYM model and calculate district-level relative risks.
#
# Model:
#   model_random
#
# Required objects from previous scripts:
#   - model_random
#   - cholera_inla
#
# Spatial structure:
#   200 districts
#   400 spatial random effects per year:
#     IDs 1-200: structured effects
#     IDs 201-400: unstructured effects
#
# Outputs:
#   - structured_spatial_effects.csv
#   - unstructured_spatial_effects.csv
#   - combined_spatial_relative_risk.csv

############################################################
# 1. Load required package
############################################################

library(dplyr)

############################################################
# 2. Extract BYM spatial random effects
############################################################

bym_eff <- model_random$summary.random$id_U

years <- sort(unique(cholera_inla$Year))

# Assign years to the spatial random effects.
# Each year contains 400 random effects.

bym_eff$Year <- rep(years, each = 400)

############################################################
# 3. Separate structured and unstructured effects
############################################################

structured_eff <- bym_eff %>%
  filter(ID <= 200)

unstructured_eff <- bym_eff %>%
  filter(ID > 200)

# Assign district identifiers to both components.

structured_eff$id_U <- rep(1:200, times = 6)

unstructured_eff$id_U <- rep(1:200, times = 6)

############################################################
# 4. Calculate structured spatial relative risks
############################################################

structured_eff <- structured_eff %>%
  mutate(
    RR = exp(mean)
  )

############################################################
# 5. Calculate combined spatial relative risks
############################################################

combined_eff <- structured_eff %>%
  select(
    id_U,
    Year,
    structured = mean
  ) %>%
  left_join(
    unstructured_eff %>%
      select(
        id_U,
        Year,
        unstructured = mean
      ),
    by = c("id_U", "Year")
  ) %>%
  mutate(
    RR_model = exp(structured + unstructured)
  )

############################################################
# 6. Attach district identifiers
############################################################

district_lookup <- cholera_inla %>%
  distinct(
    id_U,
    ID_HD
  )

# Structured spatial relative risks

structured_risk_df <- structured_eff %>%
  left_join(
    district_lookup,
    by = "id_U"
  ) %>%
  select(
    ID_HD,
    Year,
    RR_model = RR
  )

# Combined spatial relative risks

combined_risk_df <- combined_eff %>%
  left_join(
    district_lookup,
    by = "id_U"
  ) %>%
  select(
    ID_HD,
    Year,
    RR_model
  )

############################################################
# 7. Export spatial effects and relative risks
############################################################

write.csv(
  structured_eff,
  "structured_spatial_effects.csv",
  row.names = FALSE
)

write.csv(
  unstructured_eff,
  "unstructured_spatial_effects.csv",
  row.names = FALSE
)

write.csv(
  combined_risk_df,
  "combined_spatial_relative_risk.csv",
  row.names = FALSE
)
