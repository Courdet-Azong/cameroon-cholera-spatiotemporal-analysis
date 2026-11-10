# 11_hotspot_analysis.R
#
# Purpose:
#   Estimate district-level probabilities of elevated spatial risk
#   using the structured and unstructured BYM spatial effects.
#
# Required objects:
#   - model_st
#   - cholera_inla
#
# Assumptions:
#   - 200 districts
#   - 400 spatial effects per year
#   - Six years of spatial effects
#
# Note:
#   The combined probability uses a normal approximation and
#   assumes independence between the structured and unstructured
#   posterior components. It is not an exact joint posterior
#   exceedance probability.
#
# Output:
#   hotspot_probabilities.csv

############################################################
# 1. Load package
############################################################

library(dplyr)

############################################################
# 2. Extract spatial effects
############################################################

bym_eff <- model_st$summary.random$id_U

years <- sort(unique(cholera_inla$Year))

bym_eff$Year <- rep(years, each = 400)

############################################################
# 3. Separate spatial components
############################################################

structured_eff <- bym_eff %>%
  filter(ID <= 200) %>%
  mutate(
    id_U = rep(1:200, times = 6)
  )

unstructured_eff <- bym_eff %>%
  filter(ID > 200) %>%
  mutate(
    id_U = rep(1:200, times = 6)
  )

############################################################
# 4. Combine spatial effects
############################################################

combined_eff <- structured_eff %>%
  select(
    id_U,
    Year,
    structured_mean = mean,
    structured_sd = sd
  ) %>%
  left_join(
    unstructured_eff %>%
      select(
        id_U,
        Year,
        unstructured_mean = mean,
        unstructured_sd = sd
      ),
    by = c("id_U", "Year")
  ) %>%
  mutate(
    combined_mean = structured_mean + unstructured_mean,
    combined_sd = sqrt(
      structured_sd^2 + unstructured_sd^2
    )
  )

############################################################
# 5. Calculate hotspot probabilities
############################################################

hotspot_df <- combined_eff %>%
  mutate(
    hotspot_probability = 1 - pnorm(
      0,
      mean = combined_mean,
      sd = combined_sd
    )
  )

############################################################
# 6. Attach district identifiers
############################################################

district_lookup <- cholera_inla %>%
  distinct(
    id_U,
    ID_HD
  )

hotspot_df <- hotspot_df %>%
  left_join(
    district_lookup,
    by = "id_U"
  ) %>%
  select(
    ID_HD,
    Year,
    structured_mean,
    unstructured_mean,
    combined_mean,
    combined_sd,
    hotspot_probability
  )

############################################################
# 7. Export hotspot probabilities
############################################################

write.csv(
  hotspot_df,
  "hotspot_probabilities.csv",
  row.names = FALSE
)
