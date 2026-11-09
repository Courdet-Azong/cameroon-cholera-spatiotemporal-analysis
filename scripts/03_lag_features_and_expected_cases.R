# ============================================================
# 03_lag_features_and_expected_cases.R
# ============================================================
# Purpose:
# Prepare the cholera dataset for spatiotemporal modelling by
# creating lagged climate and flood variables, scaling selected
# covariates, generating temporal and regional indices, adjusting
# population estimates, and calculating expected cases.
#
# Note:
# The original surveillance and covariate datasets are not
# included in this public repository.
# ============================================================

# ------------------------------------------------------------
# 1. Load packages and prepare dataset
# ------------------------------------------------------------

library(readr)
library(dplyr)
library(tidyr)
library(lubridate)

cholera_inla <- read_csv("cholera_inla.csv")

cholera_inla <- cholera_inla %>%
  select(-E, -E2) %>%
  select(-matches("^(Rain_lag|Temp_lag|flood_lag)"))


# ------------------------------------------------------------
# 2. Create lagged climate variables
# ------------------------------------------------------------

cholera_inla <- cholera_inla %>%
  arrange(ID_HD, t_week) %>%
  group_by(ID_HD) %>%
  mutate(
    Rain_lag1 = lag(Rain_avg, 1),
    Rain_lag2 = lag(Rain_avg, 2),
    Rain_lag3 = lag(Rain_avg, 3),
    Temp_lag1 = lag(Temp, 1),
    Temp_lag2 = lag(Temp, 2),
    Temp_lag3 = lag(Temp, 3),
    Temp_lag4 = lag(Temp, 4)
  ) %>%
  ungroup()


# ------------------------------------------------------------
# 3. Scale lagged climate variables
# ------------------------------------------------------------

cholera_inla <- cholera_inla %>%
  mutate(
    Rain_lag1_sc = as.numeric(scale(Rain_lag1)),
    Rain_lag2_sc = as.numeric(scale(Rain_lag2)),
    Rain_lag3_sc = as.numeric(scale(Rain_lag3)),
    Temp_lag1_sc = as.numeric(scale(Temp_lag1)),
    Temp_lag2_sc = as.numeric(scale(Temp_lag2)),
    Temp_lag3_sc = as.numeric(scale(Temp_lag3)),
    Temp_lag4_sc = as.numeric(scale(Temp_lag4))
  )


# ------------------------------------------------------------
# 4. Create flood occurrence and lagged flood variables
# ------------------------------------------------------------

cholera_inla <- cholera_inla %>%
  mutate(
    flood_occurrence = ifelse(Rain_total > 200, 1, 0)
  ) %>%
  relocate(flood_occurrence, .after = flood_intensity_sc)

cholera_inla <- cholera_inla %>%
  arrange(ID_HD, t_week) %>%
  group_by(ID_HD) %>%
  mutate(
    flood_lag1 = lag(flood_occurrence, 1),
    flood_lag2 = lag(flood_occurrence, 2)
  ) %>%
  ungroup()

cholera_inla <- cholera_inla %>%
  mutate(
    flood_lag1_sc = as.numeric(scale(flood_lag1)),
    flood_lag2_sc = as.numeric(scale(flood_lag2))
  )


# ------------------------------------------------------------
# 5. Generate weekly and temporal indices
# ------------------------------------------------------------

week_seq <- tibble(
  week_start = seq(
    from = as.Date("2018-01-01"),
    to = as.Date("2023-12-31"),
    by = "week"
  )
) %>%
  mutate(
    Year = isoyear(week_start),
    Week = isoweek(week_start)
  ) %>%
  filter(Year >= 2018, Year <= 2023) %>%
  mutate(
    t_week = row_number(),
    t_year = Year - min(Year) + 1
  ) %>%
  select(Year, Week, t_week, t_year)

cholera_inla <- cholera_inla %>%
  select(-t_week) %>%
  left_join(week_seq, by = c("Year", "Week")) %>%
  mutate(
    week_id_reg = t_week,
    id_gamma = t_week,
    id_delta = t_week,
    yrid_alpha = t_year,
    yrid_U = t_year,
    id_U = ID_HD,
    rid_gamma = REG_ID,
    rid_delta = REG_ID,
    rid_Temp = REG_ID,
    rid_Prec = REG_ID,
    rid_Wpp = REG_ID,
    rid_Slope = REG_ID,
    rid_Wes = REG_ID,
    rid_Hosp = REG_ID,
    rid_pov = REG_ID,
    rid_sani = REG_ID,
    rid_flood = REG_ID,
    rid_conf = REG_ID
  )


# ------------------------------------------------------------
# 6. Prepare annual population estimates
# ------------------------------------------------------------

pop_df <- read_csv("population_data.csv")

pop_long <- pop_df %>%
  select(
    ID_HD,
    Pp_2023,
    `2018` = `2018Pop_sum`,
    `2019` = `2019Pop_sum`,
    `2020` = `2020Pop_sum`,
    `2021` = `2021Pop_sum`,
    `2022` = `2022Pop_sum`,
    `2023` = `2023Pop_sum`
  ) %>%
  pivot_longer(
    cols = c(`2018`, `2019`, `2020`, `2021`, `2022`, `2023`),
    names_to = "Year",
    values_to = "Pop"
  ) %>%
  mutate(
    Year = as.integer(Year)
  ) %>%
  group_by(ID_HD) %>%
  mutate(
    Pop_adj = Pp_2023 * (Pop / Pop[Year == 2023])
  ) %>%
  ungroup()

cholera_inla <- cholera_inla %>%
  left_join(pop_long, by = c("ID_HD", "Year"))


# ------------------------------------------------------------
# 7. Calculate district-level annual expected cases
# ------------------------------------------------------------

district_year <- cholera_inla %>%
  distinct(ID_HD, REG_ID, Year, Pop_adj) %>%
  left_join(
    cholera_inla %>%
      group_by(ID_HD, Year) %>%
      summarise(
        y = sum(y_casecount, na.rm = TRUE),
        .groups = "drop"
      ),
    by = c("ID_HD", "Year")
  )

annual_rate <- district_year %>%
  group_by(Year) %>%
  summarise(
    total_cases = sum(y, na.rm = TRUE),
    total_pop = sum(Pop_adj, na.rm = TRUE),
    rate = total_cases / total_pop,
    .groups = "drop"
  )

weeks_per_year <- cholera_inla %>%
  group_by(Year) %>%
  summarise(
    n_weeks_year = n_distinct(Week),
    .groups = "drop"
  )

district_year <- district_year %>%
  left_join(annual_rate, by = "Year") %>%
  left_join(weeks_per_year, by = "Year") %>%
  mutate(
    E_week = Pop_adj * rate / n_weeks_year,
    E_district = Pop_adj * rate
  )


# ------------------------------------------------------------
# 8. Join expected cases to the weekly dataset
# ------------------------------------------------------------

cholera_inla <- cholera_inla %>%
  left_join(
    district_year %>%
      select(ID_HD, Year, E_week, E_district),
    by = c("ID_HD", "Year")
  )


# ------------------------------------------------------------
# 9. Calculate regional expected cases
# ------------------------------------------------------------

region_year <- district_year %>%
  group_by(REG_ID, Year) %>%
  summarise(
    E_region = sum(E_district, na.rm = TRUE),
    .groups = "drop"
  )

cholera_inla <- cholera_inla %>%
  left_join(region_year, by = c("REG_ID", "Year"))


# ------------------------------------------------------------
# 10. Calculate Standardised Morbidity Ratio (SMR)
# ------------------------------------------------------------

cholera_inla <- cholera_inla %>%
  mutate(
    SMR = y_casecount / E_district
  )

summary(cholera_inla$SMR)


# ------------------------------------------------------------
# 11. Summarise observed and expected cases
# ------------------------------------------------------------

cat(
  "Observed cases:", sum(cholera_inla$y_casecount, na.rm = TRUE), "\n",
  "Weekly expected cases:", sum(cholera_inla$E_week, na.rm = TRUE), "\n"
)

district_year %>%
  group_by(Year) %>%
  summarise(
    observed = sum(y, na.rm = TRUE),
    expected = sum(E_district, na.rm = TRUE),
    .groups = "drop"
  )

cholera_inla %>%
  distinct(REG_ID, Year, E_region) %>%
  group_by(Year) %>%
  summarise(
    sum_region_E = sum(E_region, na.rm = TRUE),
    .groups = "drop"
  )


# ------------------------------------------------------------
# 12. Save prepared modelling dataset
# ------------------------------------------------------------

write_csv(
  cholera_inla,
  "cholera_inla_temporal_features.csv"
)
