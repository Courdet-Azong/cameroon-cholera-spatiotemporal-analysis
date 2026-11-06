
# ============================================================
# PROJECT: Spatio-Temporal Modelling of Climate-driven Cholera Risk,
#          in a Data Uncertainty Environment
# PERIOD: 2018–2023
#
# SCRIPT: 01_data_restructuring_inla.R
#
# PURPOSE:
# Clean and restructure weekly cholera surveillance data,
# harmonise health district names to match the spatial boundaries,
# integrate spatial attributes and construct a complete district-week
# dataset for subsequent Bayesian modelling using INLA.
#
# OUTPUT:
# cholera_inla - Longitudinal district-level weekly dataset
#                structured for Bayesian spatiotemporal
#                modelling using INLA.
#
# NOTE:
# Original surveillance data and spatial boundary files
# are restricted and are not included in this repository.
# ============================================================


# ============================================================
# 1. LOAD REQUIRED PACKAGES
# ============================================================

library(dplyr)   #Data manipulation
library(readxl)  #Import Excel surveillance data
library(tidyr)   #Data restructuring and missing values
library(sf)      #Spatial data processing
library(terra)   #Spatial vector operations
library(lubridate) #Epidemiological week and date handling

# Import the cholera surveillance line list.
# The original dataset is not publicly distributed. 

cases <- read_excel("data/LL cholera_Copy.xlsx")

# ============================================================
# 3. CLEAN AND STANDARDISE SURVEILLANCE DATA
# ============================================================

# Extract epidemiological week and year from the
# original symptom onset date variable.

cases_clean <- cases %>%
  separate(
    `Epiweek date symptoms`,
    into = c("Week_symp", "Year"),
    sep = "_"
  ) %>%
  mutate(
    Week = as.integer(gsub("[Ss]", "", Week_symp)), #Extract week
    Year = as.integer(Year), #Extract year
    District_clean = toupper(
      trimws(`Health district of origin`)
    )
  ) %>%
  select(-Week_symp)

# ============================================================
# 4. HARMONISE HEALTH DISTRICT NAMES
# ============================================================

# Standardise district names to ensure consistency
# between surveillance records and spatial boundaries.

cases_clean <- cases_clean %>%
  mutate(
    District_clean = recode(District_clean,
                            
                            # punctuation / underscores
                            "NEW_BELL" = "NEW-BELL",
                            "CITE_DES_PALMIERS" = "CITÉ DES PALMIERS",
                            "NJOMBE_PENJA" = "NJOMBE PENJA",
                            
                            # hyphens
                            "BIYEM ASSI" = "BIYEM-ASSI",
                            "ELIG MFOMO" = "ELIG-MFOMO",
                            
                            # GAROUA roman numerals
                            "GAROUA I"  = "GAROUA 1",
                            "GAROUA II" = "GAROUA 2",
                            
                            # Shapefile spelling
                            "MAKARY" = "MAKARI",
                            
                            # Accents
                            "ZOETELE" = "ZOÉTELE",
                            "CITE DES PALMIERS" = "CITÉ DES PALMIERS", 
                            
                            # NIÉTÉ
                            "NIETE" = "NYE'ETÉ"
    )
  )

# ============================================================
# 5. FILTER STUDY PERIOD AND HANDLE MISSING DATA
# ============================================================

# Retain observations from the study period (2018–2023).

cases_clean <- cases_clean %>%
  filter(Year >= 2018, Year <= 2023)

# Check missing district identifiers before exclusion.

missing_districts <- sum(is.na(cases_clean$District_clean))

# Exclude records without a health district identifier.

cases_clean <- cases_clean %>%
  filter(!is.na(District_clean))

# ============================================================
# 6. AGGREGATE WEEKLY CHOLERA CASES
# ============================================================

# Aggregate individual surveillance records into
# weekly case counts by health district.

cholera_weekly <- cases_clean %>%
  group_by(Year, District_clean, Week) %>%
  summarise(
    Cases = n(),
    .groups = "drop"
  )

# ============================================================
# 7. CREATE WEEKLY TEMPORAL INDEX
# ============================================================

# Generate a continuous weekly sequence (time identifiers).
# This sequence covers 2017–2024, as in the original
# analytical workflow.

week_seq <- tibble(
  week_start = seq(
    from = as.Date("2017-01-01"),
    to = as.Date("2024-12-31"),
    by = "week"
  )
) %>%
  mutate(
    Year = isoyear(week_start),
    Week = isoweek(week_start),
    Month = month(week_start)
  ) %>%
  filter(Year >= 2017) %>%
  mutate(
    t_week = row_number(),
    t_month = (Year - min(Year)) * 12 + Month,
    t_year = Year - min(Year) + 1
  ) %>%
  select(Year, Week, Month, t_week, t_month, t_year)


# ============================================================
# 8. IMPORT AND PROCESS REGIONAL BOUNDARIES
# ============================================================

# Import regional administrative boundaries.
# These boundaries are retained for the subsequent
# regional spatial modelling workflow.

reg_shp <- st_read(
  "data/boundaries/Regions_2022.shp",
  quiet = TRUE
)

# Create regional identifier.

reg_shp <- reg_shp %>%
  mutate(
    REG_ID = as.integer(Code_Reg)
  )

# Convert regional boundaries to terra vector format.

reg_vect <- vect(reg_shp)

# Convert back to sf format.

reg_sf <- st_as_sf(reg_vect)

# ============================================================
# 9. IMPORT AND PROCESS HEALTH DISTRICT BOUNDARIES
# ============================================================

# Import health district boundaries with population data.

shape <- st_read(
  "data/boundaries/Health_District_Population.shp",
  quiet = TRUE
)

# Standardise district names and generate district identifiers.

shape_clean <- shape %>%
  mutate(
    District_clean = toupper(trimws(TITRE))
  ) %>%
  distinct(
    District_clean,
    .keep_all = TRUE
  ) %>%
  arrange(District_clean) %>%
  mutate(
    ID_HD = row_number(),
    REG_ID = as.integer(factor(CODE_REG))
  ) %>%
  select(
    ID_HD,
    ID_DS,
    CODE_DS,
    District_clean,
    NOM_REGION,
    CODE_RS,
    REG_ID,
    Pop_2023,
    geometry
  )

# Convert district boundaries to terra and sf formats.

districts_vect <- vect(shape_clean)

districts_sf <- st_as_sf(districts_vect)


# ============================================================
# 10. CONSTRUCT COMPLETE DISTRICT-WEEK GRID
# ============================================================

# Create epidemiological week.

district_time_grid <- tidyr::expand_grid(
  District_clean = shape_clean$District_clean,
  week_seq
)


# ============================================================
# 11. INTEGRATE WEEKLY CHOLERA CASE COUNTS
# ============================================================

# Join observed weekly case counts to the complete
# district-week grid.
#
# Missing case counts are treated as zero, representing
# district-weeks with no recorded cases.

cholera_full <- district_time_grid %>%
  left_join(
    cholera_weekly,
    by = c("District_clean", "Year", "Week")
  ) %>%
  mutate(
    y_casecount = replace_na(Cases, 0)
  )


# ============================================================
# 12. CONSTRUCT THE INLA-READY DATASET
# ============================================================

# Integrate district-level spatial attributes into
# the weekly cholera surveillance dataset.

cholera_inla <- cholera_full %>%
  left_join(
    st_drop_geometry(shape_clean),
    by = "District_clean"
  ) %>%
  relocate(
    ID_HD,
    ID_DS,
    CODE_DS,
    .before = District_clean
  ) %>%
  relocate(
    NOM_REGION,
    CODE_RS,
    REG_ID,
    Pop_2023,
    .before = Year
  ) %>%
  select(-Cases)


# ============================================================
# 13. FINAL DATASET STRUCTURE
# ============================================================

# Inspect the resulting dataset.

glimpse(cholera_inla)

# Check dimensions.

dim(cholera_inla)

# Check for missing values in key modelling variables.

colSums(
  is.na(
    cholera_inla %>%
      select(
        ID_HD,
        REG_ID,
        Year,
        Week,
        t_week,
        y_casecount
      )
  )
)

# ============================================================
# END OF SCRIPT
# ============================================================
