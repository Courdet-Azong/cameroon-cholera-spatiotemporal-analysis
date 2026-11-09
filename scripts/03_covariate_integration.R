# ============================================================
# PROJECT: Spatiotemporal Modelling of Cholera Risk in Cameroon
# SCRIPT: 02_covariate_integration.R
# PURPOSE:
# Integrate prepared environmental, hydroclimatic,
# socioeconomic and health-related variables into cholera_inla.
#
# INPUT:
# cholera_inla - District-level weekly surveillance dataset
# Prepared district-level and weekly covariate CSV files
#
# OUTPUT:
# cholera_inla - Longitudinal district-level weekly dataset
#                with integrated study variables for Bayesian
#                spatiotemporal modelling using INLA.
# ============================================================


# 1. Load packages --------------------------------------------

library(dplyr)
library(readr)
library(tidyr)
library(sf)
library(purrr)
library(lubridate)
library(terra)
library(whitebox)


# 2. Set file paths -------------------------------------------

# Main modelling dataset (already prepared)
main_data_file <- "data/processed/cholera_inla.csv"

# Health district shapefile
district_file <- "data/spatial/health_districts.shp"

# CHIRPS daily rainfall NetCDF files
rain_dir <- "data/raw/CHIRPS"

# ERA5 daily temperature NetCDF files
temp_dir <- "data/raw/ERA5"

# SRTM elevation raster
dem_file <- "data/raw/Elevation/Elevation.tif"

# Output directory
output_dir <- "data/processed"


# 3. Load main dataset and health district boundaries ----------

cholera_inla <- read_csv(main_data_file, show_col_types = FALSE)

districts_sf <- st_read(district_file, quiet = TRUE)

# Convert district polygons to terra vector
districts_vect <- vect(districts_sf)


# 4. Calculate incidence rate and SMR -------------------------

cholera_inla <- cholera_inla %>%
  mutate(
    incidence_rate = (y_casecount / Pop_adj) * 100000,
    SMR = y_casecount / E_district
  ) %>%
  relocate(incidence_rate, .after = y_casecount) %>%
  relocate(SMR, .after = E_district)


# ============================================================
# 5. ENVIRONMENTAL AND HYDROCLIMATIC VARIABLES
# ============================================================


# 5.1 Integrate rainfall variables derived from CHIRPS --------

# Daily rainfall is extracted for each health district.
# Daily values are then aggregated into weekly totals
# and weekly averages.

precip_files <- list.files(
  rain_dir,
  pattern = "\\.nc$",
  full.names = TRUE
)

if (length(precip_files) == 0) {
  stop("No CHIRPS NetCDF files found. Check rain_dir.")
}


# Function to extract weekly rainfall from one NetCDF file

get_weekly_rain <- function(nc_file, districts_vect) {
  
  message("Processing rainfall: ", basename(nc_file))
  
  r <- rast(nc_file)
  
  # Reproject district boundaries to raster CRS
  districts_r <- project(districts_vect, crs(r))
  
  # Obtain dates associated with raster layers
  dates <- time(r)
  
  if (is.null(dates)) {
    stop(
      "No time information found in rainfall file: ",
      basename(nc_file)
    )
  }
  
  dates <- as.Date(dates)
  
  # ISO year and epidemiological week
  Year <- lubridate::isoyear(dates)
  Week <- lubridate::isoweek(dates)
  
  # Extract mean daily rainfall per district
  extracted <- terra::extract(
    r,
    districts_r,
    fun = mean,
    na.rm = TRUE
  )
  
  extracted <- as.data.frame(extracted)
  
  # Remove extraction ID column
  rainfall_matrix <- as.matrix(extracted[, -1, drop = FALSE])
  
  # Group daily observations into weeks
  week_keys <- paste(Year, Week, sep = "-")
  
  unique_weeks <- unique(week_keys)
  
  weekly_results <- lapply(unique_weeks, function(k) {
    
    idx <- which(week_keys == k)
    
    # Weekly rainfall total
    rain_total <- rowSums(
      rainfall_matrix[, idx, drop = FALSE],
      na.rm = TRUE
    )
    
    # Weekly average daily rainfall
    rain_avg <- rowMeans(
      rainfall_matrix[, idx, drop = FALSE],
      na.rm = TRUE
    )
    
    # Preserve missing values if all daily values are missing
    all_missing <- rowSums(
      !is.na(rainfall_matrix[, idx, drop = FALSE])
    ) == 0
    
    rain_total[all_missing] <- NA_real_
    rain_avg[all_missing] <- NA_real_
    
    parts <- strsplit(k, "-", fixed = TRUE)[[1]]
    
    data.frame(
      ID_HD = districts_r$ID_HD,
      Year = as.integer(parts[1]),
      Week = as.integer(parts[2]),
      Rain_total = rain_total,
      Rain_avg = rain_avg
    )
  })
  
  bind_rows(weekly_results)
}


# Extract rainfall from all NetCDF files

rain_weekly_all <- bind_rows(
  lapply(
    precip_files,
    get_weekly_rain,
    districts_vect = districts_vect
  )
)


# Save extracted weekly rainfall
write_csv(
  rain_weekly_all,
  file.path(output_dir, "rain_weekly.csv")
)


# Integrate weekly rainfall into the main dataset

cholera_inla <- cholera_inla %>%
  left_join(
    rain_weekly_all %>%
      select(ID_HD, Year, Week, Rain_total, Rain_avg),
    by = c("ID_HD", "Year", "Week")
  )


# 5.2 Integrate temperature variables derived from ERA5 -------

# ERA5 daily 2 m temperature is extracted for each health
# district, converted from Kelvin to Celsius, and aggregated
# into weekly mean temperature.

temp_files <- list.files(
  temp_dir,
  pattern = "^[0-9]{4}_era5\\.nc$",
  full.names = TRUE
)

if (length(temp_files) == 0) {
  stop("No ERA5 files found. Check temp_dir and file names.")
}


# Read and combine ERA5 raster files
temp_rasters <- lapply(temp_files, rast)

temp_stack <- rast(temp_rasters)

# Convert Kelvin to degrees Celsius
temp_stack <- temp_stack - 273.15


# Obtain dates from raster time information.
# If time information is unavailable, supply the correct
# starting date for the first layer below.

temp_dates <- time(temp_stack)

if (is.null(temp_dates)) {
  
  # Change this date if the first ERA5 layer starts on
  # a different date.
  temp_start_date <- as.Date("2017-01-01")
  
  temp_dates <- seq(
    from = temp_start_date,
    by = "day",
    length.out = nlyr(temp_stack)
  )
  
} else {
  
  temp_dates <- as.Date(temp_dates)
  
}


# Ensure date and raster layer counts match
stopifnot(length(temp_dates) == nlyr(temp_stack))


# Reproject districts to ERA5 raster CRS
districts_temp <- project(
  districts_vect,
  crs(temp_stack)
)


# Extract daily mean temperature per district
temp_extracted <- terra::extract(
  temp_stack,
  districts_temp,
  fun = mean,
  na.rm = TRUE
)

temp_matrix <- as.matrix(
  temp_extracted[, -1, drop = FALSE]
)


# Create daily district-level temperature table
temp_daily <- data.frame(
  ID_HD = rep(
    districts_temp$ID_HD,
    times = length(temp_dates)
  ),
  
  Year = rep(
    lubridate::isoyear(temp_dates),
    each = nrow(districts_temp)
  ),
  
  Week = rep(
    lubridate::isoweek(temp_dates),
    each = nrow(districts_temp)
  ),
  
  Temp = as.vector(temp_matrix)
)


# Aggregate daily temperature to weekly mean
temp_weekly <- temp_daily %>%
  group_by(ID_HD, Year, Week) %>%
  summarise(
    Temp = if (all(is.na(Temp))) {
      NA_real_
    } else {
      mean(Temp, na.rm = TRUE)
    },
    .groups = "drop"
  )


# Save weekly temperature
write_csv(
  temp_weekly,
  file.path(output_dir, "temp_weekly.csv")
)


# Integrate weekly temperature into the main dataset

cholera_inla <- cholera_inla %>%
  left_join(
    temp_weekly %>%
      select(ID_HD, Year, Week, Temp),
    by = c("ID_HD", "Year", "Week")
  )


# 5.3 Derive flood exposure indicators ------------------------

cholera_inla <- cholera_inla %>%
  mutate(
    flood1 = ifelse(Rain_avg >= 50, 1, 0),
    flood2 = ifelse(Rain_total >= 200, 1, 0),
    flood3 = pmax(Rain_avg - 50, 0)
  ) %>%
  ungroup()


# 5.4 Topographical variables: DEM, slope and TWI -------------

# Topographical variables are derived from the SRTM DEM.
#
# Outputs:
#   Elevation = Mean elevation per health district
#   Slope     = Mean terrain slope per health district
#   TWI       = Mean Topographic Wetness Index per district
#
# WhiteboxTools is used for depression filling,
# D8 flow accumulation and TWI derivation.


# Load elevation raster
Elevation <- rast(dem_file)

# Reproject district boundaries to DEM CRS
districts_dem <- project(
  districts_vect,
  crs(Elevation)
)


# Create temporary raster files
dem_path <- file.path(output_dir, "dem.tif")
filled_dem_path <- file.path(output_dir, "dem_filled.tif")
flow_path <- file.path(output_dir, "flow_accum.tif")
slope_path <- file.path(output_dir, "slope.tif")
twi_path <- file.path(output_dir, "twi.tif")


# Save DEM for WhiteboxTools
writeRaster(
  Elevation,
  dem_path,
  overwrite = TRUE
)


# Initialise WhiteboxTools
wbt_init()


# 5.4.1 Fill depressions in the DEM

wbt_fill_depressions(
  dem = dem_path,
  output = filled_dem_path
)


# Load depression-filled DEM
dem_filled <- rast(filled_dem_path)


# 5.4.2 Calculate D8 flow accumulation

wbt_d8_flow_accumulation(
  input = filled_dem_path,
  output = flow_path
)


# 5.4.3 Derive slope

slope <- terrain(
  dem_filled,
  v = "slope",
  unit = "radians"
)

writeRaster(
  slope,
  slope_path,
  overwrite = TRUE
)


# 5.4.4 Derive Topographic Wetness Index (TWI)

wbt_wetness_index(
  sca = flow_path,
  slope = slope_path,
  output = twi_path
)

twi <- rast(twi_path)


# 5.4.5 Extract district-level topographical statistics

# Elevation
elevation_extract <- terra::extract(
  Elevation,
  districts_dem,
  fun = mean,
  na.rm = TRUE
)

elevation_mean <- data.frame(
  ID_HD = districts_dem$ID_HD,
  Elevation = elevation_extract[[2]]
)


# Slope
slope_extract <- terra::extract(
  slope,
  districts_dem,
  fun = mean,
  na.rm = TRUE
)

slope_mean <- data.frame(
  ID_HD = districts_dem$ID_HD,
  Slope = slope_extract[[2]]
)


# TWI
twi_extract <- terra::extract(
  twi,
  districts_dem,
  fun = mean,
  na.rm = TRUE
)

twi_mean <- data.frame(
  ID_HD = districts_dem$ID_HD,
  TWI = twi_extract[[2]]
)


# Combine topographical variables
topography_df <- elevation_mean %>%
  left_join(slope_mean, by = "ID_HD") %>%
  left_join(twi_mean, by = "ID_HD")


# Save district-level topographical variables
write_csv(
  topography_df,
  file.path(output_dir, "topography_district.csv")
)


# Integrate topographical variables into the main dataset

cholera_inla <- cholera_inla %>%
  left_join(
    topography_df %>%
      select(ID_HD, Elevation, Slope, TWI),
    by = "ID_HD"
  )


# 5.5 Integrate water exposure --------------------------------

water_exposure <- read_csv(
  "data/processed/water_exposure.csv"
)

cholera_inla <- cholera_inla %>%
  left_join(
    water_exposure %>%
      select(ID_HD, water_exposure_scaled),
    by = "ID_HD"
  )


# ============================================================
# 6. SOCIO-ECONOMIC AND INFRASTRUCTURE VARIABLES
# ============================================================


# 6.1 Integrate population density ----------------------------

pop_density <- read_csv(
  "data/processed/pop_density.csv"
)

cholera_inla <- cholera_inla %>%
  left_join(
    pop_density %>%
      select(ID_HD, Pop_density),
    by = "ID_HD"
  ) %>%
  mutate(
    Pop_density_sc = as.numeric(scale(Pop_density))
  )


# 6.2 Integrate building density -----------------------------

building_density <- read_csv(
  "data/processed/building_density.csv"
)

cholera_inla <- cholera_inla %>%
  left_join(
    building_density %>%
      select(ID_HD, building_density),
    by = "ID_HD"
  ) %>%
  mutate(
    building_density_sc = as.numeric(scale(building_density))
  )


# 6.3 Integrate poverty index ---------------------------------

poverty <- read_csv(
  "data/processed/poverty_index.csv"
)

cholera_inla <- cholera_inla %>%
  left_join(
    poverty %>%
      select(ID_HD, poverty_index),
    by = "ID_HD"
  ) %>%
  mutate(
    poverty_index_sc = as.numeric(scale(poverty_index))
  )


# 6.4 Integrate open defecation -------------------------------

open_defecation <- read_csv(
  "data/processed/open_defecation.csv"
)

cholera_inla <- cholera_inla %>%
  left_join(
    open_defecation %>%
      select(ID_HD, open_defecation),
    by = "ID_HD"
  )


# 6.5 Integrate poor sanitation -------------------------------

wash_data <- read_csv(
  "data/processed/wash.data.csv"
)

cholera_inla <- cholera_inla %>%
  left_join(
    wash_data %>%
      select(ID_HD, Poor_sanitation),
    by = "ID_HD"
  ) %>%
  mutate(
    Poor_sanitation_sc = as.numeric(scale(Poor_sanitation))
  )


# 6.6 Integrate healthcare facilities -------------------------

hcf <- read_csv(
  "data/processed/hcf.csv"
)

cholera_inla <- cholera_inla %>%
  left_join(
    hcf %>%
      select(ID_HD, health_facilities_per1000),
    by = "ID_HD"
  )


# ============================================================
# 7. Save integrated modelling dataset
# ============================================================

write_csv(
  cholera_inla,
  "data/processed/cholera_inla_covariates.csv"
)


# 8. Final dataset checks -------------------------------------

# Dataset dimensions
dim(cholera_inla)

# Variable names
names(cholera_inla)

# Missing values in integrated variables
colSums(is.na(cholera_inla))

# Summary of selected study variables
summary(
  cholera_inla %>%
    select(
      y_casecount,
      incidence_rate,
      SMR,
      Rain_total,
      Rain_avg,
      Temp,
      flood1,
      flood2,
      flood3,
      Elevation,
      Slope,
      TWI,
      water_exposure_scaled,
      Pop_density,
      building_density,
      poverty_index,
      open_defecation,
      Poor_sanitation,
      health_facilities_per1000
    )
)
