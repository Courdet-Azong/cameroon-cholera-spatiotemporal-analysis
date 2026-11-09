# ============================================================
# Spatiotemporal Modelling of Cholera in Cameroon
# Script: 04_bym_model_comparison.R
#
# Purpose:
# Fit and compare three Bayesian spatial and spatiotemporal
# models using the INLA framework.
#
# Model 1: Baseline spatial BYM model
# Model 2: Spatial BYM with temporal random effects
# Model 3: Spatial BYM with temporal effects and covariates
#
# Model evaluation: DIC, WAIC and LCPO
# ============================================================


# ============================================================
# 1. SPATIAL GRAPH CONSTRUCTION
# ============================================================

districts_sf <- st_as_sf(districts_vect)

nb <- poly2nb(districts_sf)

nb2INLA("districts.graph", nb)

districts.adj <- paste0(getwd(), "/districts.graph")


# ============================================================
# 2. MODEL 1: BASELINE SPATIAL BYM
# ============================================================
# This model accounts for spatial variation across health
# districts using structured and unstructured random effects.
# A year-specific IID intercept is also included.

model_baseline <- inla(
  
  y_casecount ~
    f(yrid_alpha, model = "iid") +
    
    offset(log(E_district)) +
    
    f(id_U,
      model = "bym",
      graph = districts.adj,
      scale.model = TRUE,
      group = yrid_U,
      control.group = list(model = "iid")),
  
  family = "poisson",
  data = cholera_inla,
  control.predictor = list(compute = TRUE),
  control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE)
)


# Model 1: Fit statistics

DIC <- round(model_baseline$dic$dic, 2)
WAIC <- round(model_baseline$waic$waic, 2)

LCPO <- round(
  -sum(log(model_baseline$cpo$cpo), na.rm = TRUE),
  2
)

cat(
  "Model 1: Baseline BYM\n",
  "DIC :", DIC, "\n",
  "WAIC:", WAIC, "\n",
  "LCPO:", LCPO, "\n"
)


# ============================================================
# 3. MODEL 2: BYM + TEMPORAL RANDOM EFFECTS
# ============================================================
# This model extends the baseline BYM specification by
# incorporating structured RW1 and unstructured IID
# temporal random effects.

model_random <- inla(
  
  y_casecount ~
    f(yrid_alpha, model = "iid") +
    
    offset(log(E_district)) +
    
    f(id_U,
      model = "bym",
      graph = districts.adj,
      scale.model = TRUE,
      group = yrid_U,
      control.group = list(model = "iid")) +
    
    f(id_gamma,
      model = "rw1",
      group = rid_gamma,
      control.group = list(model = "iid")) +
    
    f(id_delta,
      model = "iid",
      group = rid_delta,
      control.group = list(model = "iid")),
  
  family = "poisson",
  data = cholera_inla,
  control.predictor = list(compute = TRUE),
  control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE)
)


# Model 2: Fit statistics

DIC <- round(model_random$dic$dic, 2)
WAIC <- round(model_random$waic$waic, 2)

LCPO <- round(
  -sum(log(model_random$cpo$cpo), na.rm = TRUE),
  2
)

cat(
  "Model 2: BYM + Temporal Effects\n",
  "DIC :", DIC, "\n",
  "WAIC:", WAIC, "\n",
  "LCPO:", LCPO, "\n"
)


# ============================================================
# 4. MODEL 3: BYM + TEMPORAL EFFECTS + COVARIATES
# ============================================================
# This model extends the spatiotemporal specification by
# incorporating environmental and socioeconomic covariates.
# Nonlinear effects of lagged rainfall, temperature and flood
# variables are represented using RW1 functions.


# Group covariates for RW1 effects

cholera_inla <- cholera_inla %>%
  mutate(
    Rain_rw1 = inla.group(Rain_lag2, n = 20),
    Temp_rw1 = inla.group(Temp_lag2, n = 20),
    # Hosp_rw1 = inla.group(hospital_rate_sc, n = 20),
    flood_rw1 = inla.group(flood_lag2, n = 20)
    # Slope_rw1 = inla.group(Slope, n = 20),
    # Watexp_rw1 = inla.group(water_exposure_scaled, n = 20)
  )

cholera_inla$rid_flood


# Fit spatiotemporal covariate model

model_st <- inla(
  
  y_casecount ~
    f(yrid_alpha, model = "iid") +
    
    Poor_sanitation +
    poverty_index +
    slope_percent +
    Elevation +
    # TWI +
    water_prop +
    # pop_within1km_water +
    health_facilities_per1000 +
    # Pop_density +
    conflict_events_weekly +
    
    f(Rain_rw1,
      model = "rw1",
      group = rid_Prec,
      control.group = list(model = "iid")) +
    
    f(flood_rw1,
      model = "rw1",
      group = rid_flood,
      control.group = list(model = "iid")) +
    
    f(Temp_rw1,
      model = "rw1",
      group = rid_Temp,
      control.group = list(model = "iid")) +
    
    offset(log(E_district)) +
    
    f(id_U,
      model = "bym",
      graph = districts.adj,
      scale.model = TRUE,
      group = yrid_U,
      control.group = list(model = "iid")) +
    
    f(id_gamma,
      model = "rw1",
      group = rid_gamma,
      control.group = list(model = "iid")) +
    
    f(id_delta,
      model = "iid",
      group = rid_delta,
      control.group = list(model = "iid")),
  
  family = "poisson",
  data = cholera_inla,
  control.predictor = list(compute = TRUE),
  control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE)
)


# Model 3: Fit statistics

DIC <- round(model_st$dic$dic, 2)
WAIC <- round(model_st$waic$waic, 2)

LCPO <- round(
  -sum(log(model_st$cpo$cpo), na.rm = TRUE),
  2
)

cat(
  "Model 3: BYM + Temporal Effects + Covariates\n",
  "DIC :", DIC, "\n",
  "WAIC:", WAIC, "\n",
  "LCPO:", LCPO, "\n"
)


# ============================================================
# END OF MODEL FITTING
# ============================================================
