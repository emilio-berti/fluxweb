library(tidyverse)
library(sf)
library(terra)

# load data -------------------------
# load tetradensity + mass dataset
p <- st_read("../data/interim/tetradensity-masses.shp")
# load environmental rasters
npp <- rast("../../For_the_network_symposium/Environmental data/CHELSA_EUR11_npp_norm_1981-2005_V1.1.nc")
pcv <- rast("../../For_the_network_symposium/Environmental data/CHELSA_EUR11_pr_timstd_norm_1981-2005_V1.1.nc")
temp <- rast("../../For_the_network_symposium/Environmental data/CHELSA_EUR11_tas_timmean_1981-2005_V1.1.nc")
if (!st_crs(p) == st_crs(npp)) {
   warning("CRS are not the same - reprojecting")
}
npp <- aggregate(npp, 1 / res(npp)[1], fun = "sum", na.rm = TRUE)
pcv <- aggregate(pcv, 1 / res(pcv)[1], fun = "median", na.rm = TRUE)
temp <- aggregate(temp, 1 / res(temp)[1], fun = "median", na.rm = TRUE)
pcv <- mask(pcv, npp)

# extract values for tetradensity locations
npp_vals <- extract(npp, st_coordinates(p))[[1]]
temp_vals <- extract(temp, st_coordinates(p))[[1]]
pcv_vals <- extract(pcv, st_coordinates(p))[[1]]

# add the values to the shapefile
p <- p %>% 
   mutate(NPP = npp_vals,
          Temp = temp_vals,
          PCV = pcv_vals) %>%
   drop_na()

# drop geometry for modeling; this creates a dataframe
d <- p %>% 
   st_drop_geometry() %>% 
   as_tibble() %>% 
   mutate(across(where(is.factor), ~as.character(.x)))

# make models ----------------------
# we selected models based on AIC. Run check_models.R
# from here to see how these were selected, if interested.

lm_amph <- lm(log10(Density) ~ log10(Mass) + I(log10(Mass) ^ 2) +
                 log10(NPP) + I(log10(NPP) ^ 2) +
                 PCV + I(PCV ^ 2) +
                 Order + Temp,
              data = d %>% filter(Class == "Amphibia"))

lm_bird <- lm(log10(Density) ~ log10(Mass) + I(log10(Mass) ^ 2) +
                 log10(NPP) + I(log10(NPP) ^ 2) +
                 PCV + I(PCV ^ 2) +
                 Order + Temp,
              data = d %>% filter(Class == "Aves"))

lm_mamm <- lm(log10(Density) ~ log10(Mass) + I(log10(Mass) ^ 2) +
                 log10(NPP) + I(log10(NPP) ^ 2) +
                 PCV + I(PCV ^ 2) +
                 Order + Temp, 
              data = d %>% filter(Class == "Mammalia"))

lm_rept <- lm(log10(Density) ~ log10(Mass) + I(log10(Mass) ^ 2) +
                 log10(NPP) + I(log10(NPP) ^ 2) +
                 PCV + I(PCV ^ 2) +
                 Order + Temp,
              data = d %>% filter(Class == "Reptilia"))

# list of models
lm_density <- list("Amphibia" = lm_amph,
                   "Aves" = lm_bird,
                   "Mammalia" = lm_mamm,
                   "Reptilia" = lm_rept)

# calculate the environmental part of lms ------
#amphibians
amph_env <- lm_density$Amphibia$coefficients[["(Intercept)"]] +
   log10(npp) * lm_density$Amphibia$coefficients[["log10(NPP)"]] +
   (log10(npp) ^ 2) * lm_density$Amphibia$coefficients[["I(log10(NPP)^2)"]] +
   pcv * lm_density$Amphibia$coefficients[["PCV"]] +
   pcv ^ 2 * lm_density$Amphibia$coefficients[["I(PCV^2)"]] +
   temp * lm_density$Amphibia$coefficients[["Temp"]]
#birds
bird_env <- lm_density$Aves$coefficients[["(Intercept)"]] +
   log10(npp) * lm_density$Aves$coefficients[["log10(NPP)"]] +
   (log10(npp) ^ 2) * lm_density$Aves$coefficients[["I(log10(NPP)^2)"]] +
   pcv * lm_density$Aves$coefficients[["PCV"]] +
   pcv ^ 2 * lm_density$Aves$coefficients[["I(PCV^2)"]] +
   temp * lm_density$Aves$coefficients[["Temp"]]
#mammals
mamm_env <- lm_density$Mammalia$coefficients[["(Intercept)"]] +
   log10(npp) * lm_density$Mammalia$coefficients[["log10(NPP)"]] +
   (log10(npp) ^ 2) * lm_density$Mammalia$coefficients[["I(log10(NPP)^2)"]] +
   pcv * lm_density$Mammalia$coefficients[["PCV"]] +
   pcv ^ 2 * lm_density$Mammalia$coefficients[["I(PCV^2)"]] +
   temp * lm_density$Mammalia$coefficients[["Temp"]]
#reptiles
rept_env <- lm_density$Reptilia$coefficients[["(Intercept)"]] +
   log10(npp) * lm_density$Reptilia$coefficients[["log10(NPP)"]] +
   (log10(npp) ^ 2) * lm_density$Reptilia$coefficients[["I(log10(NPP)^2)"]] +
   pcv * lm_density$Reptilia$coefficients[["PCV"]] +
   pcv ^ 2 * lm_density$Reptilia$coefficients[["I(PCV^2)"]] +
   temp * lm_density$Reptilia$coefficients[["Temp"]]

# create tables for the partial effects
px_id <- read_csv("../../For_the_network_symposium/perpixel_envVar.csv") %>% 
   dplyr::select(x, y, PxID, PAGENAME, CellNumber)

template <- rast("../../For_the_network_symposium/10k_mask2/reference_grid_10km.img")

amph_env <- project(amph_env, template, method = "near")
bird_env <- project(bird_env, template, method = "near")
mamm_env <- project(mamm_env, template, method = "near")
rept_env <- project(rept_env, template, method = "near")

xy <- xyFromCell(rept_env, 1:ncell(rept_env))

env <- px_id %>% 
   left_join(
      tibble(x = xy[, 1],
             y = xy[, 2],
             amphibia = values(amph_env)[, 1],
             aves = values(bird_env)[, 1],
             mammalia = values(mamm_env)[, 1],
             reptilia = values(rept_env)[, 1]) %>% 
         drop_na()
   )
# remove all newly created variables except the linear models
rm(npp_vals, temp_vals, pcv_vals,
   npp, temp, pcv,
   lm_amph, lm_bird, lm_mamm, lm_rept,
   d, p)
null <- gc(verbose = "FALSE")
rm(null)
