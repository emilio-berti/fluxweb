library(tidyverse)
library(sf)
library(terra)
library(MuMIn)

if (!"path" %in% ls()) path <- "/data/tib-group-share/FutureWeb/For_the_network_symposium/"

# load data -------------------------
template <- rast("../For_the_network_symposium/10k_mask2/reference_grid_10km.img")
proj4 <- crs(template, proj4 = TRUE)
# load tetradensity + mass dataset
p <- st_read("data/interim/tetradensity-masses.shp") %>% st_transform(proj4)
# load environmental rasters
npp <- rast(file.path(path, "Environmental data/CHELSA_EUR11_npp_norm_1981-2005_V1.1.nc"))
pcv <- rast(file.path(path, "Environmental data/CHELSA_EUR11_pr_timstd_norm_1981-2005_V1.1.nc"))
temp <- rast(file.path(path, "Environmental data/CHELSA_EUR11_tas_timmean_1981-2005_V1.1.nc"))
if (!st_crs(p) == st_crs(npp)) {
   warning("CRS are not the same - reprojecting")
}
# npp <- aggregate(npp, 1 / res(npp)[1], fun = "median", na.rm = TRUE) %>% project(proj4)
# pcv <- aggregate(pcv, 1 / res(pcv)[1], fun = "median", na.rm = TRUE) %>% project(proj4)
# temp <- aggregate(temp, 1 / res(temp)[1], fun = "median", na.rm = TRUE) %>% project(proj4)
npp <- npp %>% project(template) %>% mask(template)
pcv <- pcv %>% project(template) %>% mask(template)
temp <- temp %>% project(template) %>% mask(template)


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
options(na.action = "na.fail")

sub_d <- d %>% filter(Class == "Amphibia") %>% as.data.frame()
lm_amph <- lm(log10(Density) ~ 0 + log10(Mass) + I(log10(Mass) ^ 2) +
                 log10(NPP) + I(log10(NPP) ^ 2) +
                 PCV + I(PCV ^ 2) +
                 Order + Temp,
              data = sub_d)
coef_amph <- model.avg(dredge(lm_amph), delta < 4)$coefficients[1, ]

sub_d <- d %>% filter(Class == "Aves") %>% as.data.frame()
lm_bird <- lm(log10(Density) ~ 0 + log10(Mass) + I(log10(Mass) ^ 2) +
                 log10(NPP) + I(log10(NPP) ^ 2) +
                 PCV + I(PCV ^ 2) +
                 Order + Temp,
              data = sub_d)
coef_bird <- model.avg(dredge(lm_bird), delta < 4)$coefficients[1, ]

sub_d <- d %>% filter(Class == "Mammalia") %>% as.data.frame()
lm_mamm <- lm(log10(Density) ~ 0 + log10(Mass) + I(log10(Mass) ^ 2) +
                 log10(NPP) + I(log10(NPP) ^ 2) +
                 PCV + I(PCV ^ 2) +
                 Order + Temp, 
              data = sub_d)
coef_mamm <- model.avg(dredge(lm_mamm), delta < 4)$coefficients[1, ]

sub_d <- d %>% filter(Class == "Reptilia") %>% as.data.frame()
lm_rept <- lm(log10(Density) ~ 0 + log10(Mass) + I(log10(Mass) ^ 2) +
                 log10(NPP) + I(log10(NPP) ^ 2) +
                 PCV + I(PCV ^ 2) +
                 Order + Temp,
              data = sub_d)
coef_rept <- model.avg(dredge(lm_rept), delta < 4)$coefficients[1, ]

lm_density <- list(coef_amph, coef_bird, coef_mamm, coef_rept)
names(lm_density) <- c("Amphibia", "Aves", "Mammalia", "Reptilia")

# calculate the environmental part of lms --------
# set cells outside training range to NAs

#amphibians
amph_env <- log10(npp) * lm_density$Amphibia[["log10(NPP)"]] +
   (log10(npp) ^ 2) * lm_density$Amphibia[["I(log10(NPP)^2)"]] + 
   pcv * lm_density$Amphibia[["PCV"]] + 
   pcv ^ 2 * lm_density$Amphibia[["I(PCV^2)"]] +
   temp * lm_density$Amphibia[["Temp"]]
#birds
bird_env <- log10(npp) * lm_density$Aves[["log10(NPP)"]] +
   (log10(npp) ^ 2) * lm_density$Aves[["I(log10(NPP)^2)"]] +
   pcv * lm_density$Aves[["PCV"]] +
   pcv ^ 2 * lm_density$Aves[["I(PCV^2)"]] +
   temp * lm_density$Aves[["Temp"]]
#mammals
mamm_env <- log10(npp) * lm_density$Mammalia[["log10(NPP)"]] +
   (log10(npp) ^ 2) * lm_density$Mammalia[["I(log10(NPP)^2)"]] +
   pcv * lm_density$Mammalia[["PCV"]] +
   pcv ^ 2 * lm_density$Mammalia[["I(PCV^2)"]] +
   temp * lm_density$Mammalia[["Temp"]]
#reptiles
rept_env <- log10(npp) * lm_density$Reptilia[["log10(NPP)"]] +
   (log10(npp) ^ 2) * lm_density$Reptilia[["I(log10(NPP)^2)"]] + 
   pcv * lm_density$Reptilia[["PCV"]] +
   pcv ^ 2 * lm_density$Reptilia[["I(PCV^2)"]] +
   temp * lm_density$Reptilia[["Temp"]]

# create tables for the partial effects
px_id <- read_csv("../For_the_network_symposium/perpixel_envVar.csv", show_col_types = FALSE) %>% 
   dplyr::select(x, y, PxID, PAGENAME, CellNumber)

cell_size <- prod(res(npp) / 1000)
small_scale_area <- prod(res(template) / 1000)

xy <- xyFromCell(rept_env, 1:ncell(rept_env))

px_id <- read_csv("../For_the_network_symposium/perpixel_envVar.csv", show_col_types = FALSE) %>% 
   dplyr::select(x, y, PxID, PAGENAME, CellNumber)

env <- px_id %>% 
   mutate(Amphibia = extract(amph_env, CellNumber)[, 1],
          Aves = extract(bird_env, CellNumber)[, 1],
          Mammalia = extract(mamm_env, CellNumber)[, 1],
          Reptilia = extract(rept_env, CellNumber)[, 1])

# get areas of extrapolations ----------------
# amphibians
amph_extr <- c(npp, pcv, temp)
amph_extr[!is.na(amph_extr)] <- 0
extr <- npp < min(10 ^ lm_amph$model$`log10(NPP)`) | npp > max(10 ^ lm_amph$model$`log10(NPP)`)
amph_extr[[1]] <- amph_extr[[1]] + extr
extr <- pcv < min(lm_amph$model$PCV) | pcv > max(lm_amph$model$PCV)
amph_extr[[2]] <- amph_extr[[2]] + extr
extr <- temp < min(lm_amph$model$Temp) | temp > max(lm_amph$model$Temp)
amph_extr[[3]] <- amph_extr[[3]] + extr
amph_extr <- sum(amph_extr)

# birds
bird_extr <- c(npp, pcv, temp)
bird_extr[!is.na(bird_extr)] <- 0
extr <- npp < min(10 ^ lm_bird$model$`log10(NPP)`) | npp > max(10 ^ lm_bird$model$`log10(NPP)`)
bird_extr[[1]] <- bird_extr[[1]] + extr
extr <- pcv < min(lm_bird$model$PCV) | pcv > max(lm_bird$model$PCV)
bird_extr[[2]] <- bird_extr[[2]] + extr
extr <- temp < min(lm_bird$model$Temp) | temp > max(lm_bird$model$Temp)
bird_extr[[3]] <- bird_extr[[3]] + extr
bird_extr <- sum(bird_extr)

# mammals
mamm_extr <- c(npp, pcv, temp)
mamm_extr[!is.na(mamm_extr)] <- 0
extr <- npp < min(10 ^ lm_mamm$model$`log10(NPP)`) | npp > max(10 ^ lm_mamm$model$`log10(NPP)`)
mamm_extr[[1]] <- mamm_extr[[1]] + extr
extr <- pcv < min(lm_mamm$model$PCV) | pcv > max(lm_mamm$model$PCV)
mamm_extr[[2]] <- mamm_extr[[2]] + extr
extr <- temp < min(lm_mamm$model$Temp) | temp > max(lm_mamm$model$Temp)
mamm_extr[[3]] <- mamm_extr[[3]] + extr
mamm_extr <- sum(mamm_extr)

# reptiles
rept_extr <- c(npp, pcv, temp)
rept_extr[!is.na(rept_extr)] <- 0
extr <- npp < min(10 ^ lm_rept$model$`log10(NPP)`) | npp > max(10 ^ lm_rept$model$`log10(NPP)`)
rept_extr[[1]] <- rept_extr[[1]] + extr
extr <- pcv < min(lm_rept$model$PCV) | pcv > max(lm_rept$model$PCV)
rept_extr[[2]] <- rept_extr[[2]] + extr
extr <- temp < min(lm_rept$model$Temp) | temp > max(lm_rept$model$Temp)
rept_extr[[3]] <- rept_extr[[3]] + extr
rept_extr <- sum(rept_extr)

extr <- c(amph_extr, bird_extr, mamm_extr, rept_extr)
names(extr) <- c("Amphibia", "Aves", "Mammalia", "Reptilia")

#write to file
writeRaster(extr, "data/final/extrapolations.tif", datatype = "INTU1", overwrite = TRUE)

# to plot, not needed for analyses
# library(tmap)
# pal <- colorRampPalette(c("grey20", "steelblue", "gold", "tomato"))
# a <- sum(extr)
# a <- trim(a)
# a <- raster::raster(a)
# tm_shape(a) +
#   tm_graticules() +
#   tm_raster(title = "Number of extrapolations",
#             style = "cat",
#             palette = pal(12), breaks = seq(1, 12),
#             legend.reverse = TRUE) +
#   tm_layout(legend.outside = TRUE) -> tm
# 
# tmap_save(tm, "docs/plots/extrapolation-map.png", width = 6, height = 4, dpi = 600)

# remove all newly created variables except the linear models ------
rm(npp_vals, temp_vals, pcv_vals,
   npp, temp, pcv,
   lm_amph, lm_bird, lm_mamm, lm_rept,
   coef_amph, coef_bird, coef_mamm, coef_rept,
   amph_env, bird_env, mamm_env,  rept_env,
   template, px_id, sub_d,
   d, p, xy,
   amph_extr, bird_extr, mamm_extr, rept_extr)
null <- gc(verbose = "FALSE")
rm(null)
