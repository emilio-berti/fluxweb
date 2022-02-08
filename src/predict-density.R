library(tidyverse)
library(sf)
library(terra)

source("make_density_models.R")

# futureweb
futureweb <- read_delim("data/raw/futureweb.csv", delim = ";", show_col_types = FALSE) %>% 
  select(Species, `Species ID`)
# backbone
backbone <- read_csv("data/interim/backbone-masses.csv", show_col_types = FALSE)

futureweb %>% 
  left_join(backbone) %>% 
  filter(!is.na(Class)) %>% 
  mutate(missing = ifelse(is.na(Mass), TRUE, FALSE)) %>% 
  group_by(Class, missing) %>% 
  tally() %>%
  pivot_wider(names_from = Class, values_from = n)

backbone %>% 
  filter(Class == "Amphibia") %>% 
  select(Order, Family, Mass) %>% 
  ggplot() +
  aes(Family, Mass) +
  geom_boxplot() +
  scale_y_log10()

# load environmental rasters
npp <- rast("For_the_network_symposium/Environmental data/CHELSA_EUR11_npp_norm_1981-2005_V1.1.nc")
pcv <- rast("For_the_network_symposium/Environmental data/CHELSA_EUR11_pr_timstd_norm_1981-2005_V1.1.nc")
temp <- rast("For_the_network_symposium/Environmental data/CHELSA_EUR11_tas_timmean_1981-2005_V1.1.nc")
if (!st_crs(p) == st_crs(npp)) {
  warning("CRS are not the same - reprojecting")
}
npp <- aggregate(npp, 1 / res(npp)[1], fun = "sum", na.rm = TRUE)
pcv <- aggregate(pcv, 1 / res(pcv)[1], fun = "median", na.rm = TRUE)
temp <- aggregate(temp, 1 / res(temp)[1], fun = "median", na.rm = TRUE)
pcv <- mask(pcv, npp)

env <- tibble(`log10(NPP)` = log10(values(npp)[, 1]),
              PCV = values(pcv)[, 1],
              Temp = values(temp)[, 1]) %>% 
  drop_na()

predict(lm_density$Amphibia, env)
