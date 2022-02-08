library(tidyverse)
library(sf)
library(terra)

# load data -------------------------
# load tetradensity + mass dataset
p <- st_read("data/interim/tetradensity_masses.shp")
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

# check extrapolations -----------------
ranges <- d %>% 
  select(NPP, Temp, PCV) %>% 
  distinct_all() %>% 
  pivot_longer(cols = everything(), names_to = "var", values_to = "vals") %>% 
  arrange(var, vals) %>% 
  left_join(
    tibble(var = c("PCV", "Temp", "NPP"),
           EU.min = c(pcv@ptr$range_min, temp@ptr$range_min, npp@ptr$range_min),
           EU.max = c(pcv@ptr$range_max, temp@ptr$range_max, npp@ptr$range_max))
  )

ranges %>% 
  ggplot() +
  aes(var, vals) +
  geom_errorbar(aes(x = var, ymin = EU.min, ymax = EU.max), col = "steelblue") +
  geom_violin() +
  geom_jitter(width = .2) +
  facet_wrap(~var, scales = "free") +
  xlab("Environmental layer") +
  ylab("") +
  theme_light() +
  theme(axis.text.x = element_blank(),
        panel.grid.major.x = element_blank())

# visualize marginal respose to assess extrapolation issues --------------
# NPP
par(mfrow = c(2, 2), mar = c(4, 2, 2, 2))
for(cl in names(lm_density)) {
  m <- lm_density[[cl]]
  x <- 10 ^ m$model$`log10(NPP)`
  y <- 10 ^ (m$model$`log10(NPP)` * m$coefficients[[4]] + m$model$`log10(NPP)` ^ 2 * m$coefficients[[5]])
  plot(ranges %>% filter(var == "NPP") %>% select(EU.min, EU.max) %>% distinct_all() %>% as.matrix() %>% as.vector(),
       range(y),
       frame = FALSE, col = "white", main = cl, xlab = "NPP", ylab = "Density")
  rug(ranges %>% filter(var == "NPP") %>% select(EU.min, EU.max) %>% distinct_all() %>% as.matrix() %>% as.vector(),
      lwd = 3, col = "tomato")
  points(x, y, pch = 20, col = "steelblue")
}

# PCV
par(mfrow = c(2, 2), mar = c(4, 2, 2, 2))
for(cl in names(lm_density)) {
  m <- lm_density[[cl]]
  x <- m$model$PCV
  y <- 10 ^ (x * m$coefficients[[6]] + x ^ 2 * m$coefficients[[7]])
  plot(ranges %>% filter(var == "PCV") %>% select(EU.min, EU.max) %>% distinct_all() %>% as.matrix() %>% as.vector(),
       range(y),
       frame = FALSE, col = "white", main = cl, xlab = "PCV", ylab = "Density")
  rug(ranges %>% filter(var == "PCV") %>% select(EU.min, EU.max) %>% distinct_all() %>% as.matrix() %>% as.vector(),
      lwd = 3, col = "tomato")
  points(x, y, pch = 20, col = "steelblue")
}

# Temperature
par(mfrow = c(2, 2), mar = c(4, 2, 2, 2))
for(cl in names(lm_density)) {
  m <- lm_density[[cl]]
  x <- m$model$Temp
  y <- 10 ^ (x * m$coefficients[[9]])
  plot(ranges %>% filter(var == "Temp") %>% select(EU.min, EU.max) %>% distinct_all() %>% as.matrix() %>% as.vector(),
       range(y),
       frame = FALSE, col = "white", main = cl, xlab = "Temperature", ylab = "Density")
  rug(ranges %>% filter(var == "Temp") %>% select(EU.min, EU.max) %>% distinct_all() %>% as.matrix() %>% as.vector(),
      lwd = 3, col = "tomato")
  points(x, y, pch = 20, col = "steelblue")
}

# check on geographic coordinates -----------
npp[npp < min(d$NPP)] <- -9999
npp[npp > max(d$NPP)] <- -9999
npp[npp > 0] <- NA
plot(npp)

temp[temp < min(d$Temp)] <- -9999
temp[temp > max(d$Temp)] <- -9999
temp[temp > -1000] <- NA
plot(temp)

pcv[pcv < min(d$PCV)] <- -9999
pcv[pcv > max(d$PCV)] <- -9999
pcv[pcv > 0] <- NA
plot(pcv)

eu <- rast("For_the_network_symposium/Environmental data/CHELSA_EUR11_npp_norm_1981-2005_V1.1.nc")
eu <- aggregate(eu, 120, na.rm = TRUE)
eu[!is.na(eu)] <- 1
plot(eu, col = "grey65")
plot(npp, add = TRUE, col = adjustcolor("green3", alpha.f = .5))
plot(temp, add = TRUE, col = adjustcolor("tomato", alpha.f = .5))
plot(pcv, add = TRUE, col = adjustcolor("steelblue", alpha.f = .5))
