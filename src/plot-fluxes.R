library(tidyverse)
library(raster)
library(tmap)
library(sf)
library("rnaturalearth")
library("rnaturalearthdata")

# I don't think we need to load this again, as the table key is part of the 
# environmental variables table and can be extracted from there.
# maskID <- read.dbf("../10k_mask2/reference_grid_10km.img.vat.dbf") %>% 
#   as_tibble() %>% 
#   arrange(Value)
path <- "/data/tib-group-share/FutureWeb/For_the_network_symposium/"
m <- raster(file.path(path, "/10k_mask2/reference_grid_10km.img"))
# rich <- raster("../species-richness.tif")
fluxes <- read_csv("data/final/pixel-fluxes.csv", show_col_types = FALSE) 

px_id <- read_csv("../For_the_network_symposium/perpixel_envVar.csv", show_col_types = FALSE) %>%  #this is basically the key table
  dplyr::select(PAGENAME, CellNumber, x, y) %>% 
  filter(PAGENAME %in% fluxes$PAGENAME) %>% #keep only ids of cells for which we have fluxes
  left_join(fluxes) %>% # tidymagic happening
  mutate(Sum_of_fluxes = ifelse(Sum_of_fluxes == 0, NA, Sum_of_fluxes), #set 0 to NAs
         Avg_in_fluxes = ifelse(Avg_in_fluxes == 0, NA, Avg_in_fluxes),
         Avg_out_fluxes = ifelse(Avg_out_fluxes == 0, NA, Avg_out_fluxes))

# interpolate missing values ---------
# interp_sum <- mgcv::gam(log10(Sum_of_fluxes) ~ s(x) + s(y) + s(Richness), data = px_id)
# interp_in <- mgcv::gam(log10(Avg_in_fluxes) ~ s(x) + s(y) + s(Richness), data = px_id)
# interp_out <- mgcv::gam(log10(Avg_out_fluxes) ~ s(x) + s(y) + s(Richness), data = px_id)

# some shapes to plot nicely --------------
countries <- ne_countries(scale = "medium", returnclass = "sf") %>% 
  filter(continent == "Europe") %>% 
  st_transform(crs = crs(m))

grid <- st_graticule(countries)

fl <- rast(m)
crs(fl) <- crs(m)@projargs
fl[!is.na(fl)] <- NA
fl[px_id$CellNumber] <- px_id$Sum_of_fluxes
fl <- log10(fl)
plot(fl)

# plot fluxes ---------------
tmap_mode("plot")
pal <- wesanderson::wes_palette("Zissou1", 100, type = "continuous") #this is just a palette I like
# pal <- RColorBrewer::brewer.pal(9, "RdBu") #this is just a palette I like
# pal <- colorRampPalette(pal)(100)
# pal[1] <- "grey20"
bk <- m
bk[!is.na(bk)] <- 1

down <- quantile(values(fl), na.rm = TRUE, .05)
top <- quantile(values(fl), na.rm = TRUE, .95)
fl[fl < down] <- down
fl[fl > top] <- top
fl <- raster::raster(fl)

# map <- 
tm_shape(bk) + 
  tm_raster(palette = c("grey20"), legend.show = FALSE) +
  tm_shape(fl) + 
  tm_raster(
  palette = viridis::inferno(100),
  style = "cont", #continuous palette
  title = expression("Sum of fluxes (log"[10]*")"),
  legend.format = list(scientific = FALSE, digits = 2),
  # legend.is.portrait = FALSE #horizontal legend,
  legend.reverse = TRUE, 
  breaks = seq(3.5, 6, length.out = 9)
)

tmap_save(filename = "docs/plots/flux-figure.png", width = 6, height = 6, dpi = 600)

# you can do many things with the table -----------
fluxes %>% 
  slice_sample(prop = .05) %>%
  dplyr::select(Temperature:Avg_in_fluxes) %>% 
  dplyr::select(-Biomasses, -Density) %>% 
  filter(Avg_in_fluxes > -1e3, Avg_out_fluxes > -1e3) %>% 
  mutate(Avg_in_fluxes = - Avg_in_fluxes,
         Avg_out_fluxes = - Avg_out_fluxes) %>% 
  pivot_longer(cols = Sum_of_fluxes:Avg_in_fluxes,
               names_to = "what", values_to = "flux") %>%
  drop_na() %>% 
  ggplot() +
  aes(Temperature, flux) + 
  geom_point(alpha = .1) +
  geom_smooth() +
  facet_wrap(~what, scales = "free") +
  scale_y_log10() +
  theme_bw()

m <- fluxes %>% 
  dplyr::select(Temperature:Avg_in_fluxes) %>% 
  dplyr::select(-Biomasses, -Density) %>% 
  filter(Avg_in_fluxes > -1e3, Avg_out_fluxes > -1e3) %>% 
  mutate(Avg_in_fluxes = - Avg_in_fluxes,
         Avg_out_fluxes = - Avg_out_fluxes,
         flux = log10(Sum_of_fluxes)) %>%
  mutate(across(where(is.numeric), ~scale(.x)[, 1])) %>% 
  lm(flux ~ Number_of_species + Temperature, data = .)
plot(m, which = 2)
car::Anova(m, type = "III")
parameters::parameters(m)
