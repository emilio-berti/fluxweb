library(tidyverse)
library(sf)
library(rnaturalearth)
library(rnaturalearthhires)

backbone <- read_csv("data/raw/backbone.csv") #taxonomic backbone
country <- ne_countries(scale = "small", returnclass = "sf")
tetra <- read_csv("data/raw/TetraDENSITY_v.1.csv", show_col_types = FALSE) %>% 
  filter(!is.na(Longitude), !is.na(Latitude))

# multiply density to standardize as ind/km2
tetra <- tetra %>% 
  mutate(Density = modify2(Density, Density_unit, function(x, u) {
    fact <- switch (u,
      "ind/ha" = 100,
      "ind/km2" = 1,
      "males/ha" = 100 * 2,
      "pairs/km2" = 2
    )
    return ( x * fact )
  })) %>% 
  transmute(Longitude, Latitude, Density,
            Class, Order, Family,
            Species = paste(Genus, Species))

p <- tetra %>% 
  st_as_sf(coords = c("Longitude", "Latitude"), 
           crs = st_crs(country))

eu <- ne_countries(scale = "large", returnclass = "sf") %>% 
  filter(continent == "Europe") %>% 
  st_crop(xmin = -44, xmax = 66, ymin = 21, ymax = 73)

# retain only points within Europe (as defined by the CHELSAE extent)
p <- p[apply(p %>% st_intersects(eu, sparse = FALSE), 
             MARGIN = 1, FUN = function(x) any(x)), ]

st_write(p, "data/interim/tetradensity.shp", append = FALSE)
