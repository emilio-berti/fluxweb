library(tidyverse)
library(terra)

plot_distr <- function(x) {
  if (!"px_sp" %in% ls()) {
    path <- "~/idiv-mount/groups/tib/FutureWeb/For_the_network_symposium/"
    load(file.path(path, "MASTER.bin10000_opthab_tresh0.Rdata")) # master = matrix of pixels (in rows) x species (in columns)
    px_sp <- master
  }
  if (!"vert" %in% ls()) {
    # load species table
    vert <- read.csv(file.path(path, "spp_densitiesv14.csv"),
                     na.strings = c(" ", "NA", "NaN", -999)) %>%
      as_tibble() %>%
      transmute(sp_id, BM_g, Class, Order, Family, Species) %>%
      filter(!is.na(Class), !is.na(BM_g)) %>% #remove species without class and body mass information
      mutate(Species = gsub("_", " ", Species))

    # reassign taxonomic names ---------------
    vert <- read_delim(file.path(path, "Original data - FutureWeb/species_codes_and_taxonomy.csv"),
                       delim = ";", show_col_types = FALSE) %>%
      dplyr::select(`Species ID`, Species)

    tax <- read_csv("data/interim/backbone-masses.csv", show_col_types = FALSE)

    vert <- vert %>% left_join(tax, by = "Species") %>%
      dplyr::select(-`Mass imputed`, -Species)
  }

  id <- vert %>%
    filter(Accepted == x) %>%
    pull(`Species ID`)
  occ <- px_sp[, c("PAGENAME", id)] %>%
    as_tibble()

  m <- rast(file.path(path, "/10k_mask2/reference_grid_10km.img"))

  pixels = m %>%
    extract(xyFromCell(m, 1:ncell(m))) %>%
    as_tibble()

  values(m) <- left_join(pixels, occ, by = c("PageName" = "PAGENAME"))[[2]]

  plot(m, col = c("grey20", "green4"), main = x, frame = FALSE, axes = FALSE)

  return(m)

}

fox <- plot_distr("Vulpes vulpes")
bubo <- plot_distr("Bubo bubo")

plot(fox * bubo)

plot_distr("Cervus elaphus")
