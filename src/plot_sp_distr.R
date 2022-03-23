library(tidyverse)
library(terra)
library(progress)

get_range_area <- function() {
  vert %>%
    transmute(Species = `Species ID`,
              Accepted) %>%
    left_join(
      tibble(Species = colnames(px_sp[, -1]),
             Area = colSums(px_sp[, -1], na.rm = TRUE)),
      by = "Species"
    )
}

plot_distr <- function(x) {
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
}

# load data ---------------
path <- "~/idiv-mount/groups/tib/FutureWeb/For_the_network_symposium/"
load(file.path(path, "MASTER.bin10000_opthab_tresh0.Rdata")) # master = matrix of pixels (in rows) x species (in columns)
px_sp <- master
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

# plot ------------
ranges <- get_range_area()
ranges <- ranges %>% arrange(desc(Area))

S <- 500
pb <- progress_bar$new(total = S)
cairo_pdf("docs/plots/species-ranges.pdf", width = 6, height = 4, onefile = TRUE)
for (x in ranges$Accepted[1:S]) {
  pb$tick()
  plot_distr(x)
}
dev.off()

