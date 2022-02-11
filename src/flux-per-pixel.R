library(tidyverse)
library(sf)
library(terra)
library(foreign)
library(ncdf4)
library(igraph)
library(fluxweb)
library(foreach)
library(doParallel)
library(lme4)

#create lm for densities and load env rasters reprojected
source("src/make_density_models.R") #source it here as it has rm() at the end

path <- "/data/tib-group-share/FutureWeb/For_the_network_symposium/"

boltz <- 0.00008617343

create_metaweb <- function() {
  load(file.path(path, "MASTER.bin10000_opthab_tresh0.Rdata")) # master = matrix of pixels (in rows) x species (in columns)
  px_sp <- master
  rm(master); gc()
  # metaweb from dataframe to igraph
  df_g.metaweb <- read_delim(file.path(path, "TetraEU_pairwise_interactions.csv"),
                             delim = ";",
                             show_col_types = FALSE)
  # keep adults only (for sources, that correspond to adults and young, but only adults for targets?)
  metweb.adults <- df_g.metaweb %>%
    filter(sourceLifestageName == 'young and adults', targetLifestageName == 'adults')
  adults <- metweb.adults %>%
    dplyr::select(sourceTaxonId, targetTaxonId) %>%
    pivot_longer(cols = everything()) %>%
    pull(value) %>%
    unique()
  # create the igraph object
  edgelist <- as.matrix(cbind.data.frame(metweb.adults$sourceTaxonId, metweb.adults$targetTaxonId))
  g.metaweb <- graph_from_edgelist(edgelist, directed = TRUE)
  species.names = V(g.metaweb)$name
  return(list(metaweb = g.metaweb, px_sp = px_sp, adults = adults))
}

pixel_flux <- function(pagename, temp, suit_th = 0, species_list = c()) {
  # table of environmental variable
  sub_env <- env %>% filter(PAGENAME == pagename)
  # table of species per pixel
  sp_sub <- px_sp %>%
    filter(PAGENAME == pagename) %>%
    dplyr::select(-PAGENAME)
  sp_sub <- colnames(sp_sub)[which(sp_sub > suit_th)] #species with suit > suit_th

  # table of vertebrate 'traits'
  sub_vert <- vert %>%
    filter(`Species ID` %in% sp_sub) %>%
    distinct_all()
  if (nrow(sub_vert) < 2) {
    stop("Only one species in this pixel")
  }
  # calculate biomass
  biomasses <- dens %>%
    filter(`Species ID` %in% sub_vert$`Species ID`) %>%
    distinct_all() %>%
    left_join( #add effect of climate
      env %>%
        slice_sample(n = 1) %>%
        pivot_longer(cols = amphibia:reptilia) %>%
        transmute(Class = gsub("^(\\w)", "\\U\\1\\L\\2", name, perl = TRUE),
                  Env.effect = value),
      by = "Class"
    ) %>%
    mutate(Density = Density + Env.effect) %>%
    transmute(`Species ID`, Accepted, Class,
              Mass,
              Density = 10 ^ Density, #back-transform densities
              Biomass = Mass * Density) #biomasses


  # convert graph to adjancency matrix (simplify and remove isolated vertices)
  sub_g <- induced_subgraph(g.metaweb, biomasses$`Species ID`)
  isolated <- names(which(degree(sub_g, mode = "total") == 0)) #remove isolated nodes
  if (length(isolated) > 0) {
    warning("Pixel ", pagename, ": isolated vertices found - deleting them")
    sub_g <- sub_g %>% igraph::delete_vertices(isolated) %>%
      as_adjacency_matrix(sparse = FALSE)
  } else {
    sub_g <- sub_g %>% as_adjacency_matrix(sparse = FALSE)
  }
  biomasses <- biomasses %>%
    filter(`Species ID` %in% colnames(sub_g))

  # fluxing --------
  biomasses <- biomasses %>%
    mutate(`metabolic rate` = modify2(Class, Mass, function(x, y) {
      switch (which(c("Aves", "Mammalia", "Amphibia", "Reptilia") == x),
              exp((0.71 * log10(y) + 19.50) - 0.69 / (boltz * (273.15 + temp))),
              exp((0.71 * log10(y) + 19.50) - 0.69 / (boltz * (273.15 + temp))),
              exp((0.71 * log10(y) + 18.05) - 0.69 / (boltz * (273.15 + temp))),
              exp((0.71 * log10(y) + 18.05) - 0.69 / (boltz * (273.15 + temp)))
      )
    }) %>% as.numeric())
  # met.rate <- exp((0.71 * log10(biomasses$Mass) + 19.50) - 0.69 / (boltz * (273.15 + temp))) # for endotherms
  # met.rate[biomasses$Class == "Amphibia"] <- exp((0.71 * log10(biomasses$Mass[biomasses$Class == "Amphibia"]) + 18.05) - 0.69/(boltz*(273.15+temp))) # for amphibians
  # met.rate[biomasses$Class == "Reptilia"] <- exp((0.71 * log10(biomasses$Mass[biomasses$Class == "Reptilia"]) + 18.02) - 0.69/(boltz*(273.15+temp))) # for reptiles

  # pixel, e.g. max(table(env$temp)) == 1 is TRUE
  fluxes <- fluxing(mat = sub_g,
                    biomasses = biomasses$Biomass,
                    losses = biomasses$`metabolic rate`,
                    efficiencies = rep(0.906, nrow(biomasses)))

  focus_species <- which(biomasses$Accepted %in% species_list)
  if (length(focus_species) > 0) {
    focus_avg_out <- mean(rowSums(fluxes, na.rm = TRUE)[focus_species] / biomasses$Biomass[focus_species])
    focus_avg_in <- mean(colSums(fluxes, na.rm = TRUE)[focus_species] / biomasses$Biomass[focus_species])
  } else {
    focus_avg_out <- NA
    focus_avg_in <- NA
  }
  fl <- fluxes
  fl[fl == 0] <- NA
  fl <- log10(fl)
  avg_out <- mean( rowSums(fl, na.rm = TRUE) / biomasses$Biomass )
  avg_in <- mean( colSums(fl, na.rm = TRUE) / biomasses$Biomass )

  output <- tibble(PAGENAME = sub_env$PAGENAME,
                   Temperature = temp,
                   Sum_of_fluxes = sum(fluxes),
                   Number_of_species = nrow(biomasses),
                   Biomasses = list(biomasses$Biomass),
                   Avg_out_fluxes = avg_out,
                   Avg_in_fluxes = avg_in,
                   Focus_avg_out_fluxes = focus_avg_out,
                   Focus_avg_in_fluxes = focus_avg_in)
  return (output)
}

# raster to get temperatures for fluxweb
temperature <- rast("data/interim/temperature-reprojected.tif")

tmp_metaweb <- create_metaweb() #load all necessary files and create metaweb
# metaweb graph, but reducing canniabalism by 0.1
g.metaweb <- tmp_metaweb[[1]]
adj <- g.metaweb %>% as_adjacency_matrix()
diag(adj) <- diag(adj) * 0.1
g.metaweb <- graph_from_adjacency_matrix(adj)

px_sp <- tmp_metaweb[[2]]
adults <- tmp_metaweb[[3]]
rm(tmp_metaweb); gc()

# load species table
vert <- read.csv(file.path(path, "spp_densitiesv14.csv"),
                 na.strings = c(" ", "NA", "NaN", -999)) %>%
  as_tibble() %>%
  transmute(sp_id, BM_g, Class, Order, Family, Species) %>%
  filter(!is.na(Class), !is.na(BM_g)) %>% #remove species without class and body mass information
  filter(sp_id %in% adults) %>% #keep only species recorded as adults of young
  mutate(Species = gsub("_", " ", Species))

# reassign taxonomic names ---------------
vert <- read_delim(file.path(path, "Original data - FutureWeb/species_codes_and_taxonomy.csv"),
                   delim = ";", show_col_types = FALSE) %>%
  dplyr::select(`Species ID`, Species)

tax <- read_csv("data/interim/backbone-masses.csv", show_col_types = FALSE)

vert <- vert %>% left_join(tax, by = "Species") %>%
  dplyr::select(-`Mass imputed`, -Species)

modelled_orders <- sapply(1:4, function(i) {
  n <- names(lm_density[[i]])
  n <- n[grepl("Order", n)]
  return(n)
}) %>%
  unlist()

# the following is tidyverse magic - with minor headaches
dens <- vert %>%
  filter(!is.na(Order)) %>% #drop one species without Order
  mutate(Density = pmap( #this is basically a dplyr sapply
    list(Mass, Class, Order), #list of argument for anonymous function
    function(mass, cl, or) {
      if (or %in% modelled_orders) { #is the order one of those modeled
        log10(mass) * lm_density[[cl]][["log10(Mass)"]] +
          log10(mass) ^ 2 * lm_density[[cl]][["I(log10(Mass)^2)"]] +
          lm_density[[cl]][ #this finds the correct Order coefficient
            grepl(or, names(lm_density[[cl]]))
          ]
      } else { #if Order was not modeled
        log10(mass) * lm_density[[cl]][["log10(Mass)"]] +
          log10(mass) ^ 2 * lm_density[[cl]][["I(log10(Mass)^2)"]] +
          mean(lm_density[[cl]][ #average of all Order coefficient
            grepl("Order", names(lm_density[[cl]]))
          ])
      }
    }
  ) %>% unlist())

# env contains the climatic effect on densities. E.g.:
set.seed(124124)
dens %>%
  slice_sample(prop = .1) %>%
  left_join(
    env %>%
      slice_sample(n = 1) %>%
      pivot_longer(cols = amphibia:reptilia) %>%
      transmute(Class = gsub("^(\\w)", "\\U\\1\\L\\2", name, perl = TRUE),
                Env.effect = value),
    by = "Class"
  ) %>%
  mutate(Density = Density + Env.effect) %>%
  drop_na() %>%
  ggplot() +
  geom_boxplot(aes(Class, Density))

exDF <- tibble(PAGENAME = rownames(px_sp),
               Value = rep(1, nrow(px_sp)),
               PxID = as.numeric(as.factor(PAGENAME)))

# for fluxweb:
# pass pixel pagename and temperature.

# parallel comps -----------
# I append each pixel output to a data/final/pixel-fluxes.csv
# so intermediate results are saved and not needed to compute them
# again.
done <- tryCatch(read_csv("data/final/pixel-fluxes.csv", col_names = FALSE) %>%
                   pull(X1) %>%
                   unique(),
                 error = function(e) NULL)
to_do <- setdiff(env$PAGENAME, done)

RhpcBLASctl::blas_set_num_threads(1)
RhpcBLASctl::omp_set_num_threads(1)

# no focus species ------------------
registerDoParallel(cores = 6)
output_file <- "data/final/pixel-fluxes.csv"
ans <- foreach(pixel = to_do,
               .combine = "rbind",
               .errorhandling = "remove",
               .verbose = FALSE) %dopar% {
                 temp <- terra::extract(temperature,
                                        env %>%
                                          filter(PAGENAME == pixel) %>%
                                          dplyr::select(x, y))[["tas"]]
                 if (is.na(temp)) {
                   stop("No temperature recorded for this pixel")
                 }
                 res <- pixel_flux(pixel, temp)
                 write_csv(res,
                           output_file,
                           append = TRUE)
                 # write_csv(unnest(res[, c("PAGENAME", "Biomasses")], cols = c(Biomasses)),
                 #           "../pixel-biomasses.csv",
                 #           append = TRUE)
                 return(res)
               }
message(length(to_do) - nrow(ans), " of ", length(to_do), " failed.")
stopImplicitCluster()



source("controllers.R") #this load the species list of controllers of:
# - rodents: regulate_rodents
# - invertebrates: regulate_inv
# This list can be passed to pixel_flux,
# e.g. pixel_flux(pixel, species_list = regulate_inv)

# for (pixel in to_do) {
#   if (which(to_do == pixel) %% 10 == 0)
#     message("Pixel ", which(to_do == pixel), " of ", length(to_do))
#   res <- tryCatch(pixel_flux(pixel),
#                   error = function(e) NULL)
#   if (!is.null(res)) {
#     write_csv(res, "~/pixel-fluxes.csv", append = TRUE)
#   }
# }

output_file <- "../pixel-fluxes_vole.csv"
registerDoParallel(cores = 6)
ans <- foreach(pixel = to_do,
               .combine = "rbind",
               .errorhandling = "remove",
               .verbose = FALSE) %dopar% {
                 res <- pixel_flux(pixel, species_list = "Microtus arvalis")
                 write_csv(res[, c("PAGENAME", "Sum_of_fluxes", "Number_of_species",
                                   "Avg_out_fluxes", "Avg_in_fluxes")],
                           output_file,
                           append = TRUE)
                 # write_csv(unnest(res[, c("PAGENAME", "Biomasses")], cols = c(Biomasses)),
                 #           "../pixel-biomasses.csv",
                 #           append = TRUE)
                 return(res)
               }
message(length(to_do) - nrow(ans), " of ", length(to_do), " failed.")
stopImplicitCluster()
