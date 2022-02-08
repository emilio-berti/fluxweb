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
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#create lm for densities and load env rasters reprojected
source("make_density_models.R") #source it here as it has rm() at the end

boltz <- 0.00008617343

create_metaweb <- function() {
  load("../../For_the_network_symposium/MASTER.bin10000_opthab_tresh0.Rdata") # master = matrix of pixels (in rows) x species (in columns)
  px_sp <- master
  rm(master); gc()
  # metaweb from dataframe to igraph
  df_g.metaweb <- read_delim("../../For_the_network_symposium/TetraEU_pairwise_interactions.csv",
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

pixel_flux <- function(pagename, suit_th = 0, species_list = c()) {
  # table of environmental variable
  sub_env <- env %>% filter(PAGENAME == pagename)
  # table of species per pixel
  sp_sub <- px_sp %>% 
    filter(PAGENAME == pagename) %>% 
    dplyr::select(-PAGENAME)
  sp_sub <- colnames(sp_sub)[which(sp_sub > suit_th)] #species with suit > suit_th
  
  # table of vertebrate 'traits'
  sub_vert <- vert %>% 
    filter(sp_id %in% sp_sub) %>% 
    add_column(NPP = sub_env$NPP, Pcv = sub_env$prepCV,
               X_Y = paste(sub_env$x, sub_env$y, sep = "_")) %>% 
    distinct_all()
  if (nrow(sub_vert) < 2) {
    stop("Only one species in this pixel")
  }
  # model densities THIS NEEDS TO BE FIXED
  #I think density and body mass (for some taxa) in Luca's model is log10-transformed
  sub_vert <- sub_vert %>% 
    mutate(BM_g = log10(BM_g),
           NPP = log10(NPP * 80)) %>% #correction for area size in LAEA vs longlat (see luca_to_mask_area_size.R)
    mutate(across(where(is.factor), as.character)) %>% #cast as character
    mutate(density = NA)
  # Table S1
  # sub_vert[sub_vert$Class == "Amphibia", "density"] <- sub_vert %>% 
  #   filter(Class == "Amphibia") %>% 
  #   mutate(density = 2.6719991 + 0.1481797 * BM_g - 0.3279702 * BM_g) %>% 
  #   pull(density)
  # sub_vert[sub_vert$Class == "Aves", "density"] <- sub_vert %>% 
  #   filter(Class == "Aves") %>% 
  #   mutate(density = -106.2156521 - 0.3118665 * BM_g + 19.0011294 * sub_env$NPP) %>% 
  #   pull(density)
  # sub_vert[sub_vert$Class == "Mammalia", "density"] <- sub_vert %>% 
  #   filter(Class == "Mammalia") %>% 
  #   mutate(density = -57.907068317 - 0.328039394 * BM_g + 10.954635813 * sub_env$NPP - 0.004343108 * sub_env$prepCV) %>% 
  #   pull(density)
  # sub_vert[sub_vert$Class == "Reptilia", "density"] <- sub_vert %>% 
  #   filter(Class == "Reptilia") %>% 
  #   mutate(density = 2.452848017 - 0.379658597 * BM_g - 0.009144478 * sub_env$prepCV) %>% 
  #   pull(density)
  sub_vert <- sub_vert %>% 
    mutate(density = 2.2 + 2.3805962 - 0.5262883 * BM_g)
  # sub_vert %>% 
  #   filter(density < 0) %>% 
  #   arrange(density)
  # sub_vert %>% 
  #   ggplot() +
  #   aes(BM_g, density, col = Class) +
  #   geom_point()
  # sometimes we have negative densities as per Luca's model. Fix this later, for now just abs()
  if (any(sub_vert$density < 0)) {
    warning("Pixel ", pagename, ": negative densities found")
    sub_vert$density <- abs(sub_vert$density)
  }
  
  # calculate biomass
  sub_vert <- sub_vert %>% 
    mutate(BM_g = 10 ^ BM_g, density = 10 ^ density) %>% 
    mutate(sp_id = as.character(sp_id),
           biomass = density * BM_g)
  # convert graph to adjancency matrix (simplify and remove isolated vertices)
  sub_g <- induced_subgraph(g.metaweb, sub_vert$sp_id)
  isolated <- names(which(degree(sub_g, mode = "total") == 0)) #remove isolated nodes
  if (length(isolated) > 0) {
    warning("Pixel ", pagename, ": isolated vertices found - deleting them")
    sub_g <- sub_g %>% igraph::delete_vertices(isolated) %>% 
      as_adjacency_matrix(sparse = FALSE)
  } else {
    sub_g <- sub_g %>% as_adjacency_matrix(sparse = FALSE)
  }
  sub_vert <- sub_vert %>% 
    filter(sp_id %in% colnames(sub_g))
  
  # fluxing
  temp <- sub_env$temp
  met.rate <- exp((0.71 * log(sub_vert$BM_g) + 19.50) - 0.69/(boltz*(273.15+temp))) # for endotherms
  met.rate[sub_vert$Class == "Amphibia"] <- exp((0.71 * log(sub_vert$BM_g[sub_vert$Class == "Amphibia"]) + 18.05) - 0.69/(boltz*(273.15+temp))) # for amphibians
  met.rate[sub_vert$Class == "Reptilia"] <- exp((0.71 * log(sub_vert$BM_g[sub_vert$Class == "Reptilia"]) + 18.02) - 0.69/(boltz*(273.15+temp))) # for reptiles
  
  # note: memoise fluxing will not help, as temperature is unique for each 
  # pixel, e.g. max(table(env$temp)) == 1 is TRUE
  fluxes <- fluxing(mat = sub_g,
                    biomasses = sub_vert$biomass,
                    losses = met.rate,
                    efficiencies = rep(0.906, length(met.rate)))
  
  focus_species <- which(sub_vert$Species %in% species_list)
  
  if (length(focus_species) > 0) {
    avg_out <- mean(rowSums(fluxes, na.rm = TRUE)[focus_species] / sub_vert$biomass[focus_species])
    avg_in <- mean(colSums(fluxes, na.rm = TRUE)[focus_species] / sub_vert$biomass[focus_species])
  } else {
    avg_out <- NA
    avg_in <- NA
  }

  output <- tibble(PAGENAME = sub_env$PAGENAME,
                   Sum_of_fluxes = sum(fluxes),
                   Number_of_species = length(sub_vert$biomass),
                   Biomasses = list(sub_vert$biomass),
                   Avg_out_fluxes = avg_out,
                   Avg_in_fluxes = avg_in)
  return (output)
}

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
vert <- read.csv("../../For_the_network_symposium/spp_densitiesv14.csv", na.strings = c(" ", "NA", "NaN", -999)) %>% 
  as_tibble() %>% 
  transmute(sp_id, BM_g, Class, Order, Family, Species) %>% 
  filter(!is.na(Class), !is.na(BM_g)) %>% #remove species without class and body mass information
  filter(sp_id %in% adults) %>% #keep only species recorded as adults of young
  mutate(Species = gsub("_", " ", Species))

# reassign taxonomic names
vert <- read_delim("../../For_the_network_symposium/Original data - FutureWeb/species_codes_and_taxonomy.csv",
                        delim = ";", show_col_types = FALSE) %>%
  dplyr::select(`Species ID`, Species)

tax <- read_csv("../data/interim/backbone-masses.csv", show_col_types = FALSE)

vert <- vert %>% left_join(tax, by = "Species") %>% 
  dplyr::select(-`Mass imputed`, -Species)

exDF <- tibble(PAGENAME = rownames(px_sp),
               Value = rep(1, nrow(px_sp)),
               PxID = as.numeric(as.factor(PAGENAME)))

# parallel comps -----------
# I append each pixel output to a ../pixel-fluxes.csv
# so intermediate results are saved and not needed to compute them
# again.
# For ~ 100,000 pixels, this will take ~ 65 hours on a single core.
done <- tryCatch(read_csv("~/pixel-fluxes.csv", col_names = FALSE) %>% 
                   pull(X1) %>% 
                   unique(),
                 error = function(e) NULL)
to_do <- setdiff(env$PAGENAME, done)

RhpcBLASctl::blas_set_num_threads(1)
RhpcBLASctl::omp_set_num_threads(1)

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
