library(tidyverse)
library(raster)
library(foreign)
library(ncdf4)
library(igraph)
library("rstudioapi")
setwd(dirname(getActiveDocumentContext()$path))

source("make_density_models.R")

create_metaweb <- function() {
  load("../../For_the_network_symposium/MASTER.bin10000_opthab_tresh0.Rdata") # master = matrix of pixels (in rows) x species (in columns)
  px_sp <- master
  rm(master); gc()
  # metaweb from dataframe to igraph
  df_g.metaweb <- read_delim("../../For_the_network_symposium/TetraEU_pairwise_interactions.csv", delim = ";")
  # check the different lifestages
  df_g.metaweb %>% 
    distinct(sourceLifestageName, targetLifestageName)
  # keep adults only (for sources, that correspond to adults and young, but only adults for targets?)
  metweb.adults <- df_g.metaweb %>% 
    filter(sourceLifestageName == 'young and adults', targetLifestageName == 'adults')
  # create the igraph object
  edgelist <- as.matrix(cbind.data.frame(metweb.adults$sourceTaxonId, metweb.adults$targetTaxonId))
  g.metaweb <- graph_from_edgelist(edgelist, directed = TRUE)
  species.names = V(g.metaweb)$name
}

# load pixel info and mask raster
maskID <- read.dbf("../../For_the_network_symposium/10k_mask2/reference_grid_10km.img.vat.dbf")
maskID <- maskID[order(maskID$Value), ]
m <- rast("../../For_the_network_symposium/10k_mask2/reference_grid_10km.img")

# web <- create_metaweb() #load all necessary files and create metaweb
load("../../For_the_network_symposium/MASTER.bin10000_opthab_tresh0.Rdata") # master = matrix of pixels (in rows) x species (in columns)
px_sp <- master

exDF <- tibble(PAGENAME = rownames(px_sp),
               Value = rep(1, nrow(px_sp)),
               row.names = px_sp[, 1])
               
exDF <- exDF %>% mutate(PxID = as.numeric(as.factor(exDF$PAGENAME)))

PxID <- m #PxID is the raster to link exDF (key is PxID) to other rasters
PxID[!is.na(m)] <- exDF$PxID
names(PxID) <- c("PxID")

env <- c(npp, pcv, temp)
names(env) <- c("NPP", "PCV", "Temp")

env_repr <- 

env <- raster::projectRaster(raster::raster(env),
                             raster::raster(m),
                             method = "ngb")

env <- mask(env, m)
env <- stack(env, PxID)
# plot(env)

# create precipitation coefficient of variation (CV)
env$prepCV <- env[["prec.sd"]] / env[["prec.mean"]]
env = env[[c("NPP", "temp", "prepCV", "PxID")]] #keep only the one you want

if (ncell(m) != ncell(env)) {
  message("Warning: number of cells in environmental variables is not the same as in the mask.")
}
if (!all(xyFromCell(m, 1:ncell(m)) == xyFromCell(env, 1:ncell(env)))) {
  message("Warning: centroid coordiantes do not match")
}

xy <- xyFromCell(env, which(!is.na(values(m)))) #get centroids
# extract environmental values | long pipe
env_vals <- raster::extract(env, xy, cellnumbers = TRUE) %>% #extract
  as_tibble() %>% #just cause I don't want to type 'head()' every time 
  add_column(x = xy[, 1], y = xy[, 2]) %>% #add coordinates
  right_join(exDF) %>% #join with exDF
  arrange(PxID) %>% #sort by PxID
  mutate(CellNumber = cells) %>% #rename 
  dplyr::select(-Value, -cells) %>% #drop columns
  mutate(PAGENAME = as.character(PAGENAME)) #cast as char and remove levels

# env_vals has pixel ID (PxID) as in exDF, PAGENAME as in exDF and CellNumber,
# which is the number of cell in the raster associated with the PxID and
# PAGENAME

write_csv(env_vals, "perpixel_envVar.csv")

### REMOVE THIS CHUNK ##############################################################
# env_vals %>% filter(!is.na(prepCV)) %>% nrow()
# length(na.omit(values(m)))
# nrow(maskID)
# tmp <- m
# tmp[env_vals$CellNumber] <- env_vals$prepCV
# plot(tmp)
# source("create_local_foodwebs.R")
# add all info to the list of species composition and food webs (create_local_foodwebs.R)
# for (x in unique(env_vals$PAGENAME)) {
#   for (vars in c("NPP", "temp", "prepCV", "x", "y", "PxID", "CellNumber")) {
#     g.list[[x]][vars] <- env_vals[[which(env_vals$PAGENAME == x), vars]]
#   }
#   g.list[[x]]["local_sp"] <- sort(g.list[[x]]["local_sp"])
#   masses <- mass %>% 
#     filter(sp_id %in% g.list[[x]]["local_sp"]) %>% 
#     select(sp_id, BM_g)
#   masses <- masses[order(masses[["sp_id"]]), ]
#   g.list[[x]]["mass"] <- mass[["BM_g"]]
# }
# 
# write_rds(g.list, "something.rds")
##################################################

mass <- read_csv("../spp_densitiesv14.csv")

# apply Luca's model to predict abundances

# load the different models
load("../Luca`s model to predict abundance/Supplementary Data 3 - Models/Mammals/M_mod_62.Rdata")
load("../Luca`s model to predict abundance/Supplementary Data 3 - Models/Birds/B_mod_13.Rdata")
load("../Luca`s model to predict abundance/Supplementary Data 3 - Models/Amphibians/A_mod_2.Rdata")
load("../Luca`s model to predict abundance/Supplementary Data 3 - Models/Reptiles/R_mod_16.Rdata")


vert = read.csv("../spp_densitiesv14.csv", na.strings=c(" ", "NA", "NaN", -999))
head(vert)

mammals = vert$sp_id[vert$Class == "Mammalia"]
Amphibia = vert$sp_id[vert$Class == "Amphibia"]
Aves = vert$sp_id[vert$Class == "Aves"]
Reptilia = vert$sp_id[vert$Class == "Reptilia"]

set.abundances = function(pixel){
  # create the dataframe with the local species only
  local.vert = subset(vert, sp_id %in% pixel$local_sp)
  # adding in the information about local conditions
  local.vert$NPP = pixel$NPP
  local.vert$X_Y = paste(x, y, sep = "_")
  local.vert$I.NPP.2. = pixel$NPP^2
  local.vert$Pcv = pixel$NPP$prepCV
  local.vert$abundance = NA
  
  # using the models to determine species' abundances
  local.vert$abundance[local.vert$Class == "Mammalia"] = predict(M_mod_62, subset(local.vert, class == "Mammalia"))
  local.vert$abundance[local.vert$Class == "Amphibia"] = predict(A_mod_2, subset(local.vert, class == "Amphibia"))
  local.vert$abundance[local.vert$Class == "Aves"] = predict(B_mod_13, subset(local.vert, class == "Aves"))
  local.vert$abundance[local.vert$Class == "Reptilia"] = predict(R_mod_16, subset(local.vert, class == "Reptilia"))
  # calculate biomasses
  local.bioms = local.vert$abundance * local.vert$BM_g
  # as I'm not sure that vert and pixel lists of species are in the same order, reordering everything. 
  # (we could check if this last thing is actually needed)
  matches = match(local.vert$sp_id, pixel$local_sp)
  pixel$abundance = local.vert$abundance[matches]
  pixel$biomass = local.bioms[matches]
  return(pixel)
}

ready.to.flux.list = lapply(g.list, set.abundances)
# or parrallel version
# plan(multisession, gc = TRUE, workers = 16)
# ready.to.flux.list = parlapply(g.list, set.abundances, future.chunk.size = 1000)
# plan(sequential) # did not got everything here, but should close the different cluster and free some memory

# calculate the fluxes
library(fluxweb)
get.fluxes = function(pixel){
  # metabolic rates. Here I took the same parameters for every taxa. Can be refine to use
  # rgeression specific to birds, mammals, reptiles, etc. 
  # could also incorporate local temperature
  met = 19.5*pixel$BM^0.71
  # efficiencies: all species are carnivores here, so same efficiency for everybody
  efficiencies = rep(0.85, length(pixel$local_sp))
  # matrix of the network
  mat = as_adjacency_matrix(pixel$fw)
  fluxes = fluxing(mat, biomasses = pixel$biomass, losses = met, efficiencies = efficiencies)
  return(list(x = pixel$x, y = pixel$y, tot.flux = sum(fluxes)))  
}

list.fluxes = lapply(ready.to.flux.list, get.fluxes)
# or to obtain a dataframe directly (and might run faster): (to be tested)
# flux.df = vapply(ready.to.flux.list, get.fluxes, FUN.VALUE = numeric(3))

###### Load the master matrix so we can access the species presence in each pixel
# #### and then associate it to the previous data frame
# 
# load("../MASTER.bin10000_opthab_tresh0.Rdata") # master = matrix of pixels (in rows) x species (in columns)
# px_sp <- master
# rm(master)
# ###
# # # metaweb from dataframe to igraph
# df_g.metaweb = read.csv2("../TetraEU_pairwise_interactions.csv")
# # check the different lifestages
# unique(df_g.metaweb$sourceLifestageName)
# unique(df_g.metaweb$targetLifestageName)
# 
# # keep adults only (for sources, that correspond to adults and young, but only adults for targets?)
# metweb.adults = subset(df_g.metaweb, sourceLifestageName == 'young and adults' & targetLifestageName == 'adults')
# # create the igraph object
# edgelist = as.matrix(cbind.data.frame(metweb.adults$sourceTaxonId, metweb.adults$targetTaxonId))
# g.metaweb = graph_from_edgelist(edgelist, directed = TRUE)
# 
# species.names = V(g.metaweb)$name
# 
# 
# # make that the first row contains values, and rows have now row names: 
# r.names = px_sp[,1]
# px_sp = px_sp[,-1]
# row.names(px_sp) = r.names
# sp <- as_tibble(px_sp) %>% 
#   add_column(PageName = rownames(px_sp)) %>% 
#   pivot_longer(cols = 1:(ncol(.) - 1),
#                names_to = "species",
#                values_to = "suitability")
# sp <- sp %>% filter(suitability > 0)
# total.ds = left_join(sp, env_vals)
# 
# write_csv(total.ds, "Species_perPixel_var.csv")
