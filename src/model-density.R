library(tidyverse)
library(rgbif)

# species list
sp <- read_csv("data/raw/species_taxonomy_id.csv", show_col_types = FALSE)

load("For_the_network_symposium//MASTER.bin10000_opthab_tresh0.Rdata") # master = matrix of pixels (in rows) x species (in columns)
master[1:3, 1:3]

env <- read_csv("For_the_network_symposium/perpixel_envVar.csv",
                show_col_types = FALSE) %>% 
  filter(!is.na(temp)) #remove pixels without temperature

cellname <- as.character(master[which(master[, sp$id[1]] == 1), "PAGENAME"])
