library(tidyverse)
library(rgbif)

# species list
sp <- read_delim("For_the_network_symposium/TetraEU_pairwise_interactions.csv",
                 delim = ";",
                 show_col_types = FALSE) %>% 
  filter(sourceLifestageName == 'young and adults', targetLifestageName == 'adults') %>%  
  dplyr::select(sourceTaxonName, targetTaxonName) %>% 
  pivot_longer(cols = everything()) %>% 
  pull(value) %>% 
  unique()

# get taxonomy
gbif <- lapply(sp, function(x) name_backbone(x)) %>% 
  bind_rows()

# dataframe with taxonomy
sp <- tibble(
  Species = sp,
  Order = gbif[["order"]],
  Family = gbif[["family"]]
)

# dataframe of species IDs
sp_id <- read_delim("For_the_network_symposium/TetraEU_pairwise_interactions.csv",
                    delim = ";",
                    show_col_types = FALSE) %>% 
  filter(sourceLifestageName == 'young and adults', targetLifestageName == 'adults') %>%
  transmute(Species = sourceTaxonName, id = sourceTaxonId) %>% 
  bind_rows(
    read_delim("For_the_network_symposium/TetraEU_pairwise_interactions.csv",
               delim = ";",
               show_col_types = FALSE) %>% 
      filter(sourceLifestageName == 'young and adults', targetLifestageName == 'adults') %>%
      transmute(Species = targetTaxonName, id = targetTaxonId)
  ) %>% 
  distinct_all()

# combine species taxonomy with IDs
sp <- left_join(sp, sp_id)

# write to file
write_csv(sp, "data/raw/species_taxonomy_id.csv")
