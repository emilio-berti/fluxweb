library(tidyverse)
library(sf) #to read shapefile of tetradensity
library(testthat) #for testing the steps
library(rgbif)
library(foreach)
library(doParallel)

# load all lists of species names ----------
# European amphibians
amphibians <- read_csv("data/raw/EuropAmphib.csv", show_col_types = FALSE) %>% 
  select(Species) %>% #there's no order or family information
  distinct_all() %>% 
  filter(!is.na(Species))
## for birds and mammal `BodyMass-SpecLevel` == 0 -> inferred from genus and family.
# for birds, also BodyMass-Comment == MarceloR -> inferred, but marked as SpecLevel == 1.
# birds
birds <- read_delim("data/raw/elton_birds.txt", delim = "\t", show_col_types = FALSE) %>% 
  transmute(Species = Scientific, Family = BLFamilyLatin, Order = IOCOrder) %>% 
  distinct_all()
# mammals
mammals <- read_delim("data/raw/elton_mammals", delim = "\t", show_col_types = FALSE) %>% 
  transmute(Species = Scientific, Family = MSWFamilyLatin) %>% 
  distinct_all()
# reptiles
reptiles <- read_csv("data/raw/slavenko.csv", show_col_types = FALSE) %>% 
  transmute(Species = Binomial, Family) %>% #no order or family information
  distinct_all()
# tetradensity
tetra <- st_read("data/interim/tetradensity.shp") %>% 
  st_drop_geometry() %>% 
  as_tibble() %>% 
  select(Species, Class, Order, Family) %>% 
  distinct_all() %>% 
  mutate(across(everything(), ~as.character(.x)))
# species of FutureWeb
futureweb <- read_delim("For_the_network_symposium/Original data - FutureWeb/species_codes_and_taxonomy.csv",
                        delim = ";", show_col_types = FALSE) %>%
  transmute(Species, Class, Order, Family) %>% 
  distinct_all()

# harmonize --------------
# check if any species in one trait list is present in another trait list.
# If TRUE, this can potentially confuse the harmonization process
test_that("No name is present in two trait lists",
          expect_equal(
            any(amphibians$Species %in% mammals$Species,
                amphibians$Species %in% birds$Species,
                amphibians$Species %in%reptiles$Species,
                birds$Species %in% mammals$Species,
                birds$Species %in% reptiles$Species,
                mammals$Species %in% reptiles$Species
            ), FALSE))

sp <- bind_rows(
  tetra,
  amphibians %>% mutate(Class = "Amphibia"),
  birds %>% mutate(Class = "Aves"),
  mammals %>% mutate(Class = "Mammalia"),
  reptiles %>% mutate(Class = "Reptilia"),
  futureweb
) %>%
  select(Species, Class) %>% 
  distinct_all()

registerDoParallel(6)

gbif <- foreach(x = sp$Species, .combine = "rbind") %dopar%
  {
    first_search <- name_backbone(x)
    if (!"status" %in% names(first_search)) {
      ans <- tibble(Species = x,
                    GBIF = NA,
                    Class = NA,
                    Order = NA,
                    Family = NA)
      return(ans)
    }
    if (first_search$status == "ACCEPTED") {
      ans <- tibble(Species = x,
                    GBIF = first_search$canonicalName,
                    Class = first_search$class,
                    Order = first_search$order,
                    Family = first_search$family)
      
    } else {
      key <- name_backbone(x)$speciesKey
      second_search <- name_usage(key)$data
      if (second_search$taxonomicStatus == "ACCEPTED") {
        ans <- tibble(Species = x,
                      GBIF = second_search$canonicalName,
                      Class = second_search$class,
                      Order = second_search$order,
                      Family = second_search$family)
      } else {
        ans <- tibble(Species = x,
                      GBIF = NA,
                      Class = NA,
                      Order = NA,
                      Family = NA)
      }
    }
    return(ans)
  }

stopImplicitCluster()

# check names not found
gbif %>% 
  filter(is.na(Class)) %>% 
  pull(Species)

gbif %>% 
  filter(Species != GBIF) %>% 
  distinct_all()

test_that("GBIF and FutureWeb species names are identical",
          {
            d <- futureweb %>% 
              left_join(gbif, by = "Species") %>% 
              distinct_all()
            expect_identical(d$Species, d$GBIF)
          })
# NOTE: I think that futureweb was only partially harmonized using ITIS,
# as 'Bufotes belearica' is not found in ITIS, but 'Bufotes belaricus', as 
# returned by GBIF, is. I follow GBIF from here on.

test_that("GBIF and FutureWeb class names are identical",
          {
            d <- futureweb %>% 
              left_join(gbif, by = "Species") %>% 
              distinct_all()
            expect_identical(d$Class.x, d$Class.y)
          })

test_that("GBIF and FutureWeb order names are identical",
          {
            d <- futureweb %>% 
              left_join(gbif, by = "Species") %>% 
              distinct_all()
            expect_identical(d$Order.x, d$Order.y)
          })

test_that("GBIF and FutureWeb family names are identical",
          {
            d <- futureweb %>% 
              left_join(gbif, by = "Species") %>% 
              distinct_all()
            expect_identical(d$Family.x, d$Family.y)
          })

# Overall, there are minimal differences and I trust GBIF results more than the 
# original list of species names. I follow GBIF backbone from now on.
gbif %>% 
  transmute(Species, Accepted = GBIF,
            Class, Order, Family) %>% 
  write_csv("data/raw/backbone.csv")
