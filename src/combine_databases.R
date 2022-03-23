library(rstudioapi)
library(tidyverse)
library(sf)
library(RODBC)
library(mice)

"%out%" <- Negate("%in%")

# load body mass databases ---------------
message("Loading body mass datasets")
# mammals body mass (grams)
mammals <- read_delim("data/raw/elton_mammals", delim = "\t", show_col_types = FALSE) %>% 
  transmute(Species = Scientific,
            Mass = `BodyMass-Value`) #no NA
# birds body mass (grams)
birds <- read_delim("data/raw/elton_birds.txt", delim = "\t", show_col_types = FALSE) %>% 
  transmute(Species = Scientific,
            Mass = `BodyMass-Value`)
# reptiles body mass (grams)
reptiles <- read_csv("data/raw/slavenko.csv", show_col_types = FALSE) %>%
  transmute(Species = Binomial,
            Mass = 10 ^ `Maximum mass (log10(g))`)
# amphibian body mass (grams)
amphibians <- read_csv("data/raw/EuropAmphib.csv", show_col_types = FALSE) %>% 
  transmute(Species = Species,
            Mass = Adult.body.mass,
            Length = Adult.total.length)

backbone <- read_csv("data/raw/backbone.csv", show_col_types = FALSE) #load taxonomic backbone

# add amphibians body mass
backbone <- backbone %>% 
  left_join(amphibians, by = c("Species" = "Species"), na_matches = "never") %>%
  left_join(amphibians, by = c("Accepted" = "Species"), na_matches = "never") %>%
  mutate(Mass = modify2(Mass.x, Mass.y, function(x, y) ifelse(is.na(x), y, x)),
         Length = modify2(Length.x, Length.y, function(x, y) ifelse(is.na(x), y, x))) %>% 
  dplyr::select(-Mass.x, -Mass.y, -Length.x, -Length.y) %>% 
  distinct_all()

# add bird body mass
backbone <- backbone %>% 
  left_join(birds, by = c("Species" = "Species"), na_matches = "never") %>% 
  left_join(birds, by = c("Accepted" = "Species"), na_matches = "never") %>% 
  mutate(Mass = modify2(Mass.x, Mass.y, function(x, y) {
    if (!is.na(x)) {
      x
    } else if (!is.na(y)) {
      y
    } else {
      NA
    }
  })) %>% 
  dplyr::select(-Mass.x, - Mass.y) %>% 
  distinct_all()

# add reptile body mass
backbone <- backbone %>% 
  left_join(reptiles, by = c("Species" = "Species"), na_matches = "never") %>% 
  left_join(reptiles, by = c("Accepted" = "Species"), na_matches = "never") %>%
  mutate(Mass = modify2(Mass.x, Mass.y, function(x, y) {
    if (!is.na(x)) {
      x
    } else if (!is.na(y)) {
      y
    } else {
      NA
    }
  })) %>% 
  select(-Mass.x, - Mass.y)

# add bird and mammal body mass
backbone <- backbone %>% 
  left_join(mammals, by = c("Species" = "Species"), na_matches = "never") %>% 
  left_join(mammals, by = c("Accepted" = "Species"), na_matches = "never") %>%
  mutate(Mass = modify2(Mass.x, Mass.y, function(x, y) {
    if (!is.na(x)) {
      x
    } else if (!is.na(y)) {
      y
    } else {
      NA
    }
  })) %>%
  dplyr::select(-Mass.x, - Mass.y)

missing <- backbone %>%
  filter(is.na(Class)) %>% 
  distinct_all() %>% 
  filter(!is.na(Species))

backbone %>% write_csv("data/interim/backbone-masses.csv")

backbone %>% 
  filter(Species %in% futureweb$Species) %>% 
  mutate(Mass = ifelse(is.na(Mass), FALSE, TRUE)) %>% 
  group_by(Class, Mass) %>% 
  tally() %>% 
  pivot_wider(names_from = Class, values_from = n)

# multiimputation -------------
Y <- backbone %>% 
  filter(!is.na(Order)) %>% 
  mutate(across(where(is.numeric), ~log10(.x))) %>% 
  distinct_all() %>% 
  dplyr::select(Species, Family, Mass) %>% 
  as.data.frame()

mi <- mice(data = Y, m = 10, maxit = 10, method = "norm") #impute only mass
# plot(mi) # for convergence

### overall good imputations
fit <- with(mi, lm(Mass ~ Family))
pooled <- pool(fit$analyses)
summary(pooled)

gamma <- pooled$pooled$fmi
message("Influence of missing data on estimate uncertainty (proportion): ", round(mean(gamma, na.rm = TRUE), 3)) #aim for 0.2 or < 0.5
epsilon <- ( 1 + gamma / mi$m ) ^ - 1
message("Relative efficiency: ", round(mean(epsilon, na.rm = TRUE), 3)) #aim for over 99%

backbone <- backbone %>% 
  left_join(
    complete(mi) %>% 
      as_tibble() %>% 
      mutate(Mass = 10 ^ Mass),
    by = c("Species", "Family")
  ) %>% 
  mutate(`Mass imputed` = ifelse(is.na(Mass.x), TRUE, FALSE),
         Mass = ifelse(is.na(Mass.x), Mass.y, Mass.x)) %>% 
  select(-Mass.x, -Mass.y) %>% 
  relocate(Mass, .before = `Mass imputed`) %>% 
  select(-Length)

backbone %>% write_csv("data/interim/backbone-masses.csv")

sum(backbone$`Mass imputed`) / nrow(backbone)

# load tetra and futureweb --------
# tetradensity
tetra <- st_read("data/interim/tetradensity.shp") %>%
  select(Species, Density) %>% 
  left_join(backbone, by = "Species")

st_write(tetra, "data/interim/tetradensity-masses.shp", append = FALSE) #write p as shapefile

# species of FutureWeb
futureweb <- read_delim("data/raw/futureweb.csv",
                        delim = ";", show_col_types = FALSE) %>%
  select(Species, `Species ID`) %>% 
  left_join(backbone, by = "Species") %>% 
  distinct_all()

futureweb %>% write_csv("data/interim/futureweb-masses.csv")
