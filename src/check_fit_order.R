library(tidyverse)
library(sf)
library(terra)

source("make_density_models.R")

futureweb <- read_csv("data/interim/futureweb-masses.csv", show_col_types = FALSE)

modelled <- sapply(lm_density, function(x) unique(x$model$Order)) %>% 
  unlist() %>% 
  as.vector()

futureweb %>% filter(!Order %in% modelled)

for (cl in names(lm_density)) {
  m <- lm_density[[cl]]
  if (cl == names(lm_density)[1]) {
    lm_data <- m$model %>% 
      as_tibble() %>%
      mutate(fitted = m$fitted.values,
             Class = cl)
  } else {
    lm_data <- lm_data %>% bind_rows(m$model %>% 
                                       as_tibble() %>%
                                       mutate(fitted = m$fitted.values,
                                              Class = cl))
  }
}

lm_data %>% 
  pivot_longer(cols = c(`log10(Density)`, fitted)) %>% 
  mutate(name = ifelse(name == "fitted", "Fitted", "Recorded")) %>%
  ggplot() +
  geom_boxplot(aes(Order, value, col = name)) +
  ggsci::scale_color_nejm() +
  theme_light() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))

missing <- futureweb %>% 
  filter(!Order %in% lm_data$Order) %>% 
  select(`Species ID`, Class, Order) %>% 
  distinct(Class, Order)

read_delim("data/raw/futurweb-interactions.csv", delim = ";", show_col_types = FALSE) %>% 
  filter(sourceLifestageName == 'young and adults', targetLifestageName == 'adults') %>% 
  filter(
    sourceTaxonId %in% (futureweb %>% 
                          filter(Order %in% missing$Order) %>% 
                          pull(`Species ID`)) |
      targetTaxonId %in% (futureweb %>% 
                            filter(Order %in% missing$Order) %>% 
                            pull(`Species ID`))
  ) %>% 
  select(sourceTaxonName, targetTaxonName) %>% 
  distinct_all() %>% 
  left_join(backbone %>% select(Species, Class, Order), by = c("sourceTaxonName" = "Species")) %>%
  left_join(backbone %>% select(Species, Class, Order), by = c("targetTaxonName" = "Species"),
            suffix = c(".source", ".target")) %>% 
  select(starts_with("Class") | starts_with("Order")) %>% 
  relocate(Order.source, .before = Class.target) %>% 
  distinct_all()


(read_delim("data/raw/futurweb-interactions.csv", delim = ";", show_col_types = FALSE) %>% 
  filter(sourceLifestageName == 'young and adults', targetLifestageName == 'adults') %>% 
  filter(
    sourceTaxonId %in% (futureweb %>% 
                          filter(Order %in% missing$Order) %>% 
                          pull(`Species ID`)) |
      targetTaxonId %in% (futureweb %>% 
                            filter(Order %in% missing$Order) %>% 
                            pull(`Species ID`))
  ) %>% nrow()) / 
  read_delim("data/raw/futurweb-interactions.csv", delim = ";", show_col_types = FALSE) %>% nrow()


read_delim("data/raw/futurweb-interactions.csv", delim = ";", show_col_types = FALSE) %>% 
  filter(sourceLifestageName == 'young and adults', targetLifestageName == 'adults') %>% 
  mutate(missing = ifelse (sourceTaxonId %in% (futureweb %>% 
                                 filter(Order %in% missing$Order) %>% 
                                 pull(`Species ID`)) | 
             targetTaxonId %in% (futureweb %>% 
                                 filter(Order %in% missing$Order) %>% 
                                 pull(`Species ID`)), "missing", "complete")) %>% 
  left_join(backbone %>% select(Species, Class, Order), by = c("sourceTaxonName" = "Species")) %>%
  left_join(backbone %>% select(Species, Class, Order), by = c("targetTaxonName" = "Species"),
            suffix = c(".source", ".target")) %>%
  group_by(missing, Class.source, Class.target) %>% 
  tally() %>% 
  pivot_wider(names_from = missing, values_from = n) %>% 
  replace_na(list(complete = 0, missing = 0)) %>% 
  arrange(desc(missing)) %>% 
  mutate(`complete fraction` = complete / (complete + missing)) %>% 
  write_csv("docs/proportion-of-interactions-modelled.csv")

  