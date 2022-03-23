library(tidyverse)
library(sf)
library(terra)

m <- lm_density$Mammalia
d <- m$model %>% as_tibble()

d <- d %>%
  mutate(predMass = `log10(Mass)` * m$coefficients["log10(Mass)"] +
           `log10(Mass)` ^ 2 * m$coefficients["I(log10(Mass)^2)"]) %>%
  mutate(predPCV = `PCV` * m$coefficients["PCV"] +
           `PCV` ^ 2 * m$coefficients["I(PCV^2)"]) %>%
  mutate(predNPP = `log10(NPP)` * m$coefficients["log10(NPP)"] +
           `log10(NPP)` ^ 2 * m$coefficients["I(log10(NPP)^2)"]) %>%
  pivot_longer(cols = predMass : predNPP, names_to = "fit") %>%
  pivot_longer(cols = c(2, 4, 6), names_to = "covariate", values_to = "x")

d <- d %>%
  mutate(fit = gsub("pred", "", fit),
         covariate = gsub("[log10(|)]", "", covariate)) %>%
  filter(covariate == fit)

d %>%
  ggplot() +
  geom_line(aes(x, value)) +
  geom_rug(aes(x)) +
  facet_wrap(~covariate, scales = "free") +
  theme_bw()
