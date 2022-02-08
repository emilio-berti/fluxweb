library(tidyverse)
library(sf)
library(terra)
library(easystats)
library(cowplot)
library(MuMIn)

options(na.action = "na.fail")

# model selection amphibians --------------
sub_d <- d %>%
  filter(Class == "Amphibia") %>% 
  as.data.frame()
m <- lm(log10(Density) ~ log10(Mass) + I(log10(Mass) ^ 2) + 
          log10(NPP) + I(log10(NPP) ^ 2) + 
          PCV + I(PCV ^ 2) + 
          Temp + Order,
        data = sub_d)
dredge(m)[1:10, ]
car::Anova(m, type = 3)

# model selection birds --------------
sub_d <- d %>%
  filter(Class == "Aves") %>% 
  as.data.frame()
m <- lm(log10(Density) ~ log10(Mass) + I(log10(Mass) ^ 2) + 
          log10(NPP) + I(log10(NPP) ^ 2) + 
          PCV + I(PCV ^ 2) + 
          Temp + Order,
        data = sub_d)
dredge(m)[1:10, ]
car::Anova(m, type = 3)

# model selection mammals --------------
sub_d <- d %>%
  filter(Class == "Mammalia") %>% 
  as.data.frame()
m <- lm(log10(Density) ~ log10(Mass) + I(log10(Mass) ^ 2) +
          log10(NPP) + I(log10(NPP) ^ 2) +
          PCV + I(PCV ^ 2) +
          Order + Temp,
        data = sub_d)
dredge(m)[1:10, ]
car::Anova(m, type = 3)

# model selection reptiles --------------
sub_d <- d %>%
  filter(Class == "Reptilia") %>% 
  as.data.frame()
m <- lm(log10(Density) ~ log10(Mass) + I(log10(Mass) ^ 2) +
          log10(NPP) + I(log10(NPP) ^ 2) +
          PCV + I(PCV ^ 2) +
          Order + Temp,
        data = sub_d)
dredge(m)[1:10, ]
car::Anova(m, type = 3)
