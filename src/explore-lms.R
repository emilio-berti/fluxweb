m <- lm_mamm

plot(m)

tibble(obs = m$model$`log10(Density)`,
       fitted = m$fitted.values) %>% 
  ggplot() +
  geom_point(aes(obs, fitted)) +
  geom_abline(intercept = 0, slope = 1) +
  xlim(-2, 5) +
  ylim(-2, 5) +
  theme_bw()

m <- lm_density$Mammalia

d %>%
  mutate(predMass = log10(Mass) * m["log10(Mass)"] + log10(Mass) ^ 2 * m["I(log10(Mass)^2)"]) %>%
  mutate(predPCV = PCV * m["PCV"] + PCV ^ 2 * m["I(PCV^2)"]) %>%
  mutate(predNPP = log10(NPP) * m["log10(NPP)"] + log10(NPP) ^ 2 * m["I(log10(NPP)^2)"]) %>%
  pivot_longer(cols = predMass : predNPP, names_to = "fit") %>%
  pivot_longer(cols = c(PCV, NPP, Mass), names_to = "covariate", values_to = "x") %>%
  mutate(fit = gsub("pred", "", fit),
         covariate = gsub("[log10(|)]", "", covariate)) %>%
  filter(covariate == fit) %>%
  ggplot() +
  geom_line(aes(x, value)) +
  geom_rug(aes(x)) +
  facet_wrap(~covariate, scales = "free") +
  scale_x_log10() +
  theme_bw()
