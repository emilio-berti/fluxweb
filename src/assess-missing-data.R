library(tidyverse)
library(sf)

birds <- read_delim("data/raw/lislevand.csv", delim = "\t", show_col_types = FALSE) %>% 
  select(F_mass, F_wing) %>% 
  mutate(across(where(is.numeric), function(.x) {ifelse(.x <= 0, NA, log10(.x))}))

# p <- st_read("data/interim/tetradensity_masses.shp")
# 
# d %>% 
#   mutate(Mass = ifelse(!is.na(Mass) & Mass > 0, "Valid data", "Missing data")) %>%
#   group_by(Class, Mass) %>%
#   tally() %>%
#   ungroup() %>% 
#   ggplot() +
#   aes(x = Class, y = n, fill = Mass) +
#   geom_col(position = "stack", width = .5) +
#   scale_y_log10() +
#   scale_fill_manual(values = c("blue4", "red4")) +
#   xlab("") +
#   ylab("Number of data points") +
#   theme_classic()
# 
# d %>% 
#   mutate(Mass = ifelse(!is.na(Mass) & Mass > 0, "Valid data", "Missing data")) %>%
#   group_by(Class, Order, Mass) %>%
#   tally() %>%
#   ungroup() %>% 
#   ggplot() +
#   aes(x = Order, y = n, fill = Mass) +
#   geom_col(position = "stack", width = .5) +
#   scale_y_log10() +
#   scale_fill_manual(values = c("blue4", "red4")) +
#   xlab("") +
#   ylab("Number of data points") +
#   theme_classic() +
#   facet_wrap(~Class, scales = "free", drop = TRUE) +
#   theme(axis.text.x = element_text(angle = 90, hjust= 1, vjust = .5),
#         strip.background = element_blank())
# 
# d %>% 
#   filter(Class == "Aves") %>% 
#   mutate(Mass = ifelse(!is.na(Mass) & Mass > 0, "Valid data", "Missing data")) %>%
#   group_by(Order, Family, Mass) %>%
#   tally() %>%
#   ungroup() %>% 
#   ggplot() +
#   aes(x = Family, y = n, fill = Mass) +
#   geom_col(position = "stack", width = .5) +
#   scale_y_log10() +
#   scale_fill_manual(values = c("blue4", "red4")) +
#   xlab("") +
#   ylab("Number of data points") +
#   theme_classic() +
#   facet_wrap(~Order, scales = "free", drop = TRUE) +
#   theme(axis.text.x = element_text(angle = 90, hjust= 1, vjust = .5),
#         strip.background = element_blank())

# Missing data affects particularly birds and reptiles, with few orders
# specifically affected. This is a 'missing not at random' scenario and
# imputation should be performed.

# multi-imputation using `mice`.
Y <- birds %>% 
  distinct_all() %>% 
  as.data.frame()

md.pairs(birds)

mi <- mice(data = Y, m = 20, maxit = 10, method = c("norm", "norm")) #impute only mass

# stripplot(mi, pch = 20, cex = 2) #red points are imputed
# xyplot(mi, F_mass ~ F_wing | .imp, pch = 20, cex = 1)

plot(mi) # for convergence

densityplot(mi, scales = list(x = list(relation = "free")), thicker = 5)

mass.na <- is.na(Y$F_mass)
fit.mass <- with(mi, lm(F_mass ~ F_wing))
ps <- rep(rowMeans(sapply(fit.mass$analyses, fitted.values)), mi$m + 1)
xyplot(mi, F_mass ~ F_wing | .imp, pch = c(1, 20), cex = c(.5, .5))

mass <- complete(mi, "long", TRUE)$F_mass
fit <- lm(mass[!is.na(mass)] ~ ps[!is.na(mass)])
densityplot(~residuals(fit), group = is.na(mass), plot.points = FALSE, #show overlap largely
            ref = TRUE, scales = list(y = list(draw = FALSE)),
            xlab = "Residuals mass ~ wing",
            lwd = 2)

fit <- with(mi, lm(F_mass ~ F_wing))
pooled <- pool(fit$analyses)

### very good imputation
gamma <- pooled$pooled$fmi
message("Influence of missing data on estimate uncertainty (proportion): ", paste(round(gamma, 3), collapse = "; "))
epsilon <- ( 1 + gamma / mi$m ) ^ - 1
message("Relative efficiency: ", paste(round(epsilon, 3), collapse = "; ")) #aim for over 99%
### 


analysis <- with (mi, lm (Mass ~ Order)) # the analysis step with a GLM (see Chapters 6 and 12)
pooling <- pool (analysis) # the pooling step

ps <- rep(rowMeans(sapply(analysis$analyses, fitted.values)), 6)
xyplot(mi, Mass ~ ps | .imp, pch = 20, cex = 3)

mi$imp$Mass



stripplot(pooling, pch = 20, cex = 3) #red points are imputed

summary(pooling)
