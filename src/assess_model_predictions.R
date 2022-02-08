library(tidyverse)
library(sf)

source("make_density_models.R")

# check assumptions of a single model ---------
cairo_pdf("docs/plots/lm_density_residuals.pdf", width = 8, height = 6, onefile = TRUE)
par(mfrow = c(2, 2), mar = c(4, 4, 2, 2))
for (x in c("Amphibia", "Aves", "Mammalia", "Reptilia")) {
  m <- lm_density[[x]]
  plot(m, pch = 20)
}
dev.off()
  
# plot pred vs obs densities ---------
cairo_pdf("docs/plots/lm_density_predictions.pdf", width = 10, height = 6, onefile = TRUE)
par(mfrow = c(2, 2), mar = c(4, 4, 2, 2))
for (x in c("Amphibia", "Aves", "Mammalia", "Reptilia")) {
  m <- lm_density[[x]]
  scatter.smooth(m$model$`log10(Density)`,
                 m$fitted.values,
                 frame = FALSE,
                 pch = 20,
                 col = factor(m$model$Order),
                 main = x,
                 xlab = "log10(density) from TetraDENSITY",
                 ylab = "Predicted log10(density)")
  abline(0, 1, lty = 2, col = "blue", lw = 2)
}
dev.off()

x = "Mammalia"
m = lm_density[[x]]
par(mfrow = c(1, 2), mar = c(8, 4, 2, 2))
boxplot(m$model$`log10(Density)` ~ m$model$Order, main = "TetraDENSITY",
        xlab = "", ylab = "Density", las = 2)
grid(lty = 5)
boxplot(m$fitted.values ~ m$model$Order, main = "Fitted",
        xlab = "", ylab = "", las = 2)
grid(lty = 5)

