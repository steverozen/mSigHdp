# example of using PlotExposureByRange

pdf(path, width = 8.2677, height = 11.6929, # for A4
    onefile = TRUE, useDingbats = FALSE)
PlotExposureByRange(
  exposures       = exposures,
  num.per.line    = 50,
  plot.proportion = FALSE,
  main = "Just an example")
dev.off()
