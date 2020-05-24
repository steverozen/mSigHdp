context("Exposure plotting functions")

test_that("PlotExposure function", {
ex <-
  SynSigGen::ReadExposure("SBS96.ground.truth/ground.truth.syn.exposures.csv")

par(mfrow = c(2, 1), mar = c(3, 1, 2, 1), oma = c(2, 3, 0, 0))
PlotExposure(SortExp(ex[, 1:30]))
PlotExposure(SortExp(ex[ , 1:30]), plot.proportion = TRUE)
})

test_that("PlotExposureByRange function", {
  ex <-
    SynSigGen::ReadExposure("SBS96.ground.truth/ground.truth.syn.exposures.csv")

  par(mfcol = c(2, 2), mar = c(2, 3, 4, 2), oma = c(4, 1, 1, 1))
  PlotExposureByRange(exp = SortExp(ex[, 1:43 ]),
                      num.per.line = 20, main = "xxxxx")

  PlotExposureByRange(exp = SortExp(ex[, 3:6 ]),
                      num.per.line = 20)

  par(mfcol = c(2, 2), mar = c(2, 3, 4, 2), oma = c(4, 1, 1, 1))
  PlotExposureByRange(exp = SortExp(ex[, 1:43 ]),
                      num.per.line = 20,
                      plot.proportion = TRUE)

  PlotExposureByRange(exp = SortExp(ex[, 3:6 ]),
                      num.per.line = 20,
                      plot.proportion = TRUE, col = c("red", "white"))
})
