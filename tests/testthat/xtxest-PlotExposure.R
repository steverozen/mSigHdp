

ex <- SynSigGen::ReadExposure("tests/testthat/SBS96.ground.truth/ground.truth.syn.exposures.csv")

par(mfrow = c(2, 1))
PlotExposure(SortExp(ex[, 1:30]))
PlotExposure(SortExp(ex[ , 1:30]), plot.proportion = TRUE)

par(mfcol = c(2, 3))
PlotExposureByRange(exp = SortExp(ex[, 1:43 ]),
                    num.per.line = 20, main = "xxxxx")

PlotExposureByRange(exp = SortExp(ex[, 3:6 ]),
                    num.per.line = 20)

par(mfcol = c(2, 3))
PlotExposureByRange(exp = SortExp(ex[, 1:43 ]),
                    num.per.line = 20,
                    plot.proportion = TRUE)

PlotExposureByRange(exp = SortExp(ex[, 3:6 ]),
                    num.per.line = 20,
                    plot.proportion = TRUE, col = c("red", "white"))

