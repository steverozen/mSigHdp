context("Diagnostic Plotting functions")

test_that("AnalyzeAndPlotretval", {

  reg <- new.env()
  load("RunhdpInternal.testdata/test.AnalyzeAndPlotretval.Rdata",
       envir = reg)


  AnalyzeAndPlotretval(retval = reg$retvalx,
                       out.dir = "./test_AnalyzeAndPlotretval_out_dir")
  unlink("Rplots.pdf")
})

