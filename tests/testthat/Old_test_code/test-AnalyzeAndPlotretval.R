context("Diagnostic Plotting functions")

test_that("AnalyzeAndPlotretval", {

  reg <- new.env()
  load("RunhdpInternal.testdata/test.AnalyzeAndPlotretval.Rdata",
       envir = reg)

  input.catalog <-
    ICAMS::ReadCatalog("SBS96.ground.truth/2.type.ground.truth.syn.catalog.csv")


  AnalyzeAndPlotretval(retval = reg$retvalx,ground.truth.catalog = input.catalog[,1:15],
                       out.dir = "./test_AnalyzeAndPlotretval_out_dir")
  unlink("Rplots.pdf")
})

