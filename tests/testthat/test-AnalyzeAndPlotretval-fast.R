
test_that("AnalyzeAndPlotretval", {

  input <- new.env()
  load("RunhdpInternal.testdata/test.CombineChainsAndExtractSigs.Rdata",
       envir=input)

  reg <- new.env()
  load("RunhdpInternal.testdata/test.AnalyzeAndPlotretval.Rdata",
       envir = reg)

  input.catalog <-
    ICAMS::ReadCatalog("SBS96.ground.truth/ground.truth.syn.catalog.csv")


  retvalx <-
    AnalyzeAndPlotretval(retval = input$retvalx,
                         out.dir = "./tmp",
                         diagnostic.plot = T,
                         overwrite = T,
                         input.catalog = input.catalog[1:10 , 1:15]
                               )
  unlink("./tmp",recursive = T) #remove the out.dir to avoid warning

  #save(retvalx, file = "RunhdpInternal.testdata/test.AnalyzeAndPlotretval.Rdata")

  #expect_equal(retvalx, reg$retvalx)

})
