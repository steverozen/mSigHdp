
test_that("SaveAnalysis", {

  input <- new.env()
  load("RunhdpInternal.testdata/test.CombineChainsAndExtractSigs.Rdata",
       envir=input)

  input.catalog <-
    ICAMS::ReadCatalog("SBS96.ground.truth/ground.truth.syn.catalog.csv")

  retvalx <-
    SaveAnalysis(retval = input$retvalx,
                 out.dir = "./tmp",
                 diagnostic.plot = T,
                 overwrite = T,
                 input.catalog = input.catalog[1:10 , 1:15]
    )
  unlink("./tmp",recursive = T) #remove the out.dir to avoid warning

  expect_null(retvalx)

})
