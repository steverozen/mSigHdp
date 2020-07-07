
test_that("CombinePosteriorChains", {

  input <- new.env()
  load("RunhdpInternal.testdata/test.CleanChlist.Rdata",
       envir=input)

  reg <- new.env()
  load("RunhdpInternal.testdata/test.CombinePosteriorChains.Rdata",
       envir = reg)

  input.catalog <-
    ICAMS::ReadCatalog("SBS96.ground.truth/ground.truth.syn.catalog.csv")


  retvalx <-
    CombinePosteriorChains(clean.chlist = input$retvalx,
                           input.catalog = input.catalog[1:10 , 1:15],
                           cluster.method = "kmedians",
                           multi.types = FALSE,
                           verbose = TRUE)

  #save(retvalx, file = "RunhdpInternal.testdata/test.CombinePosteriorChains.Rdata")

  expect_equal(retvalx, reg$retvalx)

})
