
test_that("CombinePosteriorChainsAndExtractSigs", {

  input <- new.env()
  load("RunhdpInternal.testdata/test.CleanChlist.Rdata",
       envir = input)

  reg <- new.env()
  load("RunhdpInternal.testdata/test.CombineChainsAndExtractSigs.Rdata",
       envir = reg)

  input.catalog <-
    ICAMS::ReadCatalog("SBS96.ground.truth/ground.truth.syn.catalog.csv")


  retvalx <-
    CombineChainsAndExtractSigs(clean.chlist = input$retvalx,
                                input.catalog = input.catalog[1:10 , 1:15],
                                verbose = TRUE)

  ##Test to take output from test-GibbsSamplingAfterBurnin.R
  input2 <- new.env()
  load("RunhdpInternal.testdata/test.GibbsSamplingAfterBurnin.Rdata",
       envir = input2)

  retvalx2 <-
    CombineChainsAndExtractSigs(clean.chlist = c(list(input2$sample.chain.1),
                                                 list(input2$sample.chain.2)),
                                input.catalog = input.catalog[1:10 , 1:15],
                                verbose = TRUE)
  #save(retvalx, retvalx2,file = "RunhdpInternal.testdata/test.CombineChainsAndExtractSigs.Rdata")

  expect_equal(retvalx, reg$retvalx)
  expect_equal(retvalx2, reg$retvalx2)

})
