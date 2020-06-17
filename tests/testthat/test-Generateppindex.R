test_that("ActivateAndBurnin-fast", {

  input.catalog <-
    ICAMS::ReadCatalog("SBS96.ground.truth/ground.truth.syn.catalog.csv")

  reg <- new.env()
  load("RunhdpInternal.testdata/test.Generateppindex.Rdata",
       envir = reg)

  retvalx <- Generateppindex(multi.types     = F,
                             one.parent.hack = T,
                             input.catalog   = input.catalog
  )

  #save(retvalx, file = "RunhdpInternal.testdata/test.Generateppindex.Rdata")

  expect_equal(retvalx$ppindex, reg$retvalx$ppindex)
  expect_equal(retvalx$num.tumor.types, reg$retvalx$num.tumor.types)

})
