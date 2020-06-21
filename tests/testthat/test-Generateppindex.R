test_that("Generateppindex", {

  input.catalog <-
    ICAMS::ReadCatalog("SBS96.ground.truth/1.type.ground.truth.syn.catalog.csv")

  input.catalog2 <-
    ICAMS::ReadCatalog("SBS96.ground.truth/2.type.ground.truth.syn.catalog.csv")

  reg <- new.env()
  load("RunhdpInternal.testdata/test.Generateppindex.Rdata",
       envir = reg)
  load("RunhdpInternal.testdata/test.Generateppindex2.Rdata",
       envir = reg)
  load("RunhdpInternal.testdata/test.Generateppindex3.Rdata",
       envir = reg)

  retvalx <- Generateppindex(multi.types   = TRUE,
                             input.catalog = input.catalog)

  # save(retvalx, file = "RunhdpInternal.testdata/test.Generateppindex.Rdata")

  expect_equal(retvalx, reg$retvalx)

  retvalx2 <- Generateppindex(multi.types  = FALSE,
                             input.catalog = input.catalog)

  # save(retvalx2, file = "RunhdpInternal.testdata/test.Generateppindex2.Rdata")

  expect_equal(retvalx2, reg$retvalx2)

  retvalx3 <- Generateppindex(multi.types   = TRUE,
                              input.catalog = input.catalog2)

  # save(retvalx3, file = "RunhdpInternal.testdata/test.Generateppindex3.Rdata")

  expect_equal(retvalx3, reg$retvalx3)

})
