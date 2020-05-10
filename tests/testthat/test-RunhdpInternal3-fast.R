
test_that("RunhdpInternal3-fast", {

  input.catalog <-
    ICAMS::ReadCatalog("SBS96.ground.truth/ground.truth.syn.catalog.csv")

  load("RunhdpInternal.testdata/test.RunhdpInternal.2.Rdata")

  retvalx <- RunhdpInternal3(
    input.catalog = input.catalog[1:10 , 1:15],
    CPU.cores     = 1,
    seedNumber    = 44,
    K.guess       = 5,
    multi.types   = FALSE,
    verbose       = TRUE,
    post.burnin   = 50,
    num.posterior = 1
  )

  testthat::expect_equal(retvalx$signature, test.RunhdpInternal.2$signature)
  testthat::expect_equal(retvalx$exposure, test.RunhdpInternal.2$exposure)
})
