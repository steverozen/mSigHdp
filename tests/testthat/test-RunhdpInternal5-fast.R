
test_that("RunhdpInternal5-fast", {

  input.catalog <-
    ICAMS::ReadCatalog("SBS96.ground.truth/ground.truth.syn.catalog.csv")

  reg <- new.env()
  load("RunhdpInternal.testdata/test.RunhdpInternal.2.Rdata",
       envir = reg)

  retvalx <- RunhdpInternal5(
    input.catalog = input.catalog[1:10 , 1:15],
    CPU.cores     = 1,
    seedNumber    = 44,
    K.guess       = 5,
    multi.types   = FALSE,
    verbose       = TRUE,
    post.burnin   = 50,
    num.posterior = 1
  )

  expect_equal(retvalx$signature, reg$retvalx$signature)

  expect_equal(retvalx$exposure, reg$retvalx$exposure)


})
