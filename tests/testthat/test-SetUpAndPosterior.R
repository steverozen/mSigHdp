
test_that("RunhdpInternal4-fast", {

  input.catalog <-
    ICAMS::ReadCatalog("SBS96.ground.truth/ground.truth.syn.catalog.csv")

  reg <- new.env()
  load("RunhdpInternal.testdata/test.RunhdpInternal.2.Rdata",
       envir = reg)

  retvalx <- SetUpAndPosterior(
    input.catalog = input.catalog[1:10 , 1:15],
    seedNumber    = 44,
    K.guess       = 5,
    multi.types   = FALSE,
    verbose       = TRUE,
    post.burnin   = 50
  )

  # save(retvalx, file = "RunhdpInternal.testdata/test.RunhdpInternal.2.Rdata")

  expect_equal(retvalx$exposure, reg$retvalx$exposure)
})
