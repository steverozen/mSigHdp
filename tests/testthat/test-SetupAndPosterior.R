
test_that("SetupAndPosterior-fast", {

  input.catalog <-
    ICAMS::ReadCatalog("SBS96.ground.truth/ground.truth.syn.catalog.csv")

  reg <- new.env()
  load("RunhdpInternal.testdata/test.SetupAndPosterior.Rdata",
       envir = reg)

  retvalx <- SetupAndPosterior(
    input.catalog = input.catalog[1:10 , 1:15],
    seedNumber    = 44,
    K.guess       = 5,
    multi.types   = FALSE,
    verbose       = TRUE,
    burnin   = 50
  )

  #save(retvalx, file = "RunhdpInternal.testdata/test.SetupAndPosterior.Rdata")

  expect_equal(retvalx, reg$retvalx)
})
