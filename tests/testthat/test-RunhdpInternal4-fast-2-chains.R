
test_that("RunhdpInternal4-fast-2-chains", {
  # skip_on_cran()   # Because uses multiple cores
  # skip_on_travis() # Because uses multiple cores

  input.catalog <-
    ICAMS::ReadCatalog("SBS96.ground.truth/ground.truth.syn.catalog.csv")

  reg <- new.env()
  load("RunhdpInternal.testdata/test.RunhdpInternal3-fast-2-chains.Rdata",
       envir = reg)

  retvalx <- RunhdpInternal4(
    input.catalog = input.catalog[1:10 , 1:15],
    CPU.cores     = 2,
    seedNumber    = 44,
    K.guess       = 5,
    multi.types   = FALSE,
    verbose       = TRUE,
    post.burnin   = 50,
    num.posterior = 2
  )

  # save(retvalx, file = "RunhdpInternal.testdata/test.RunhdpInternal3-fast-2-chains.Rdata")

  testthat::expect_equal(retvalx, reg$retvalx)
})

