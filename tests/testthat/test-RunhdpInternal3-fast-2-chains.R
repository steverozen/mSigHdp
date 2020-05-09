
test_that("RunhdpInternal3-fast", {
  skip_on_cran()   # Because uses multiple cores
  skip_on_travis() # Because uses multiple cores

  input.catalog <-
    ICAMS::ReadCatalog(
      system.file(
        "tests/SBS96.ground.truth/ground.truth.syn.catalog.csv",
        package = "SynSigRun",
        mustWork = TRUE))

    regress <- new.env()
    load(
      system.file("tests/RunhdpInternal.testdata/test.RunhdpInternal3-fast-2-chains.Rdata",
                  package = "SynSigRun",
                  mustWork = TRUE),
      envir = regress)

  retvalx <- RunhdpInternal3(
    input.catalog = input.catalog[1:10 , 1:15],
    CPU.cores     = 2,
    seedNumber    = 44,
    K.guess       = 5,
    multi.types   = FALSE,
    verbose       = TRUE,
    post.burnin   = 50,
    num.posterior = 2
  )

  # save(retvalx, file = "tests/RunhdpInternal.testdata/test.RunhdpInternal3-fast-2-chains.Rdata")

  testthat::expect_equal(retvalx$signature, regress$retvalx$signature)
  testthat::expect_equal(retvalx$exposure,  regress$retvalx$exposure)
  testthat::expect_equal(retvalx, regress$retvalx)
})

