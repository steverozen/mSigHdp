
test_that("RunhdpInternal-fast96", {

  input.catalog <-
    ICAMS::ReadCatalog(
      system.file(
        "tests/SBS96.ground.truth/ground.truth.syn.catalog.csv",
        package = "SynSigRun",
        mustWork = TRUE))

  regression <- new.env()
  load(
    system.file("tests/RunhdpInternal.testdata/RunhdpInternal-fast96.Rdata",
                package = "SynSigRun",
                mustWork = TRUE),
    envir = regression)

  retvalx <- RunhdpInternal(
    input.catalog = input.catalog[, 1:10],
    CPU.cores     = 1,
    seedNumber    = 44,
    K.guess       = 5,
    multi.types   = FALSE,
    verbose       = TRUE,
    num.posterior = 1,
    post.burnin   = 50, # Super low for fast testing
    post.space    = 5,  # Low for fast testing
    post.cpiter   = 1   # Low for fast testing
  )

  testthat::expect_equal(retvalx, regression$retvalx)
})
