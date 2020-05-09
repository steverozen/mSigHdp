
test_that("RunhdpInternal-slow", {
input.catalog <-
  ICAMS::ReadCatalog(
    system.file(
      "tests/SBS96.ground.truth/ground.truth.syn.catalog.csv",
      package = "SynSigRun",
      mustWork = TRUE))

load(
  system.file("tests/RunhdpInternal.testdata/t2.out.Rdata",
              package = "SynSigRun",
              mustWork = TRUE))

retvalx <- RunhdpInternal(
  input.catalog = input.catalog[ , 1:15],
  CPU.cores     = 1,
  seedNumber    = 44,
  K.guess       = 5,
  multi.types   = FALSE,
  verbose       = TRUE,
  num.posterior = 1
)

foo <- t2.out$signature
class(foo) <- "matrix"
testthat::expect_equal(retvalx$signature, foo, check.attributes = FALSE)
testthat::expect_equal(retvalx$exposure, t2.out$exposure)
})
