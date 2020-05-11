
test_that("RunhdpInternal4-slow", {
  skip_if_not(Sys.getenv("MSIGHDP_LONG_TEST") == "y")

  input.catalog <-
    ICAMS::ReadCatalog(
      "SBS96.ground.truth/ground.truth.syn.catalog.csv")

  reg <- new.env()
  load("RunhdpInternal.testdata/t2.out.Rdata", envir = reg)

  retval <- RunhdpInternal4(
    input.catalog = input.catalog[ , 1:15],
    CPU.cores     = 1,
    seedNumber    = 44,
    K.guess       = 5,
    post.burnin   = 20,
    multi.types   = FALSE,
    verbose       = TRUE,
    num.posterior = 1
  )

  save(retval, file = "RunhdpInternal.testdata/t2.out.Rdata")

  testthat::expect_equal(retval, reg$retval)
})
