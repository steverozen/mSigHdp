
test_that("RunHdpxParallel-slow-multiT", {
  if (Sys.getenv("MSIGHDP_LONG") == "") {
    skip("Sys.setenv(MSIGHDP_LONG=\"Y\") to run test-RunHdpxParallel-multiT-slow")
    # Sys.setenv(MSIGHDP_LONG = "Y")
  }
  require(ICAMS)
  input.catalog <- PCAWG7::spectra$PCAWG$ID
  input.catalog <-
    input.catalog[ , colSums(input.catalog) > 50, drop = FALSE]
  set.seed(1)
  input.catalog <-
    input.catalog[ , sort(sample(ncol(input.catalog), 500, replace = FALSE))]

  reg <- new.env()
  load("tdata/RunHdpParallel-slow-multiT.Rdata", envir = reg)

  retvalx <- RunHdpxParallel(
    input.catalog     = input.catalog,
    CPU.cores         = 20,
    seedNumber        = 44,
    K.guess           = 5,
    multi.types       = TRUE,
    verbose           = FALSE,
    num.child.process = 20,
    burnin            = 500,
    post.space        = 50,
    post.n            = 25,
    overwrite         = TRUE,
    checkpoint        = FALSE,
    out.dir           = tempfile()

  )

  # To save a new version of the output:
  # save(retvalx, file = "tdata/RunHdpParallel-slow-multiT.Rdata")
  expect_equal(retvalx, reg$retvalx)

})
