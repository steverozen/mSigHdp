
test_that("RunAndEvalHdp5-fast96", {
  # skip_on_cran() # Uses mulitple cores
  # skip_on_travis() # Uses multiple cores

  input.catalog.file <- "SBS96.ground.truth/ground.truth.syn.catalog.csv"

  input.exposure.file <- "SBS96.ground.truth/ground.truth.syn.exposures.csv"

  input.signature.file <- "SBS96.ground.truth/ground.truth.syn.sigs.csv"

  regression <- new.env()
  load("RunhdpInternal.testdata/RunhdpInternal-fast96-2-cores.Rdata",
    envir = regression)

  retvalx <- RunAndEvalHdp5(
    input.catalog.file         = input.catalog.file, # The spectra
    ground.truth.exposure.file = input.exposure.file,
    ground.truth.sig.file      = input.signature.file,
    test.only     = 10, # Only use columns 1:10 of input.catalog
    CPU.cores     = 2,
    seedNumber    = 44,
    K.guess       = 5,
    multi.types   = FALSE,
    verbose       = TRUE,
    num.posterior = 2,
    post.burnin   = 50, # Super low for fast testing
    post.space    = 5,  # Low for fast testing
    post.cpiter   = 1,  # Low for fast testing
    out.dir       = "test-RunAnd..etc5_out_dir",
    overwrite     = TRUE
  )

  expect_equal(retvalx$signature, regression$retvalx$signature)
  expect_equal(retvalx$exposure, regression$retvalx$exposure)
})
