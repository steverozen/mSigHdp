
test_that("RunAndEvalHdp4-fast96", {
  # skip_on_cran() # Uses mulitple cores
  # skip_on_travis() # Uses multiple cores

  input.catalog.file <- "SBS96.ground.truth/ground.truth.syn.catalog.csv"

  input.exposure.file <- "SBS96.ground.truth/ground.truth.syn.exposures.csv"

  input.signature.file <- "SBS96.ground.truth/ground.truth.syn.sigs.csv"

  regression <- new.env()
  load("RunhdpInternal.testdata/RunhdpInternal-fast96-2-cores.Rdata",
    envir = regression)

  out.dir.root <- "Runhdp3-fast96-2-cores"

  retvalx <- RunAndEvalHdp4(
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
    out.dir       = file.path(out.dir.root, "test_Runhdp2-fast96_out_dir"),
    overwrite     = TRUE
  )

  # save(retvalx, file = "RunhdpInternal.testdata/RunhdpInternal-fast96-2-cores.Rdata")

  testthat::expect_equal(retvalx, regression$retvalx)
})
