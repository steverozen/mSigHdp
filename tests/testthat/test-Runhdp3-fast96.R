
test_that("Runhdp3-fast96", {
  skip_on_cran() # Uses mulitple cores
  skip_on_travis() # Uses multiple cores

  input.catalog.file <-
    system.file(
      "tests/SBS96.ground.truth/ground.truth.syn.catalog.csv",
      package = "SynSigRun",
      mustWork = TRUE)

  regression <- new.env()
  load(
    system.file("tests/RunhdpInternal.testdata/RunhdpInternal-fast96-2-cores.Rdata",
                package = "SynSigRun",
                mustWork = TRUE),
    envir = regression)

  out.dir.root <- system.file("tests/Runhdp3-fast96-2-cores",
                              package = "SynSigRun",
                              mustWork = TRUE)

  retvalx <- Runhdp3(
    input.catalog.file = input.catalog.file,
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

  #   save(retvalx, file = "tests/RunhdpInternal.testdata/RunhdpInternal-fast96-2-cores.Rdata")

  testthat::expect_equal(retvalx, regression$retvalx)
})
