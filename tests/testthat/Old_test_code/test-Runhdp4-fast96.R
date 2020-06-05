
test_that("Runhdp4-fast96", {

  input.catalog.file <- "SBS96.ground.truth/ground.truth.syn.catalog.csv"

  regression <- new.env()
  load("RunhdpInternal.testdata/RunhdpInternal-fast96-2-cores.Rdata",
    envir = regression)

  out.dir.root <- "Runhdp3-fast96-2-cores"

  retvalx <- Runhdp4(
    input.catalog = input.catalog.file,
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
    out.dir       = "test_Runhdp4-fast96_out_dir",
    overwrite     = TRUE
  )

  #   save(retvalx, file = "tests/RunhdpInternal.testdata/RunhdpInternal-fast96-2-cores.Rdata")

  expect_equal(retvalx$signature, regression$retvalx$signature)
  expect_equal(retvalx$exposure,  regression$retvalx$exposure)

  in.catalog <- ICAMS::ReadCatalog(input.catalog.file)

  retvalx <- Runhdp4(
    input.catalog = in.catalog, # Test passing in the catalog itself
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
    out.dir       = "test_Runhdp4-fast96_out_dir",
    overwrite     = TRUE
  )

  #   save(retvalx, file = "tests/RunhdpInternal.testdata/RunhdpInternal-fast96-2-cores.Rdata")

  expect_equal(retvalx$signature, regression$retvalx$signature)
  expect_equal(retvalx$exposure,  regression$retvalx$exposure)


})
