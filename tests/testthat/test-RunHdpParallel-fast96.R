
test_that("RunHdpParallel-fast96", {

  input.catalog <-
    ICAMS::ReadCatalog("SBS96.ground.truth/ground.truth.syn.catalog.csv")

  regression <- new.env()
  load("RunhdpInternal.testdata/RunHdpParallel-fast96-2-cores.Rdata",
       envir = regression)

  retvalx <- RunHdpParallel(
    input.catalog = input.catalog[1:10,1:15],
    CPU.cores     = 2,
    seedNumber    = 44,
    K.guess       = 5,
    multi.types   = FALSE,
    verbose       = TRUE,
    num.child.process =  2,
    post.burnin   = 50, # Super low for fast testing
    post.space    = 5,  # Low for fast testing
    post.cpiter   = 1,  # Low for fast testing
    overwrite     = TRUE
  )

  # save(retvalx, file = "RunhdpInternal.testdata/RunHdpParallel-fast96-2-cores.Rdata")
  expect_equal(retvalx$signature, regression$retvalx$signature)
  expect_equal(retvalx$exposure, regression$retvalx$exposure)
  expect_equal(retvalx$multi.chains, regression$retvalx$multi.chains)
})
