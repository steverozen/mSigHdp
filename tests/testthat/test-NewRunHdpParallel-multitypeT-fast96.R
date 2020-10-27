
test_that("RunHdpParallel-multi-type-T-fast96", {

  input.catalog <-
    ICAMS::ReadCatalog("SBS96.ground.truth/2.type.ground.truth.syn.catalog.csv")

  reg <- new.env()
  load("RunhdpInternal.testdata/NewRunHdpParallel-fast96-2-cores-multi-type-T.Rdata",
       envir = reg)

  retvalx1 <- RunHdpxParallel(
    input.catalog = input.catalog[1:10,1:15],
    CPU.cores         = 2,
    seedNumber        = 44,
    K.guess           = 5,
    multi.types       = TRUE,
    verbose           = TRUE,
    num.child.process = 2,
    post.burnin       = 100, # Super low for fast testing
    post.n            = 2,
    post.space        = 5,  # Low for fast testing
    post.cpiter       = 1,  # Low for fast testing
    overwrite         = TRUE,
    burnin.multiplier = 2,
    burnin.checkpoint = T
  )

  retvalx2 <- RunHdpxParallel(
    input.catalog = input.catalog[1:10,1:15],
    CPU.cores         = 2,
    seedNumber        = 44,
    K.guess           = 5,
    multi.types       = TRUE,
    verbose           = TRUE,
    num.child.process = 2,
    post.burnin       = 200, # Super low for fast testing
    post.n            = 2,
    post.space        = 5,  # Low for fast testing
    post.cpiter       = 1,  # Low for fast testing
    overwrite         = TRUE,
    burnin.multiplier = 1,
    burnin.checkpoint = T
  )

  #save(retvalx1, file = "RunhdpInternal.testdata/NewRunHdpParallel-fast96-2-cores-multi-type-T.Rdata")
  expect_equal(retvalx1, reg$retvalx1)
  expect_equal(retvalx2, reg$retvalx1)

})
