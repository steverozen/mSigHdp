
test_that("RunHdpxParallel-fast96-downsample", {

  input.catalog <-
    ICAMS::ReadCatalog("SBS96.ground.truth/ground.truth.syn.catalog.csv")

  reg <- new.env()
  load("RunhdpInternal.testdata/NewRunHdpParallel-fast96-2-cores.Rdata",
       envir = reg)

  retvalx <- RunHdpxParallel(
    input.catalog         = input.catalog[1:10,1:15],
    CPU.cores             = 2,
    seedNumber            = 44,
    K.guess               = 5,
    multi.types           = FALSE,
    verbose               = FALSE,
    num.child.process     =  2,
    burnin                = 50, # Super low for fast testing
    post.space            = 5,  # Low for fast testing
    post.cpiter           = 1,  # Low for fast testing
    overwrite             = TRUE,
    checkpoint            = FALSE,
    downsample_threshold  = 1e6, # Very high; should have no effect
    out.dir               = tempfile()
  )

  if (FALSE) { # To regenerate test data
    save(retvalx,
         file = "RunhdpInternal.testdata/NewRunHdpParallel-fast96-2-cores.Rdata")
  }
  expect_equal(retvalx, reg$retvalx)


  retvalx2 <- RunHdpxParallel(
    input.catalog         = input.catalog[1:10,1:15],
    CPU.cores             = 2,
    seedNumber            = 44,
    K.guess               = 5,
    multi.types           = FALSE,
    verbose               = FALSE,
    num.child.process     =  2,
    burnin                = 50, # Super low for fast testing
    post.space            = 5,  # Low for fast testing
    post.cpiter           = 1,  # Low for fast testing
    overwrite             = TRUE,
    checkpoint            = FALSE,
    downsample_threshold  = 1e3, # Unrealistically low because this is toy data
    out.dir               = tempfile()
  )

  if (FALSE) { # To regenerate test data
    save(retvalx2,
         file = "RunhdpInternal.testdata/NewRunHdpParallel-fast96-2-cores-downsample.Rdata")
  }

  reg2 <- new.env()
  load("RunhdpInternal.testdata/NewRunHdpParallel-fast96-2-cores-downsample.Rdata",
       envir = reg2)
  expect_equal(retvalx2, reg2$retvalx2)


})
