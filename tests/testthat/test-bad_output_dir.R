
test_that("Bad output dir", {

  input.catalog <-
    ICAMS::ReadCatalog("SBS96.ground.truth/ground.truth.syn.catalog.csv")

  expect_error(
    retvalx <- RunHdpxParallel(input.catalog = input.catalog[1:10 , 1:15],
                                out.dir      = foobar)
  )

  expect_error(
    retvalx <- RunHdpxParallel(input.catalog = input.catalog[1:10 , 1:15],
                               out.dir       = "/yyy/")
  )

  retvalx <- RunHdpxParallel(
    input.catalog     = input.catalog[1:10,1:15],
    CPU.cores         = 2,
    seedNumber        = 44,
    K.guess           = 5,
    multi.types       = FALSE,
    verbose           = FALSE,
    num.child.process =  2,
    burnin            = 50, # Super low for fast testing
    post.space        = 5,  # Low for fast testing
    post.cpiter       = 1,  # Low for fast testing
    overwrite         = TRUE,
    checkpoint        = FALSE)

  new.dirs <- dir(".", pattern = "^RunHdpxParallel_out_\\d", full.names = TRUE)
  expect_equal(length(new.dirs), 1)
  unlink(new.dirs, recursive = TRUE)


})
