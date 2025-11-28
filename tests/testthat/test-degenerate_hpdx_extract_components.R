test_that("degenerate_hdpx_extract_components", {
  testthat::local_edition(3)

  input.catalog <-
    read.csv("SBS96.ground.truth/ground.truth.syn.catalog.csv")[1:10, c(3, 3)]

  retvalx <- suppressWarnings({
    RunHdpxParallel(
      input.catalog = input.catalog,
      CPU.cores = 2,
      seedNumber = 44,
      K.guess = 5,
      multi.types = FALSE,
      verbose = FALSE,
      num.child.process = 2,
      burnin = 50, # Super low for fast testing
      post.space = 5, # Low for fast testing
      post.cpiter = 1, # Low for fast testing
      overwrite = TRUE,
      checkpoint = FALSE,
      downsample_threshold = 1e6, # Very high; should have no effect
      out.dir = tempfile()
    )
  })

  expect_snapshot(retvalx)
})

test_that("degenerate_hdpx_extract_components_warnings", {
  testthat::local_edition(3)

  input.catalog <-
    read.csv("SBS96.ground.truth/ground.truth.syn.catalog.csv")[1:10, c(3, 3)]

  expect_snapshot_warning(
    RunHdpxParallel(
      input.catalog = input.catalog,
      CPU.cores = 2,
      seedNumber = 44,
      K.guess = 5,
      multi.types = FALSE,
      verbose = FALSE,
      num.child.process = 2,
      burnin = 50, # Super low for fast testing
      post.space = 5, # Low for fast testing
      post.cpiter = 1, # Low for fast testing
      overwrite = TRUE,
      checkpoint = FALSE,
      downsample_threshold = 1e6, # Very high; should have no effect
      out.dir = tempfile()
    )
  )
})
