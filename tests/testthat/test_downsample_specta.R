test_that("downsample_spectra", {
  reg <- new.env()
  load("tdata/downsample_spectra_and_result.Rdata", envir = reg)
  res <- downsample_spectra(reg$test.ex, downsample_threshold = 3000)
  expect_equal(res, reg$res)
})


