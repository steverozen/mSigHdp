
test_that("ActivateAndBurnin-one-multitype", {

  input.catalog <-
    ICAMS::ReadCatalog(
      "SBS96.ground.truth/1.type.ground.truth.syn.catalog.csv")

  reg <- new.env()
  load("RunhdpInternal.testdata/test.ActivateAndBurnin.Rdata",
       envir = reg)

  retvalx <- ActivateAndBurnin(input.catalog = input.catalog[1:10 , 1:15],
                             seedNumber          = (44 + 3e6),
                             K.guess             = 5,
                             multi.types         = TRUE,
                             verbose             = TRUE,
                             burnin              = 100,
                             cpiter              = 3,
                             burnin.verbosity    = 0,
                             gamma.alpha         = 1,
                             gamma.beta          = 1,
                             one.parent.hack     = TRUE
  )

  #save(retvalx, file = "RunhdpInternal.testdata/test.ActivateAndBurnin.Rdata")

  expect_equal(retvalx$hdplist, reg$retvalx$hdplist)
  expect_equal(retvalx$likelihood, reg$retvalx$likelihood)

})
