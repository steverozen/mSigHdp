
test_that("BurninIteration-fast", {

  input.catalog <-
    ICAMS::ReadCatalog("SBS96.ground.truth/ground.truth.syn.catalog.csv")

  reg <- new.env()
  load("RunhdpInternal.testdata/test.BurninIteration.Rdata",
       envir = reg)

  retvalx <- BurninIteration(input.catalog = input.catalog[1:10 , 1:15],
             seedNumber          = 44,
             K.guess             =5,
             multi.types         = FALSE,
             verbose             = TRUE,
             post.burnin         = 100,
             post.cpiter         = 3,
             post.verbosity      = 0

            )

  #save(retvalx, file = "RunhdpInternal.testdata/test.BurninIteration.Rdata")

  expect_equal(retvalx, reg$retvalx)
})
