
test_that("CombinePosteriorChains", {

  input.catalog <-
    ICAMS::ReadCatalog("SBS96.ground.truth/ground.truth.syn.catalog.csv")

  reg <- new.env()
  load("RunhdpInternal.testdata/test.CombinePosteriorChains.Rdata",
       envir = reg)

  clean.chlist <- MultipleSetupAndPosterior( input.catalog = input.catalog[1:10 , 1:15],
                                        seedNumber          = 44,
                                        K.guess             = 5,
                                        multi.types         = FALSE,
                                        verbose             = TRUE,
                                        post.burnin         = 50,
                                        post.n              = 5,
                                        post.space          = 50,
                                        post.cpiter         = 3,
                                        post.verbosity      = 0,
                                        CPU.cores           = 2,
                                        num.child.process   = 2)


  retvalx <- CombinePosteriorChains(clean.chlist = clean.chlist,
                                    input.catalog = input.catalog[1:10 , 1:15],
                                    multi.types = FALSE,
                                    verbose = TRUE)



  #save(retvalx, file = "RunhdpInternal.testdata/test.CombinePosteriorChains.Rdata")

  expect_equal(retvalx$signature, reg$retvalx$signature)
  expect_equal(retvalx$exposure, reg$retvalx$exposure)
  expect_equal(retvalx$multi.chains, reg$retvalx$multi.chains)
})
