test_that("SetupAndPosterior-fast", {

  reg <- new.env()
  load("RunhdpInternal.testdata/test.CleanChlist.Rdata",
       envir = reg)

  if (FALSE) { # For generating the test data
    input.catalog <-
      ICAMS::ReadCatalog("SBS96.ground.truth/ground.truth.syn.catalog.csv")
    input.catalog = input.catalog[1:10,1:15]

    CPU.cores           = 2
    num.child.process   = 2

    run.setup.and.posterior <- function(seedNumber) {

      if (verbose) message("Runing run.setup.and.posterior on ", seedNumber)
      sample.chain <-SetupAndPosterior(
        input.catalog,
        seedNumber     = 44,
        K.guess        = 5,
        multi.types    = FALSE,
        verbose        = TRUE,
        post.burnin    = 50,
        post.n         = 50,
        post.space     = 5,
        post.cpiter    = 1,
        post.verbosity = 0)
      return(sample.chain)
    }

    chlist <- parallel::mclapply(
      # Must choose a different seed for each of the chains
      X = (seedNumber + 1:num.child.process * 10^6) ,
      FUN = run.setup.and.posterior,
      mc.cores = CPU.cores)

    save(chlist, file = "misc.test.data/chlist.Rdata")
  }

  c.env <- new.env()
  load("misc.test.data/chlist.Rdata", envir = c.env)

  retvalx <- CleanChlist(c.env$chlist)

  # save(retvalx, file = "RunhdpInternal.testdata/test.CleanChlist.Rdata")

  expect_equal(retvalx,reg$retvalx)

  v2 <- CleanChlist(chlist[[1]])
  expect_equal(v2, chlist[1])
  expect_error(expect_warning(CleanChlist(NULL)))

})
