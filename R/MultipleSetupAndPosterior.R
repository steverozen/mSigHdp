#' Activate hierarchical Dirichlet processes and run posterior sampling in parallel.
#'
#' @param CPU.cores Number of CPUs to use; there is no
#'    point in making this larger than \code{num.child.process}.
#'
#' @param num.child.process Number of posterior sampling chains; can set to
#'   1 for testing.
#'
#' @inheritParams SetupAndPosterior
#'
#' @return Invisibly,
#'    the clean \code{chlist} (output of the hdp_posterior calls).
#'
#' @export


MultipleSetupAndPosterior <- function(input.catalog,
                                      seedNumber          = 1,
                                      K.guess,
                                      multi.types         = FALSE,
                                      verbose             = TRUE,
                                      post.burnin         = 4000,
                                      post.n              = 50,
                                      post.space          = 50,
                                      post.cpiter         = 3,
                                      post.verbosity      = 0,
                                      CPU.cores           = 1,
                                      num.child.process   = 4){

  run.setup.and.posterior <- function(seedNumber) {
<<<<<<< HEAD
    if (verbose) message("Runing run.setup.and.posterior on ", my.seed)
=======
    if (verbose) message("Runing run.setup.and.posterior on ", seedNumber)
>>>>>>> bdece76336b40ca7c7f180f6249597b98bfa79b6
    sample.chain <-SetupAndPosterior(
      input.catalog,
      seedNumber     = seedNumber,
      K.guess        = K.guess,
      multi.types    = multi.types,
      verbose        = verbose,
      post.burnin    = post.burnin,
      post.n         = post.n,
      post.space     = post.space,
      post.cpiter    = post.cpiter,
      post.verbosity = post.verbosity)
    return(sample.chain)
  }

  chlist <- parallel::mclapply(
    # Must choose a different seed for each of the chains
    X = (seedNumber + 1:num.child.process * 10^6) ,
    FUN = run.setup.and.posterior,
    mc.cores = CPU.cores)

  clean.chlist <- CleanChlist(chlist)

  return(invisible(clean.chlist))
}
