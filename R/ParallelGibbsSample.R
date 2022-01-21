#' Setup hierarchical Dirichlet processes and run parallel Gibbs sampling chains
#'
#' @inheritParams SetupAndPosterior
#'
#' @param CPU.cores Number of CPUs to use; this should be no
#'   more than \code{num.child.process}.
#'
#' @param num.child.process Number of posterior sampling chains; can set to
#'   1 for testing. We recommend 20 for real data analysis
#'
#' @return Invisibly,
#'    the clean \code{chlist} (output of \code{CleanChlist}).
#'    This is a list of \code{\link[hdpx]{hdpSampleChain-class}}
#'    objects (see package hdpx).
#'
#' @export


ParallelGibbsSample <- function(input.catalog,
                                seedNumber          = 1,
                                K.guess,
                                multi.types         = FALSE,
                                verbose             = TRUE,
                                burnin              = 5000,
                                burnin.multiplier   = 2,
                                post.n              = 200,
                                post.space          = 100,
                                post.cpiter         = 3,
                                post.verbosity      = 0,
                                CPU.cores           = 20,
                                num.child.process   = 20,
                                gamma.alpha         = 1,
                                gamma.beta          = 20,
                                checkpoint          = TRUE,
                                prior.sigs          = NULL,
                                prior.pseudoc       = NULL) {


  run.setup.and.posterior <- function(seedNumber) {

    if (verbose) message("Runing run.setup.and.posterior on ", seedNumber)
    sample.chain <- SetupAndPosterior(
      input.catalog,
      seedNumber          = seedNumber,
      K.guess             = K.guess,
      multi.types         = multi.types,
      verbose             = verbose,
      burnin              = burnin,
      post.n              = post.n,
      post.space          = post.space,
      post.cpiter         = post.cpiter,
      post.verbosity      = post.verbosity,
      gamma.alpha         = gamma.alpha,
      gamma.beta          = gamma.beta,
      prior.sigs          = prior.sigs,
      prior.pseudoc       = prior.pseudoc,
      burnin.multiplier   = burnin.multiplier,
      checkpoint          = checkpoint)
    return(sample.chain)
  }

  chlist <- parallel::mclapply(
    # Must choose a different seed for each of the chains
    X = (seedNumber + 1:num.child.process * 10^6) ,
    FUN = run.setup.and.posterior,
    mc.cores = CPU.cores)

  if (FALSE) { # For debugging
    save(chlist, file = paste0("initial.chlist.Rdata"))
  }

  clean.chlist <- CleanChlist(chlist, verbose)

  return(invisible(clean.chlist))
}
