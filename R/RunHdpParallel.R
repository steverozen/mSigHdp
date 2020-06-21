#' Extract mutational signatures and optionally compare them to existing signatures and exposures.
#'
#' @inheritParams MultipleSetupAndPosterior
#'
#' @inheritParams AnalyzeAndPlotretval
#'
#' @inheritParams CombinePosteriorChains
#'
#' @return Invisibly, a list with the following elements:\describe{
#' \item{signature}{The extracted signature profiles as a matrix;
#'             rows are mutation types, columns are
#'             samples (e.g. tumors).}
#'
#' \item{exposure}{The inferred exposures as a matrix of mutation counts;
#'            rows are signatures, columns are samples (e.g. tumors).}
#'
#' \item{multi.chains}{A \code{\link[hdpx]{hdpSampleMulti-class}} object.
#'     This object has the method \code{\link[hdpx]{chains}} which returns
#'     a list of \code{\link[hdpx]{hdpSampleChain-class}} objects. Each of these
#'     sample chains objects has a method \code{\link[hdpx]{final_hdpState}}
#'     (actually the methods seems to be just \code{hdp})
#'     that returns the \code{hdpState} from which it was generated.}
#' }
#'
#'
#' @export

RunHdpParallel <- function(input.catalog,
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
                           num.child.process   = 4,
                           cos.merge           = 0.9,
                           min.sample          = 1,
                           ground.truth.sig    = NULL,
                           ground.truth.exp    = NULL,
                           overwrite           = TRUE,
                           out.dir             = NULL,
                           gamma.alpha         = 1,
                           gamma.beta          = 1){

  # Step 1: Activate hierarchical Dirichlet processes and
  # run posterior sampling in parallel;
  # chlist is a list of hdpSampleChain-class objects.

  chlist <-
    MultipleSetupAndPosterior(input.catalog,
                              seedNumber          = seedNumber,
                              K.guess             = K.guess,
                              multi.types         = multi.types,
                              verbose             = verbose,
                              post.burnin         = post.burnin,
                              post.n              = post.n,
                              post.space          = post.space,
                              post.cpiter         = post.cpiter,
                              post.verbosity      = post.verbosity,
                              CPU.cores           = CPU.cores,
                              num.child.process   = num.child.process,
                              gamma.alpha         = gamma.alpha,
                              gamma.beta          = gamma.beta)

  # Step 2: Combine the posterior chains and extract
  # signatures and exposures;
  # retval has signatures, exposures, and multi.chains, a
  # hdpSampleMulti-class object.

  retval <-
    CombinePosteriorChains(chlist,
                           input.catalog = input.catalog,
                           multi.types   = multi.types,
                           verbose       = verbose,
                           cos.merge     = cos.merge,
                           min.sample    = min.sample)

  # Step 3: Plot diagnostic plots, signatures, exposures
  # and compare with ground truth signature and exposures.

  if(!is.null(out.dir)) {
    AnalyzeAndPlotretval(retval,
                         out.dir,
                         ground.truth.sig,
                         ground.truth.exp,
                         verbose,
                         overwrite)
  }
  return(invisible(retval))
}
