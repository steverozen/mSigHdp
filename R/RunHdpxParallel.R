#' Extract mutational signatures and optionally compare them to existing signatures and exposures.
#'
#' @inheritParams MultipleSetupAndPosterior
#'
#' @inheritParams AnalyzeAndPlotretval
#'
#' @inheritParams CombineChainsAndExtractSigs
#'
#'
#' @return Invisibly, a list with the following elements:\describe{
#' \item{signature}{The extracted signature profiles as a matrix;
#'                 rows are mutation types, columns are signatures including
#'                 high confident signatures -'hdp' signatures and moderate
#'                 confident signatures - 'potential hdp' signatures.}
#'
#' \item{signature.post.samp.number}{A dataframe with two columns. The first
#'                                   column corresponds to each signature in \code{signature}
#'                                   and the second columns contains the number of posterior
#'                                   samples that found the raw clusters contributing to the signature.}
#'
#' \item{signature.cdc}{A \code{\link[hdpx]{comp_dp_counts}} like dataframe.
#'                      Each column corresponds to the sum of all \code{\link[hdpx]{comp_dp_counts}}
#'                      matrices of the raw clusters contributing to each signature in code{signature}}
#'
#' \item{exposureProbs}{The inferred exposures as a matrix of mutation probabilities;
#'                      rows are signatures, columns are samples (e.g. tumors).}
#'
#' \item{noise.signature}{The extracted signature profiles as a matrix; rows are mutation types,
#'                        columns are signatures with less than \code{noise.prop} of posterior samples }
#'
#' \item{noise.post.samp.number}{A data frame with two columns. The first column corresponds
#'                               to each signature in \code{noise.signature} and the second
#'                               column contains the number of posterior samples that found
#'                               the raw clusters contributing to the signature.}
#'
#' \item{noise.cdc}{A \code{\link[hdpx]{comp_dp_counts}} like data frame. Each column corresponds
#'                  to the sum of all \code{\link[hdpx]{comp_dp_counts}} matrices of the raw clusters
#'                  contributing to each signature in code{noise.signature}}
#'
#' \item{extracted.retval}{A list object returned from code{\link[hdpx]{interpret_components}}.}
#'
#' }
#'
#' @export

RunHdpxParallel <- function(input.catalog,
                           seedNumber          = 123,
                           K.guess,
                           multi.types         = FALSE,
                           verbose             = TRUE,
                           burnin              = 5000,
                           burnin.multiplier   = 2,
                           burnin.checkpoint   = FALSE,
                           post.n              = 200,
                           post.space          = 100,
                           post.cpiter         = 3,
                           post.verbosity      = 0,
                           CPU.cores           = 20,
                           num.child.process   = 20,
                           confident.prop      = 0.9,
                           noise.prop          = 0.5,
                           hc.cutoff           = 0.10,
                           hc.method           = "average",
                           hc                  = "agglomerative",
                           ground.truth.sig    = NULL,
                           ground.truth.exp    = NULL,
                           overwrite           = TRUE,
                           out.dir             = NULL,
                           gamma.alpha         = 1,
                           gamma.beta          = 20,
                           gamma0.alpha        = gamma.alpha,
                           gamma0.beta         = gamma.beta,
                           checkpoint.chlist   = TRUE,
                           checkpoint.1.chain  = TRUE,
                           prior.sigs          = NULL,
                           prior.pseudoc       = NULL,
                           posterior.checkpoint= FALSE){

  # Step 1: Activate hierarchical Dirichlet processes and
  # run posterior sampling in parallel;
  # chlist is a list of hdpSampleChain-class objects.

  chlist <-
    MultipleSetupAndPosterior(input.catalog,
                              seedNumber          = seedNumber,
                              K.guess             = K.guess,
                              multi.types         = multi.types,
                              verbose             = verbose,
                              burnin              = burnin,
                              post.n              = post.n,
                              post.space          = post.space,
                              post.cpiter         = post.cpiter,
                              post.verbosity      = post.verbosity,
                              CPU.cores           = CPU.cores,
                              num.child.process   = num.child.process,
                              gamma.alpha         = gamma.alpha,
                              gamma.beta          = gamma.beta,
                              gamma0.alpha        = gamma0.alpha,
                              gamma0.beta         = gamma0.beta,
                              checkpoint.chlist   = checkpoint.chlist,
                              checkpoint.1.chain  = checkpoint.1.chain,
                              prior.sigs          = prior.sigs,
                              prior.pseudoc       = prior.pseudoc,
                              burnin.multiplier   = burnin.multiplier,
                              burnin.checkpoint   = burnin.checkpoint,
                              posterior.checkpoint= posterior.checkpoint)

  # Step 2: Combine the posterior chains and extract
  # signatures and exposures;
  # multi.chains.etc has signatures, exposures, and multi.chains, a
  # hdpSampleMulti-class object.

  multi.chains.etc <-
    CombineChainsAndExtractSigs(chlist,
                           input.catalog  = input.catalog,
                           multi.types    = multi.types,
                           verbose        = verbose,
                           confident.prop = confident.prop,
                           noise.prop     = noise.prop,
                           hc.cutoff      = hc.cutoff,
                           hc.method      = hc.method,
                           hc             = hc)

  # Step 3: Plot diagnostic plots, signatures, exposures
  # and compare with ground truth signature and exposures.

  if(!is.null(out.dir)) {

    AnalyzeAndPlotretval(retval                = multi.chains.etc,
                         input.catalog         = input.catalog,
                         out.dir               = out.dir,
                         ground.truth.sig      = ground.truth.sig,
                         ground.truth.exp      = ground.truth.exp,
                         verbose               = verbose,
                         overwrite             = overwrite)
  }
  return(invisible(multi.chains.etc))
}

