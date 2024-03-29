#' Deprecated, extract mutational signatures and optionally compare them to existing signatures and exposures.
#'
#' Deprecated, This functions uses the original method of combining raw clusters into "components".
#' Use \code{\link{RunHdpxParallel}} instead.
#'
#' @inheritParams MultipleSetupAndPosterior
#'
#' @inheritParams AnalyzeAndPlotretval
#'
#' @inheritParams CombinePosteriorChains
#'
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
#'
#' \item{sum_raw_clusters_after_cos_merge}{A matrix containing aggregated
#'       spectra of raw clusters after cosine
#'       similarity merge step in \code{\link[hdpx]{hdp_merge_and_extract_components}}.}
#'
#' \item{sum_raw_clusters_after_nonzero_categ}{A matrix containing aggregated
#'       spectra of raw clusters after non-zero category selecting
#'       step in \code{\link[hdpx]{hdp_merge_and_extract_components}}.}
#'
#' \item{clust_hdp0_ccc4}{A matrix containing aggregated spectra of
#'       raw clusters moving to hdp.0 after non-zero category selection step
#'       in \code{\link[hdpx]{hdp_merge_and_extract_components}}.}
#'
#' \item{clust_hdp0_ccc5}{A matrix containing aggregated spectra
#'       of raw clusters moving to hdp.0 after non-zero observation
#'       selection step in \code{\link[hdpx]{hdp_merge_and_extract_components}}.}
#'
#' }
#'
#' @export

OldRunHdpParallel <- function(input.catalog,
                           seedNumber          = 1,
                           K.guess,
                           multi.types         = FALSE,
                           verbose             = TRUE,
                           burnin         = 5000,
                           post.n              = 200,
                           post.space          = 100,
                           post.cpiter         = 3,
                           post.verbosity      = 0,
                           CPU.cores           = 20,
                           num.child.process   = 20,
                           cos.merge           = 0.9,
                           min.sample          = 1,
                           categ.CI            = 0.95,
                           exposure.CI         = 0.95,
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
                           burnin.multiplier   = 2,
                           burnin.checkpoint   = TRUE){

  # warning("This function is deprecated; use RunHdp xParallel instead.")

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
                              burnin.checkpoint   = burnin.checkpoint)

  # Step 2: Combine the posterior chains and extract
  # signatures and exposures;
  # multi.chains.etc has signatures, exposures, and multi.chains, a
  # hdpSampleMulti-class object.

  multi.chains.etc <-
    CombinePosteriorChains(chlist,
                           input.catalog  = input.catalog,
                           multi.types    = multi.types,
                           verbose        = verbose,
                           cos.merge      = cos.merge,
                           min.sample     = min.sample,
                           categ.CI       = categ.CI,
                           exposure.CI    = exposure.CI)

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
