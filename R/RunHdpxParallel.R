#' Extract mutational signatures and optionally generate diagnostic plots to help
#' understand the results: e.g. the stability each extracted signature and
#' the tumors that drive the extraction of each signature.
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
#'                 rows are mutation types, columns are signatures with
#'                 high confidence.}
#'
#' \item{signature.post.samp.number}{A data frame with two columns. The first
#'                                   column corresponds to each signature in \code{signature}
#'                                   and the second columns contains the number of posterior
#'                                   samples that found the raw clusters contributing to the signature.}
#'
#' \item{signature.cdc}{A numeric data frame. Each column corresponds
#'                  to the sum of all mutations contributing to each signature in \code{signature}}
#'
#' \item{exposureProbs}{The inferred exposures as a matrix of mutation probabilities;
#'                      rows are signatures, columns are samples (e.g. tumors). This is
#'                      similar to \code{signature.cdc} but every column was normalized to sum of 1}
#'
#' \item{low.confidence.signature}{The profiles of signatures extracted with low confidence as a matrix; rows are mutation types,
#'                        columns are signatures with less than \code{high.confidence.prop} of posterior samples }
#'
#' \item{low.confidence.post.samp.number}{A data frame with two columns. The first column corresponds
#'                               to each signature in \code{low.confidence.signature} and the second
#'                               column contains the number of posterior samples that found
#'                               the raw clusters contributing to the signature.}
#'
#' \item{low.confidence.cdc}{A numeric data frame. Each column corresponds
#'                  to the sum of all mutations contributing to each
#'                  signature in \code{low.confidence.signature}}
#'
#' \item{extracted.retval}{A list object returned from code{\link[hdpx]{extract_components_from_clusters}}.}
#'
#' }
#'
#' @export

RunHdpxParallel <- function(input.catalog,
                            seedNumber           = 123,
                            K.guess,
                            multi.types          = TRUE,
                            verbose              = TRUE,
                            burnin               = 1000,
                            burnin.multiplier    = 10,
                            post.n               = 200,
                            post.space           = 100,
                            post.cpiter          = 3,
                            post.verbosity       = 0,
                            CPU.cores            = 20,
                            num.child.process    = 20,
                            high.confidence.prop = 0.9,
                            hc.cutoff           = 0.10,
                            overwrite           = TRUE,
                            out.dir             = NULL,
                            gamma.alpha         = 1,
                            gamma.beta          = 20,
                            gamma0.alpha        = gamma.alpha,
                            gamma0.beta         = gamma.beta,
                            checkpoint          = TRUE,
                            prior.sigs          = NULL,
                            prior.pseudoc       = NULL) {

  # Step 0: Get the input.catalog and keeping track of
  # whether it is an ICAMS catalog (encoded as an
  # additional class).
  input.catalog <- GetPossibleICAMSCatalog(input.catalog)

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
                              prior.sigs          = prior.sigs,
                              prior.pseudoc       = prior.pseudoc,
                              burnin.multiplier   = burnin.multiplier,
                              checkpoint          = checkpoint)

  # Step 2: Combine the posterior chains and extract
  # signatures and exposures;
  # retval has signatures, exposures, and multi.chains, a
  # hdpSampleMulti-class object.

  retval <-
    CombineChainsAndExtractSigs(chlist,
                                input.catalog  = input.catalog,
                                verbose        = verbose,
                                high.confidence.prop = high.confidence.prop,
                                hc.cutoff      = hc.cutoff)

  # Step 3: Save and plot signatures, exposures, diagnostics

  if(!is.null(out.dir)) {
    if (dir.exists(out.dir)) {
      if (!overwrite) stop(out.dir, " already exists")
      if (verbose) message("Using existing out.dir ", out.dir)
    } else {
      dir.create(out.dir, recursive = T)
      if (verbose) message("Created new out.dir ", out.dir)
    }
    save(retval, input.catalog, file = file.path(out.dir, "hdp.retval.Rdata"))
    AnalyzeAndPlotretval(retval                = retval,
                         input.catalog         = input.catalog,
                         out.dir               = out.dir,
                         verbose               = verbose,
                         overwrite             = overwrite)
  }
  return(invisible(retval))
}

