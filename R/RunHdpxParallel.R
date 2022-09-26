#' Extract (discover) mutational signatures from a matrix of mutational spectra
#'
#' Please see the vignette for an example.
#'
#' @inheritParams ParallelGibbsSample
#'
#' @inheritParams SaveAnalysis
#'
#' @inheritParams CombineChainsAndExtractSigs
#'
#' @param hc.cutoff Deprecated, use \code{merge.raw.cluster.args}.
#'
#' @param downsample_threshold See \code{\link{downsample_spectra}}
#'        and \code{link{show_downsample_curves}}.
#'
#' @inherit CombineChainsAndExtractSigs return
#'
#' @details Please see our paper at
#' https://www.biorxiv.org/content/10.1101/2022.01.31.478587v1
#' for suggestions on argument values to use.
#'
#' @export

RunHdpxParallel <- function(input.catalog,
                            seedNumber             = 123,
                            K.guess,
                            multi.types            = FALSE,
                            verbose                = FALSE,
                            burnin                 = 1000,
                            burnin.multiplier      = 10,
                            post.n                 = 200,
                            post.space             = 100,
                            post.cpiter            = 3,
                            post.verbosity         = 0,
                            CPU.cores              = 20,
                            num.child.process      = 20,
                            high.confidence.prop   = 0.9,
                            hc.cutoff              = NULL,
                            merge.raw.cluster.args =
                              hdpx::default_merge_raw_cluster_args(),
                            overwrite              = TRUE,
                            out.dir                = NULL,
                            gamma.alpha            = 1,
                            gamma.beta             = 20,
                            checkpoint             = TRUE,
                            downsample_threshold   = NULL) {

  # Check for suitable version of hdpx
  if (utils::packageVersion("hdpx") < "1.0.3.0009") {
    stop("hdpx version must be >= 1.0.3.0009")
  }

  if (!is.null(hc.cutoff)) {
    merge.raw.cluster.args <- hdpx::default_merge_raw_cluster_args()
    merge.raw.cluster.args$clustering.cutoff <- 1 - hc.cutoff
    warning("hc.cutoff is deprecated, use merge.raw.cluster.args")
  }
  rm(hc.cutoff)

  # Step 0: Get the input.catalog and keeping track of
  # whether it is an ICAMS catalog (encoded as an
  # additional class).
  input.catalog <- GetPossibleICAMSCatalog(input.catalog)

  if (!is.null(downsample_threshold)) {
    tmp.catalog <-
      downsample_spectra(input.catalog,
                         downsample_threshold = downsample_threshold)$down_spec
    if (ICAMS::IsICAMSCatalog(input.catalog)) {
      tmp.catalog <-
        ICAMS::as.catalog(tmp.catalog,
                          ref.genome   = attr(input.catalog, "ref.genome"),
                          region       = attr(input.catalog, "region"),
                          catalog.type = attr(input.catalog, "catalog.type"),
                          abundance    = attr(input.catalog, "abundance"))
    }
    input.catalog <- tmp.catalog
    rm(tmp.catalog)
  }

  # Step 1: Activate hierarchical Dirichlet processes and
  # run posterior sampling in parallel;
  # chlist is a list of hdpSampleChain-class objects.

  chlist <-
    ParallelGibbsSample(input.catalog,
                        seedNumber            = seedNumber,
                        K.guess               = K.guess,
                        multi.types           = multi.types,
                        verbose               = verbose,
                        burnin                = burnin,
                        post.n                = post.n,
                        post.space            = post.space,
                        post.cpiter           = post.cpiter,
                        post.verbosity        = post.verbosity,
                        CPU.cores             = CPU.cores,
                        num.child.process     = num.child.process,
                        gamma.alpha           = gamma.alpha,
                        gamma.beta            = gamma.beta,
                        burnin.multiplier     = burnin.multiplier,
                        checkpoint            = checkpoint)

  # For preparing test data
  if (FALSE) {
    save(chlist, file = "big.chlist.from.ParallelGibbsSample.Rdata")
  }

  # Step 2: Combine the posterior chains and extract
  # signatures and exposures;
  # retval has signatures, exposures, and multi.chains, a
  # hdpSampleMulti-class object.

  retval <-
    CombineChainsAndExtractSigs(chlist,
                                input.catalog          = input.catalog,
                                verbose                = verbose,
                                high.confidence.prop   = high.confidence.prop,
                                merge.raw.cluster.args = merge.raw.cluster.args)

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
    SaveAnalysis(retval                = retval,
                 input.catalog         = input.catalog,
                 out.dir               = out.dir,
                 verbose               = verbose,
                 overwrite             = overwrite)
  }
  return(invisible(retval))
}

