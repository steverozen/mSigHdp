#' Evaluate and plot retval from CombinePosteriorChains

#' @param retval the output from function CombinePosteriorChains
#'
#' @param out.dir Directory that will be created for the output;
#'   if \code{overwrite} is \code{FALSE} then
#'   abort if \code{out.dir} already exits.
#' @param overwrite If \code{TRUE} overwrite \code{out.dir} if it exists, otherwise
#'  raise an error.
#'
#' @param verbose If \code{TRUE} then \code{message} progress information.
#'
#' @param ground.truth.catalog Ground truth catalog matrix or
#'                            path to file with ground truth catalog
#'
#' @param ground.truth.exp Optional. Ground truth exposure matrix or
#'   path to file with ground truth exposures.
#'   If \code{NULL} skip checks that need this information.
#'
#' @param ground.truth.sig Optional. Either a string with the
#'   path to file with ground truth signatures or and
#'   \code{\link[ICAMS]{ICAMS}} catalog with the ground truth signatures.
#'   These are the signatures used to construct the ground truth spectra.
#'
#' @param diagnostic.plot If \code{TRUE} plot diagnostic plot.
#'    This is optional because there are cases having error
#'
#' @export
#'
#'
AnalyzeAndPlotretval <- function(retval,
                                 ground.truth.catalog,
                                 out.dir          = NULL,
                                 ground.truth.sig = NULL,
                                 ground.truth.exp = NULL,
                                 verbose          = TRUE,
                                 overwrite        = TRUE,
                                 diagnostic.plot  = TRUE) {

  if (dir.exists(out.dir)) {
    if (!overwrite) stop(out.dir, " already exits")
    if (verbose) message("Using existing out.dir ", out.dir)
  } else {
    dir.create(out.dir, recursive = T)
    if (verbose) message("Created new out.dir", out.dir)
  }

  save(retval, file = file.path(out.dir, "hdp.retval.Rdata"))

  # Plot the diagnostics of sampling chains.
  if(diagnostic.plot){


    if (mode(ground.truth.catalog) == "character") {
      if (verbose) message("Reading input catalog file ", ground.truth.catalog)
      ground.truth.catalog <- ICAMS::ReadCatalog(ground.truth.catalog, strict = FALSE)
    } else {
      ground.truth.catalog <- ground.truth.catalog
    }


    ChainsDiagnosticPlot(retval  = retval,
                         ground.truth.catalog = ground.truth.catalog,
                         out.dir = out.dir,
                         verbose = verbose)
  }


  if (verbose) message("Writing signatures")
  extractedSignatures <- ICAMS::as.catalog(retval$signature,
                                           region       = "unknown",
                                           catalog.type = "counts.signature")
  ICAMS::WriteCatalog(extractedSignatures,
                      file.path(out.dir,"extracted.signatures.csv"))

  ICAMS::PlotCatalogToPdf(extractedSignatures,
                          file.path(out.dir, "extracted.signature.pdf"))

  if (verbose) message("Writing exposures")
  ICAMSxtra::WriteExposure(retval$exposure,
                           file.path(out.dir,"inferred.exposures.csv"))

  ICAMSxtra::PlotExposureToPdf(ICAMSxtra::SortExposure(retval$exposure),
                               file.path(out.dir,"inferred.exposure.count.pdf"))

  ICAMSxtra::PlotExposureToPdf(ICAMSxtra::SortExposure(retval$exposure),
                               file.path(out.dir,"inferred.exposure.proportion.pdf"),
                               plot.proportion = TRUE)

  ###here is optional.


  # Do this early to catch any possible error before we do a lot
  # of computation
  # Exposure related plotting

  if(!is.null(ground.truth.exp)){
    ##read ground.truth.exp
    if (mode(ground.truth.exp) == "character") {
      if (verbose) {
        message("Reading ground truth exposures from ",
                ground.truth.exp)
      }
      ground.truth.exp <- ICAMSxtra::ReadExposure(ground.truth.exp)
    }
    #stopifnot(is.matrix(ground.truth.exp))

    ICAMSxtra::PlotExposureToPdf(ICAMSxtra::SortExposure(ground.truth.exp),
                                 file.path(out.dir,"ground.truth.exposure.count.pdf"))

    ICAMSxtra::PlotExposureToPdf(ICAMSxtra::SortExposure(ground.truth.exp),
                                 file.path(out.dir,"ground.truth.exposure.proportion.pdf"),
                                 plot.proportion = TRUE)
  }

  # Compare with ground truth
  # Only proceed if both ground.truth.sig and ground.truth.exp are provided

  if(!is.null(ground.truth.sig)){
    ##read ground.truth.exp
    if (mode(ground.truth.sig) == "character") {
      if (verbose) {
        message("Reading ground truth signatures from ",
                ground.truth.sig)
      }
      ground.truth.sig <- ICAMS::ReadCatalog(ground.truth.sig)
    }
    stopifnot(is.matrix(ground.truth.sig))

    sigAnalysis0 <- SynSigEval::MatchSigsAndRelabel(
      ex.sigs  = retval$signature,
      gt.sigs  = ground.truth.sig,
      exposure = ground.truth.exp)

    # Writes bi-directional matching and cos.sim calculation
    utils::write.csv(sigAnalysis0$match1, file = file.path(out.dir, "match1.w0.csv"))
    utils::write.csv(sigAnalysis0$match2, file = file.path(out.dir, "match2.w0.csv"))

    ICAMS::PlotCatalogToPdf(ICAMS::as.catalog(sigAnalysis0$gt.sigs,
                                              catalog.type = "counts.signature"), # Need to fix this
                            file.path(out.dir, "ground.truth.sigs.w0.pdf"))

    ICAMS::PlotCatalogToPdf(ICAMS::as.catalog(sigAnalysis0$ex.sigs,
                                              catalog.type = "counts.signature"),
                            file.path(out.dir, "extracted.sigs.w0.pdf"))



    utils::capture.output(
      cat("Call\n"),
      match.call(),
      cat("\nAverage cosine similarity\n"),
      sigAnalysis0$averCosSim,
      cat("\nAverage cosine similarity to each ground-truth signature\n"),
      sigAnalysis0$cosSim,
      cat("\nNumber of ground-truth signatures\n"),
      ncol(sigAnalysis0$gt.sigs),
      cat("\nNumber of extracted signatures\n"),
      ncol(sigAnalysis0$ex.sigs),
      cat("\nsigAnalysis0$extracted.with.no.best.match\n"),
      sigAnalysis0$extracted.with.no.best.match,
      cat("\nsigAnalysis0$ground.truth.with.no.best.match\n"),
      sigAnalysis0$ground.truth.with.no.best.match,
      file = file.path(out.dir,"other.results.txt"))
  }


}


