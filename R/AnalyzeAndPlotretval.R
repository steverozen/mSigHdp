#' Evaluate and plot retval from CombinePosteriorChains
#' This function now calls for NR's pipeline or Mo's pipeline

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
#' @param input.catalog input catalog matrix or
#'                            path to file with input catalog
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
                                 input.catalog,
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
    if (verbose) message("Created new out.dir ", out.dir)
  }

  save(retval, file = file.path(out.dir, "hdp.retval.Rdata"))

  if (mode(input.catalog) == "character") {
    if (verbose) message("Reading input catalog file ", input.catalog)
    input.catalog <- ICAMS::ReadCatalog(input.catalog, strict = FALSE)
  } else {
    input.catalog <- input.catalog
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

  if(!("extracted.retval" %in% names(retval))){

    # Plot the diagnostics of sampling chains.
    if(diagnostic.plot){

      dir.create(paste0(out.dir,"/Diagnostic_Plots"), recursive = T)

      mSigHdp::ChainsDiagnosticPlot(retval  = retval,
                                    input.catalog = input.catalog,
                                    out.dir = paste0(out.dir,"/Diagnostic_Plots"),
                                    verbose = verbose)
    }



  }else{
    # Plot the diagnostics of sampling chains.
    if(diagnostic.plot){

      dir.create(paste0(out.dir,"/Diagnostic_Plots"), recursive = T)

      mSigHdp::ChainsDiagnosticPlotMo(retval  = retval,
                                      input.catalog = input.catalog,
                                      out.dir = paste0(out.dir,"/Diagnostic_Plots"),
                                      verbose = verbose)
    }


    moderate.spectrum <- retval$extracted.retval$moderate.spectrum
    noise.spectrum <- retval$extracted.retval$noise.spectrum

    extracted.stats <- retval$extracted.retval$high.confident.stats
    moderate.stats <- retval$extracted.retval$moderate.stats
    noise.stats <- retval$extracted.retval$noise.stats

    utils::write.csv(t(extracted.stats),file = file.path(out.dir, "high.confidence.sig.stats.csv"))


    if(!is.null(ncol(moderate.spectrum))){
      moderate.spectrum <- moderate.spectrum[,order(moderate.stats,decreasing=T)]
      moderate.stats <- moderate.stats[order(moderate.stats,decreasing=T)]

      colnames(moderate.spectrum) <- paste0("potential hdp.",1:ncol(moderate.spectrum))
      row.names(moderate.spectrum) <- row.names(extractedSignatures)



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

        sigAnalysis.moderate <- ICAMSxtra::MatchSigsAndRelabel(
          ex.sigs  = moderate.spectrum,
          gt.sigs  = ground.truth.sig,
          exposure = ground.truth.exp)

        moderate.catalog <- ICAMS::as.catalog(sigAnalysis.moderate$ex.sigs)
        moderate.catalog <- ICAMS::TransformCatalog(moderate.catalog,target.catalog.type = "counts.signature")


        ICAMS::PlotCatalogToPdf(moderate.catalog,
                                file.path(out.dir, "moderate.sigs.w0.pdf"))

        utils::write.csv(t(moderate.stats),file = file.path(out.dir, "moderate.sig.stats.csv"))
      }
      if(!is.null(ncol(noise.spectrum))){

        noise.spectrum <- noise.spectrum[,order(noise.stats,decreasing=T)]
        noise.stats <- noise.stats[order(noise.stats,decreasing=T)]

        noise.spectrum <- data.frame(noise.spectrum)
        colnames(noise.spectrum) <- paste0("noise hdp.",1:ncol(noise.spectrum))
        row.names(noise.spectrum) <- row.names(extractedSignatures)

        ICAMS::PlotCatalogToPdf(ICAMS::as.catalog(noise.spectrum),
                                file.path(out.dir, "clusters.with.low.confidence.pdf"))
        utils::write.csv(t(noise.stats),file = file.path(out.dir, "noise.stats.csv"))


      }

    }

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
  }
  ###here is optional.


  # Do this early to catch any possible error before we do a lot
  # of computation
  # Exposure related plotting
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

    sigAnalysis0 <- ICAMSxtra::MatchSigsAndRelabel(
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


