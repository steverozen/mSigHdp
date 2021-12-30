#' Evaluate and plot retval from \code{CombinePosteriorChains} or \code{CombineChainsAndExtractSigs}
#' This function now works for both NR's pipeline and Mo's pipeline

#' @param retval the output from function  \code{CombinePosteriorChains} or \code{CombineChainsAndExtractSigs}
#'
#' @param out.dir Directory that will be created for the output;
#'   if \code{overwrite} is \code{FALSE} then
#'   abort if \code{out.dir} already exits.
#' @param overwrite If \code{TRUE} overwrite \code{out.dir} if it exists, otherwise raise an error.
#'
#' @param verbose If \code{TRUE} then \code{message} progress information.
#'
#' @param input.catalog input catalog matrix or
#'                            path to file with input catalog
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
                                 verbose          = TRUE,
                                 overwrite        = TRUE,
                                 diagnostic.plot  = TRUE) {

  if (dir.exists(out.dir)) {
    if (!overwrite) stop(out.dir, " already exists")
    if (verbose) message("Using existing out.dir ", out.dir)
  } else {
    dir.create(out.dir, recursive = T)
    if (verbose) message("Created new out.dir ", out.dir)
  }

  # Fragile, to be replaced by an ICAMS function
  IS.ICAMS <- any(grepl("Catalog", class(input.catalog)))

  if (verbose) message("Writing signatures")
  if(IS.ICAMS){
    extractedSignatures <- ICAMS::as.catalog(retval$signature,
                                             region       = "unknown",
                                             catalog.type = "counts.signature")
    ICAMS::WriteCatalog(extractedSignatures,
                        file.path(out.dir,"extracted.signatures.csv"))

    ICAMS::PlotCatalogToPdf(extractedSignatures,
                            file.path(out.dir, "extracted.signature.pdf"))

  }else{
    extractedSignatures <- retval$signature
    utils::write.csv(extractedSignatures,file.path(out.dir,"extracted.signatures.csv"),
                     row.names=F,quote=F)
  }



  if (verbose) message("Writing exposures")
  exposureCounts <- retval$exposureProbs %*% diag(colSums(input.catalog))
  colnames(exposureCounts) <- colnames(input.catalog)

  ICAMSxtra::WriteExposure(exposureCounts,
                           file.path(out.dir,"inferred.exposures.csv"))

  ICAMSxtra::PlotExposureToPdf(ICAMSxtra::SortExposure(exposureCounts),
                               file.path(out.dir,"inferred.exposure.count.pdf"))

  ICAMSxtra::PlotExposureToPdf(ICAMSxtra::SortExposure(retval$exposureProbs),
                               file.path(out.dir,"inferred.exposure.proportion.pdf"),
                               plot.proportion = TRUE)

    signature.post.samp.number <- retval$signature.post.samp.number

    signature.post.samp.number[,1] <- colnames(retval$signature)
    utils::write.csv(data.frame(signature.post.samp.number)
                     ,file = file.path(out.dir, "signature.post.samp.number.csv"),row.names = F,quote=F)

    low.confidence.signature <- retval$low.confidence.signature

    if(!is.null(ncol(low.confidence.signature)) &&
       (ncol(data.frame(low.confidence.signature))>0)) {

      low.confidence.signature <- apply(low.confidence.signature,2,function(x)x/sum(x))
      low.confidence.signature.post.samp.number <-
        retval$low.confidence.post.samp.number
      low.confidence.signature <- data.frame(low.confidence.signature)
      colnames(low.confidence.signature) <-
        paste0("low confidence hdp.",1:ncol(low.confidence.signature))
      low.confidence.signature.post.samp.number[,1] <- colnames(low.confidence.signature)
      row.names(low.confidence.signature) <- NULL

      if(IS.ICAMS){

        ICAMS::PlotCatalogToPdf(ICAMS::as.catalog(
          low.confidence.signature,
          infer.rownames = T,catalog.type = "counts.signature"),
          file.path(out.dir, "low.confidence.signatures.pdf"))
        ICAMS::WriteCatalog(ICAMS::as.catalog(
          low.confidence.signature,infer.rownames = T,catalog.type = "counts.signature"),
          file.path(out.dir,"low.confidence.signatures.csv"))
      }else{
        utils::write.csv(low.confidence.signature,
                         file.path(out.dir,"low.confidence.signatures.csv"),
                         row.names=F,quote=F)
      }


      utils::write.csv(data.frame(low.confidence.signature.post.samp.number),
                       file = file.path(out.dir,
                                        "low.confidence.signature.post.samp.number.csv"),
                       row.names = F,quote=F)

    }

  if(diagnostic.plot){

      dir.create(paste0(out.dir,"/Diagnostic_Plots"), recursive = T)
      ##this calls the diagnostic plotting function
      ##compatible with ML's component extraction

      mSigHdp::ComponentDiagnosticPlotting(retval        = retval,
                                           input.catalog = input.catalog,
                                           IS.ICAMS      = IS.ICAMS,
                                           out.dir       = paste0(out.dir,"/Diagnostic_Plots"),
                                           verbose       = verbose)
  }
}


