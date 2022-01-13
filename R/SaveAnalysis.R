#' Save the analysis from \code{\link{CombineChainsAndExtractSigs}} to files.
#'
#' @param retval The output from \code{CombineChainsAndExtractSigs}.
#'
#' @param out.dir Directory that will be created for the output, including
#'   csv files and plots (pdfs) of extracted signatures and their exposures.
#'
#' @param overwrite If \code{TRUE} overwrite \code{out.dir}
#'    if it exists, otherwise raise an error.
#'
#' @param verbose If \code{TRUE} then \code{message} progress information.
#'
#' @param input.catalog input catalog matrix or
#'                            path to file with input catalog
#'
#' @param diagnostic.plot If \code{TRUE} plot diagnostic plots
#'   in a subdirectory \code{Diagnostic_Plots} of \code{out.dir}.
#'
#' @export
#'
#'
SaveAnalysis <- function(retval,
                         input.catalog,
                         out.dir          = NULL,
                         verbose          = TRUE,
                         overwrite        = TRUE,
                         diagnostic.plot  = TRUE) {

  IS.ICAMS <- ICAMS::IsICAMSCatalog(input.catalog)
  if (dir.exists(out.dir)) {
    if (!overwrite) stop(out.dir, " already exists")
    if (verbose) message("Using existing out.dir ", out.dir)
  } else {
    dir.create(out.dir, recursive = T)
    if (verbose) message("Created new out.dir ", out.dir)
  }

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

  mSigAct::WriteExposure(exposureCounts,
                           file.path(out.dir,"inferred.exposures.csv"))

  mSigAct::PlotExposureToPdf(mSigAct::SortExposure(exposureCounts),
                               file.path(out.dir,"inferred.exposure.count.pdf"))

  mSigAct::PlotExposureToPdf(mSigAct::SortExposure(retval$exposureProbs),
                               file.path(out.dir,"inferred.exposure.proportion.pdf"),
                               plot.proportion = TRUE)

    signature.post.samp.number <- retval$signature.post.samp.number

    signature.post.samp.number[,1] <- colnames(retval$signature)

    signature.post.samp.number <- data.frame(signature.post.samp.number)
    colnames(signature.post.samp.number) <- c("Signature","NumberOfPostSamples")
    utils::write.csv(
      signature.post.samp.number,
      file = file.path(out.dir,
                       "extracted.signatures.post.samp.number.csv"),
      row.names = F,quote=F)

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

      low.confidence.signature.post.samp.number <- data.frame(low.confidence.signature.post.samp.number)
      colnames(low.confidence.signature.post.samp.number) <- c("Signature","NumberOfPostSamples")
      utils::write.csv(low.confidence.signature.post.samp.number,
                       file = file.path(out.dir,
                                        "low.confidence.signatures.post.samp.number.csv"),
                       row.names = F,quote=F)

    }

  if (diagnostic.plot) {

    dir.create(paste0(out.dir,"/Diagnostic_Plots"), recursive = T)
    ##this calls the diagnostic plotting function
    ##compatible with ML's component extraction

    ComponentDiagnosticPlotting(retval        = retval,
                                input.catalog = input.catalog,
                                out.dir       = paste0(out.dir,"/Diagnostic_Plots"),
                                verbose       = verbose)
  }

  return(NULL)
}


