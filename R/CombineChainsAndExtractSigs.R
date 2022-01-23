#' Extract signatures etc. from multiple Gibbs sample chains
#
#' @param clean.chlist A list of \code{\link[hdpx]{hdpSampleChain-class}}
#'  objects (from package hdpx), typically returned from
#'  \code{ParallelGibbsSample}.
#'  Each element must be the result of one posterior sample chain.
#'
#'
#' @param input.catalog Input spectra catalog as a matrix.
#'
#' @param verbose If \code{TRUE} then \code{message} progress information.
#'
#' @param high.confidence.prop
#'                       Raw clusters of mutations
#'                       found in \eqn{>=} \code{high.confidence.prop} proportion of posterior
#'                       samples are signatures with high confidence.
#'
#                       (Passed to \code{\link[hdpx]{interpret_components}}).
#'
#' @param hc.cutoff The cutoff of
#'                  height of the hierarchical clustering dendrogram used in combining
#'                  raw clusters of mutations into aggregated clusters.
#'
#                  (Passed to \code{hdpx::\link[hdpx]{extract_components_from_clusters}}.)
#'
#' @return Invisibly, a list with the following elements:\describe{
#'\item{signature}{The extracted signature profiles as a matrix;
#'                 rows are mutation types, columns are signatures with
#'                 high confidence.}
#'
#' \item{signature.post.samp.number}{A data frame with two columns. The first
#'            column corresponds to each signature in \code{signature}
#'            and the second columns contains the number of posterior
#'            samples that found the raw clusters contributing to the signature.}
#'
#' \item{signature.cdc}{
#'      A numeric data frame. Columns correspond to signatures
#'      as in \code{signature}. Rows correspond to either biological
#'      samples or to parent and grandparent Dirichlet processes.}
#'
#' \item{exposureProbs}{The inferred exposures as a matrix
#'        of mutation probabilities;
#'        rows are signatures, columns are samples (e.g. tumors). This is
#'        similar to \code{signature.cdc}, but every column was normalized
#'        to sum to 1.}
#'
#' \item{low.confidence.signature}{The profiles of signatures extracted
#'     with low confidence as a matrix; rows are mutation types,
#'     columns are signatures with <
#'     \code{high.confidence.prop} of posterior samples.}
#'
#' \item{low.confidence.post.samp.number}{
#'      Analogous to \code{signature.post.samp.number}, except that
#'      this one is for
#'      signatures in \code{low.confidence.signature}.}
#'
#' \item{low.confidence.cdc}{Analogous to
#'      \code{signature.cdc}, except that columns in this
#'      matrix correspond to columns
#'      in \code{low.confidence.signature}.}
#'
#' \item{extracted.retval}{A list object returned from
#'          \code{\link[hdpx]{extract_components_from_clusters}}
#'          in package hdpx.}
#'
#' }
#'
#' @export
#'
CombineChainsAndExtractSigs <-
  function(clean.chlist,
           input.catalog,
           verbose              = FALSE,
           high.confidence.prop = 0.9,
           hc.cutoff            = 0.10
  ) {

    input.catalog <- GetPossibleICAMSCatalog(input.catalog)
    if (FALSE) {
      if (mode(input.catalog) == "character") {
        if (verbose) message("Reading input catalog file ", input.catalog)
        input.catalog <- ICAMS::ReadCatalog(input.catalog)
      } else {
        input.catalog <- input.catalog
      }
    }
    input.catalog <- input.catalog[,colSums(input.catalog)>0]
    convSpectra <- t(input.catalog)
    number.channels <- nrow(input.catalog)
    number.samples  <- ncol(input.catalog)

    multi.chains <- hdpx::hdp_multi_chain(clean.chlist)
    if (verbose) message("calling extract_components_from_clusters ", Sys.time())
    # Group raw "clusters" into "components" (i.e. signatures).

    extract.time <- system.time(
      multi.chains.retval <-
        hdpx::extract_components_from_clusters(multi.chains,
                                               hc.cutoff = hc.cutoff)
    )

    if (verbose) {
      message("extract_sigs_from_clusters time: ")
      for (xn in names(extract.time)) {
        message(" ", xn, " ", extract.time[[xn]])
      }
    }

    intepret.comp.retval <-
      hdpx::interpret_components(multi.chains.retval  = multi.chains.retval,
                                 high.confidence.prop = high.confidence.prop,
                                 verbose              = verbose)

    confidentSignatures <-
      data.frame(intepret.comp.retval$high_confidence_components)

    rownames(confidentSignatures) <- rownames(input.catalog)
    # Set signature names to "hdp.0","hdp.1","hdp.2", ...
    colnames(confidentSignatures) <-
      paste("hdp", c(1:ncol(confidentSignatures)), sep = ".")

    combinedSignatures <- confidentSignatures

    combinedSignatures <- apply(combinedSignatures,2,function(x)x/sum(x))
    combined.stats <- intepret.comp.retval$high_confidence_components_post_number

    combined.cdc  <- intepret.comp.retval$high_confidence_components_cdc

    if (verbose) message("extracting signatures exposures ", Sys.time())

    exposureProbs <- t(apply(combined.cdc,1,function(x){x/sum(x)}))
    if(nrow(exposureProbs)==1){
      exposureProbs <- t(exposureProbs)
    }
    # Remove columns corresponding to parent or grandparent nodes
    # (leaving only columns corresponding to samples.)
    # Transpose so it conforms to SynSigEval format
    exposureProbs <- t(exposureProbs[-c(1:(nrow(exposureProbs)-ncol(input.catalog))), ])

    colnames(exposureProbs) <- colnames(input.catalog)

    row.names(exposureProbs) <-
      colnames(combined.cdc) <-
      combined.stats[,1] <-
      colnames(combinedSignatures) # These are the

    colnames(combined.stats) <- c("Signature","NumberOfPostSamples")

    low.confidence.signature <-
      data.frame(intepret.comp.retval$low_confidence_components)

    low.confidence.post.samp.number <-
      data.frame(intepret.comp.retval$low_confidence_components_post_number)

    # browser()

    if(!is.null(ncol(low.confidence.signature)) &&
       (ncol(data.frame(low.confidence.signature))>0)) {

      if (FALSE) {
      low.confidence.cdc <-
        data.frame(intepret.comp.retval$low_confidence_components_cdc[,1:ncol(low.confidence.signature)])
      } else {
        low.confidence.cdc <-
          data.frame(intepret.comp.retval$low_confidence_components_cdc)
      }
      # browser() # The next line is hack for testing -- previously these rownames
      # were sometimes integer and sometimes character.
      rownames(low.confidence.cdc) <- as.character(rownames(low.confidence.cdc))

      colnames(low.confidence.cdc) <-
        low.confidence.post.samp.number[,1] <-
        colnames(low.confidence.signature) <-
        paste("low confidence hdp", 1:ncol(low.confidence.signature), sep = ".")

      colnames(low.confidence.post.samp.number) <-
        c("Signature","NumberOfPostSamples")

    } else {
      low.confidence.signature <-
        low.confidence.post.samp.number <-
        low.confidence.cdc <- NULL
    }

    return(invisible(list(signature                   = combinedSignatures,
                          signature.post.samp.number  = combined.stats,
                          signature.cdc               = combined.cdc,
                          exposureProbs               = exposureProbs,
                          low.confidence.signature             = low.confidence.signature,
                          low.confidence.post.samp.number      = low.confidence.post.samp.number,
                          low.confidence.cdc                   = low.confidence.cdc,
                          extracted.retval            = multi.chains.retval)))

  }
