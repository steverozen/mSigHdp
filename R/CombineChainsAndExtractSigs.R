#' Extract components and exposures from multiple posterior sample chains
#' This function returns signatures with high confidence (found in more than 90% #' posterior samples)
#'
#' @param clean.chlist A list of \code{\link[hdpx]{hdpSampleChain-class}} objects.
#'  Each element is the result of one posterior sample chain.
#'
#' @param input.catalog Input spectra catalog as a matrix or
#' in \code{\link[ICAMS]{ICAMS}} format.
#'
#' @param multi.types A logical scalar or
#' a character vector.
#' If \code{FALSE}, The HDP analysis
#'   will regard all input spectra as one tumor type.
#'   HDP structure as one parent node for all tumors
#'
#' If \code{TRUE}, the HDP analysis
#'   will infer tumor types based on the string before "::" in their names.
#'   e.g. tumor type for "SA.Syn.Ovary-AdenoCA::S.500" would be "SA.Syn.Ovary-AdenoCA"
#'   HDP structure as a grandparent node for whole data and
#'   one parent node for each tumor type
#'
#' If \code{multi.types} is a character vector, then it should be of the same length
#' as the number of columns in \code{input.catalog}, and each value is the
#' name of the tumor type of the corresponding column in \code{input.catalog}.
#'
#' e.g. \code{c("SA.Syn.Ovary-AdenoCA", "SA.Syn.Kidney-RCC")}.
#'
#' @param verbose If \code{TRUE} then \code{message} progress information.
#'
#' @param cos.merge The cosine similarity threshold for merging raw clusters
#'      from the posterior sampling chains into "components" i.e. signatures;
#'      passed to \code{\link[hdpx]{extract_components_from_clusters}}.
#'
#' @param confident.prop Pass to \code{\link[hdpx]{interpret_components}}.
#'                       clusters with at least \code{confident.prop} of posterior
#'                       samples are high confident signatures
#' @param noise.prop Pass to \code{\link[hdpx]{interpret_components}}.
#'                   Clusters with less than \code{noise.prop} of posterior samples
#'                  are noise signatures
#' @param hc.cutoff Pass to \code{\link[hdpx]{extract_components_from_clusters}}. The cutoff of
#'                  height of hierarchical clustering dendrogram
#' @param hc.method Pass to \code{\link[hdpx]{extract_components_from_clusters}}. The agglomeration method
#'
#' @param hc 'agglomerative' or 'divisive'.
#'            Pass to \code{\link[hdpx]{extract_components_from_clusters}} For testing purpose now.
#'                  for hierarchical clustering. Default is "ward.D2".
#' @return Invisibly, a list with the following elements:\describe{
#' \item{signature}{The extracted signature profiles as a matrix;
#'             rows are mutation types, columns are signatures including high confident signatures -'hdp' signatures
#'            and moderate confident signatures - 'potential hdp' signatures.}
#'
#' \item{signature.post.samp.number}{A dataframe with two columns. The first column corresponds to each signature
#'                                   in \code{signature} and the second columns contains the number of posterior
#'                                   samples that found the raw clusters contributing to the signature.}
#'
#' \item{signature.cdc}{A \code{\link[hdpx]{comp_dp_counts}} like dataframe. Each column corresponds to the sum of
#'                      all \code{\link[hdpx]{comp_dp_counts}} matrices of the raw clusters contributing to each
#'                      signature in code{signature}}
#'
#' \item{exposureProbs}{The inferred exposures as a matrix of mutation probabilities;
#'            rows are signatures, columns are samples (e.g. tumors).}
#'
#' \item{noise.signature}{The extracted signature profiles as a matrix; rows are mutation types,
#'                        columns are signatures with less than \code{noise.prop} of posterior samples.}
#'
#' \item{noise.post.samp.number}{A data frame with two columns. The first column corresponds to each signature
#'                                   in \code{noise.signature} and the second columns contains the number of posterior
#'                                   samples that found the raw clusters contributing to the signature.}
#'
#' \item{noise.cdc}{A \code{\link[hdpx]{comp_dp_counts}} like data frame. Each column corresponds to the sum of
#'                      all \code{\link[hdpx]{comp_dp_counts}} matrices of the raw clusters contributing to each
#'                      signature in code{noise.signature}}
#'
#' \item{extracted.retval}{A list object returned from code{\link[hdpx]{interpret_components}}.}
#' }
#'
#' @export
#'
CombineChainsAndExtractSigs <-
  function(clean.chlist,
           input.catalog,
           multi.types,
           verbose             = TRUE,
           cos.merge           = 0.9,
           confident.prop      = 0.9,
           noise.prop          = 0.1,
           hc.cutoff           = 0.10,
           hc.method           = "ward.D2",
           hc                  = "agglomerative"
  ) {
    if (mode(input.catalog) == "character") {
      if (verbose) message("Reading input catalog file ", input.catalog)
      input.catalog <- ICAMS::ReadCatalog(input.catalog)
    } else {
      input.catalog <- input.catalog
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
                                               cos.merge = cos.merge,
                                               hc.cutoff = hc.cutoff,
                                               hc.method = hc.method,
                                               hc        = hc
        )
    )

    if (verbose) {
      message("extract_sigs_from_clusters time: ")
      for (xn in names(extract.time)) {
        message(" ", xn, " ", extract.time[[xn]])
      }
    }

    intepret.comp.retval <-  hdpx::interpret_components(multi.chains.retval = multi.chains.retval,
                                                       confident.prop      = confident.prop,
                                                       noise.prop          = noise.prop,
                                                       verbose            = verbose)


    confidentSignatures <- data.frame(intepret.comp.retval$high_confident_components)

    rownames(confidentSignatures) <- rownames(input.catalog)
    # Set signature names to "hdp.0","hdp.1","hdp.2", ...
    colnames(confidentSignatures) <-
      paste("hdp", c(1:ncol(confidentSignatures)), sep = ".")

    combinedSignatures <- confidentSignatures

    potentialSignatures <- data.frame(intepret.comp.retval$moderate_components)

    if(!is.null(potentialSignatures) && ncol(potentialSignatures)>0){
      potentialSignatures <- apply(potentialSignatures,2,function(x)x/sum(x))

      rownames(potentialSignatures) <- rownames(input.catalog)
      # Set signature names to "potential hdp.0","potential hdp.1","potential hdp.2", ...
      colnames(potentialSignatures) <-
        paste("potential hdp", c(1:ncol(potentialSignatures)), sep = ".")
      combinedSignatures <- cbind(combinedSignatures,potentialSignatures)

    }

    combinedSignatures <- apply(combinedSignatures,2,function(x)x/sum(x))
    combined.stats <- rbind(intepret.comp.retval$high_confident_components_post_number,
                            intepret.comp.retval$moderate_components_post_number)
    combined.cdc  <- cbind(intepret.comp.retval$high_confident_components_cdc,
                           intepret.comp.retval$moderate_components_cdc)

    if (verbose) message("extracting signatures exposures ", Sys.time())

    exposureProbs <- t(apply(combined.cdc,1,function(x){x/sum(x)}))
    if(nrow(exposureProbs)==1){
      exposureProbs <- t(exposureProbs)
    }
    # Remove columns corresponding to parent or grandparent nodes
    # (leaving only columns corresponding to samples.
    # Transpose so it conforms to SynSigEval format
    exposureProbs <- t(exposureProbs[-c(1:(nrow(exposureProbs)-ncol(input.catalog))), ])

    colnames(exposureProbs) <- colnames(input.catalog)

    row.names(exposureProbs) <- colnames(combinedSignatures)

    return(invisible(list(signature                   = combinedSignatures,
                          signature.post.samp.number  = combined.stats,
                          signature.cdc               = combined.cdc,
                          exposureProbs               = exposureProbs,
                          noise.signature             = intepret.comp.retval$noise_components,
                          noise.post.samp.number      = intepret.comp.retval$noise_components_post_number,
                          noise.cdc                   = intepret.comp.retval$noise_components_cdc,
                          extracted.retval            = multi.chains.retval)))

  }
