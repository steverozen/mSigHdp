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
#'
#' If \code{TRUE}, the HDP analysis
#'   will infer tumor types based on the string before "::" in their names.
#' e.g. tumor type for "SA.Syn.Ovary-AdenoCA::S.500" would be "SA.Syn.Ovary-AdenoCA"
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
#'      passed to \code{\link[hdpx]{extract_sigs_from_clusters}}.
#'
#' @param confident.prop clusters with at least \code{confident.prop} of posterior samples are high confident signatures
#' @param noise.prop clusters with less than \code{noise.prop} of posterior samples are noise signatures
#' @param hc.cutoff passed to \code{\link[hdpx]{extract_sigs_from_clusters}}. The cutoff of height of hierarchical clustering                          dendrogram(default is 0.12)
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
#'
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
           hc.cutoff           = 0.12
  ) {
    if (mode(input.catalog) == "character") {
      if (verbose) message("Reading input catalog file ", input.catalog)
      input.catalog <- ICAMS::ReadCatalog(input.catalog, strict = FALSE)
    } else {
      input.catalog <- input.catalog
    }

    convSpectra <- t(input.catalog)
    number.channels <- nrow(input.catalog)
    number.samples  <- ncol(input.catalog)

    multi.chains <- hdpx::hdp_multi_chain(clean.chlist)
    if (verbose) message("calling extract_sigs_from_clusters ", Sys.time())
    # Group raw "clusters" into "components" (i.e. signatures).

    extract.time <- system.time(
      multi.chains.retval <-
        hdpx::extract_sigs_from_clusters(multi.chains,
                                         cos.merge      = cos.merge,
                                         hc.cutoff = hc.cutoff
        )
    )


    if (verbose) {
      message("extract_sigs_from_clusters time: ")
      for (xn in names(extract.time)) {
        message(" ", xn, " ", extract.time[[xn]])
      }
    }
    if (verbose) message("extracting signatures ", Sys.time())
    spectrum.df <- multi.chains.retval$clustered.spectrum
    spectrum.stats <- multi.chains.retval$stats.post.samples
    nsamp <-  multi.chains.retval$nsamp
    spectrum.cdc <- multi.chains.retval$spectrum.cdc

    spectrum.df <- spectrum.df[,order(spectrum.stats[,2],decreasing=T)]
    spectrum.cdc <- spectrum.cdc[,order(spectrum.stats[,2],decreasing=T)]
    spectrum.stats <- spectrum.stats[order(spectrum.stats[,2],decreasing=T),]

    high.confident.spectrum <- spectrum.df[,which(spectrum.stats[,2]>=(confident.prop*nsamp))]
    high.confident.stats <- spectrum.stats[which(spectrum.stats[,2]>=(confident.prop*nsamp)),]
    high.confident.cdc <- spectrum.cdc[,which(spectrum.stats[,2]>=(confident.prop*nsamp))]

    moderate.spectrum <- spectrum.df[,intersect(which(spectrum.stats[,2]>=(noise.prop*nsamp)),which(spectrum.stats[,2]<(confident.prop*nsamp)))]
    moderate.stats <- spectrum.stats[intersect(which(spectrum.stats[,2]>=(noise.prop*nsamp)),which(spectrum.stats[,2]<(confident.prop*nsamp))),]
    moderate.cdc <- spectrum.cdc[,intersect(which(spectrum.stats[,2]>=(noise.prop*nsamp)),which(spectrum.stats[,2]<(confident.prop*nsamp)))]

    noise.spectrum <- spectrum.df[,which(spectrum.stats[,2]<(noise.prop*nsamp))]
    noise.stats <- spectrum.stats[which(spectrum.stats[,2]<(noise.prop*nsamp)),]
    noise.cdc <- spectrum.cdc[,which(spectrum.stats[,2]<(noise.prop*nsamp))]

    noise.spectrum <- cbind(noise.spectrum,multi.chains.retval$each.chain.noise.spectrum)
    noise.cdc <- cbind(noise.cdc,multi.chains.retval$each.chain.noise.cdc)

    extractedSignatures <- high.confident.spectrum

    rownames(extractedSignatures) <- rownames(input.catalog)
    # Set signature names to "hdp.0","hdp.1","hdp.2", ...
    colnames(extractedSignatures) <-
      paste("hdp", c(1:ncol(extractedSignatures)), sep = ".")
    combinedSignatures <- extractedSignatures

    potentialSignatures <- moderate.spectrum
    if(!is.null(potentialSignatures) && ncol(potentialSignatures)>0){
      potentialSignatures <- apply(potentialSignatures,2,function(x)x/sum(x))

      rownames(potentialSignatures) <- rownames(input.catalog)
      # Set signature names to "hdp.0","hdp.1","hdp.2", ...
      colnames(potentialSignatures) <-
        paste("potential hdp", c(1:ncol(potentialSignatures)), sep = ".")
      combinedSignatures <- cbind(extractedSignatures,potentialSignatures)

    }

    combinedSignatures <- apply(combinedSignatures,2,function(x)x/sum(x))
    combined.stats <- rbind(high.confident.stats,moderate.stats)
    combined.cdc  <- cbind(high.confident.cdc,moderate.cdc)

    #sigmatchretval <- apply(combinedSignatures,2,function(x){
    #  hdpx::extract_ccc_cdc_from_hdp(x,
    ##                                 ccc_0 = multi.chains.retval$ccc_0,
    #                                 cdc_0 = multi.chains.retval$cdc_0,
    #                                 cos.merge = 0.95)})

    ## Calculate the exposure probability of each signature (component) for each
    ## tumor sample (posterior sample corresponding to a Dirichlet process node).
    ## This is the probability distribution of signatures (components) for all
    ## tumor samples (DP nodes).

    if (verbose) message("extracting signatures exposures ", Sys.time())
    #exposureProbs <- do.call(cbind,lapply(sigmatchretval,function(x)x[["cdc_mean"]]))
    exposureProbs <- t(apply(combined.cdc,1,function(x){x/sum(x)}))
    # Remove columns corresponding to parent or grandparent nodes
    # (leaving only columns corresponding to samples.
    # Transpose so it conforms to SynSigEval format
    exposureProbs <- t(exposureProbs[-c(1:(nrow(exposureProbs)-ncol(input.catalog))), ])

    # Now rows are signatures, columns are samples
    # Calculate exposure counts from exposure probabilities and total mutation
    # counts

    #if sample doesn't have mutation count, prob is NA
    if(sum(is.na(exposureProbs))>0){
      exposureProbs[which(is.na(exposureProbs))]<-0
    }
   # exposureProbs <- exposureProbs[1:ncol(extractedSignatures),]
    exposureCounts <- exposureProbs %*% diag(rowSums(convSpectra))

    colnames(exposureCounts) <- colnames(input.catalog)

    row.names(exposureCounts) <- colnames(combinedSignatures)

    colnames(exposureProbs) <- colnames(input.catalog)

    row.names(exposureProbs) <- colnames(combinedSignatures)



    return(invisible(list(signature       = combinedSignatures,
                          post.stats      = combined.stats,
                          exposureCounts  = exposureCounts,
                          exposure        = exposureProbs,
                          noise.spectrum  = noise.spectrum,
                          noise.stats     = noise.stats,
                          noise.cdc       = noise.cdc,
                          all.cdc         = spectrum.cdc,
                          extracted.retval = multi.chains.retval)))

  }
