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
#'      passed to \code{\link[hdpx]{hdp_extract_components}}.
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
#' #'
#' }

#'
#' @export
#'
CombineChainsAndExtractSigs <-
  function(clean.chlist,
           input.catalog,
           multi.types,
           verbose             = TRUE,
           cos.merge           = 0.9
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
    if (verbose) message("calling hdp_extract_components ", Sys.time())
    # Group raw "clusters" into "components" (i.e. signatures).

    extract.time <- system.time(
      multi.chains.retval <-
        hdpx::extract_sigs_from_clusters(multi.chains,
                                         cos.merge      = cos.merge
        )
    )


    if (verbose) {
      message("hdp_extract_components time: ")
      for (xn in names(extract.time)) {
        message(" ", xn, " ", extract.time[[xn]])
      }
    }
    if (verbose) message("calling hdpx::comp_categ_distn ", Sys.time())
    extractedSignatures <- multi.chains.retval$high.confident.spectrum
    extractedSignatures <- apply(extractedSignatures,2,function(x)x/sum(x))

    rownames(extractedSignatures) <- rownames(input.catalog)
    # Set signature names to "hdp.0","hdp.1","hdp.2", ...
    colnames(extractedSignatures) <-
      paste("hdp", c(1:ncol(extractedSignatures)), sep = ".")

    sigmatchretval <- apply(extractedSignatures,2,function(x){
      hdpx::extract_ccc_cdc_from_hdp(x,ccc_0 = multi.chains.retval$ccc_0,
                                     cdc_0 = multi.chains.retval$cdc_0,cos.merge = 0.9)})

    ## Calculate the exposure probability of each signature (component) for each
    ## tumor sample (posterior sample corresponding to a Dirichlet process node).
    ## This is the probability distribution of signatures (components) for all
    ## tumor samples (DP nodes).

    if (verbose) message("Calling hdpx::comp_dp_distn ", Sys.time())
    exposureProbs <- do.call(cbind,lapply(sigmatchretval,function(x)x[["cdc_mean"]]))
    exposureProbs <- t(apply(exposureProbs,1,function(x){x/sum(x)}))
    # Remove columns corresponding to parent or grandparent nodes
    # (leaving only columns corresponding to samples.
    # Transpose so it conforms to SynSigEval format
    exposureProbs <- t(exposureProbs[-c(1:(nrow(exposureProbs)-ncol(input.catalog))), ])

    # Now rows are signatures, columns are samples
    # Calculate exposure counts from exposure probabilities and total mutation
    # counts
    exposureCounts <- exposureProbs %*% diag(rowSums(convSpectra))

    colnames(exposureCounts) <- colnames(input.catalog)

    row.names(exposureCounts) <- colnames(extractedSignatures)

    return(invisible(list(signature       = extractedSignatures,
                          exposure        = exposureCounts,
                          multi.chains    = multi.chains.retval$multi.chains,
                          extracted.retval = multi.chains.retval,
                          diagnostic.retval = sigmatchretval)))

  }
