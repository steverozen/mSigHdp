#' Extract components and exposures from multiple posterior sample chains
#' A test function combine new extraction code
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
#' @param prop.samp proportion of samples for a signature to be selected
#'
#' @param cos.merge cosine similarity cutoff
#'
#' @param min.sample passed to another function
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
#'}
#'
#'
#'
#'
#' @export
#'
NewCombinePosteriorChains <-
  function(clean.chlist,
           input.catalog,
           multi.types,
           verbose             = TRUE,
           cos.merge           = 0.9,
           prop.samp           = 0.5,
           min.sample          = 1
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
      multi.chains <-
        hdpx::extract_sigs_from_clusters(multi.chains,
                                               prop.samp       = prop.samp,
                                               cos.merge      = cos.merge,
                                               min.sample     = min.sample)
    )

    if (verbose) {
      message("hdp_extract_components time: ")
      for (xn in names(extract.time)) {
        message(" ", xn, " ", extract.time[[xn]])
      }
    }
    if (verbose) message("calling hdpx::comp_categ_distn ", Sys.time())
    extractedSignatures <- hdpx::comp_categ_distn(multi.chains)$comp_categ_distn
    noise.individual.clusters <- hdpx::comp_categ_distn(multi.chains)$noise.individual.clusters
    nonselected.matched.samples <- hdpx::comp_categ_distn(multi.chains)$nonselected.matched.samples


    rownames(extractedSignatures) <- rownames(input.catalog)
    # Set signature names to "hdp.0","hdp.1","hdp.2", ...
    colnames(extractedSignatures) <-
      paste("hdp", c(0:(ncol(extractedSignatures)-1)), sep = ".")

    ## Calculate the exposure probability of each signature (component) for each
    ## tumor sample (posterior sample corresponding to a Dirichlet process node).
    ## This is the probability distribution of signatures (components) for all
    ## tumor samples (DP nodes).

    if (verbose) message("Calling hdpx::comp_dp_distn ", Sys.time())
    exposureProbs <- hdpx::comp_dp_distn(multi.chains)$mean

    # Remove columns corresponding to parent or grandparent nodes
    # (leaving only columns corresponding to samples.
    # Transpose so it conforms to SynSigEval format
    #exposureProbs <- t(exposureProbs[-c(1:(nrow(exposureProbs)-ncol(input.catalog))), ])

    # Now rows are signatures, columns are samples
    # Calculate exposure counts from exposure probabilities and total mutation
    # counts
    exposureCounts <- diag(rowSums(convSpectra)) %*% exposureProbs

    rownames(exposureCounts) <- colnames(input.catalog)
    colnames(exposureCounts) <- colnames(extractedSignatures)

    #if(!(any(grepl("N",colnames(extractedSignatures)))||any(grepl("P",colnames(extractedSignatures))))){
     # rownames(t(exposureCounts)) <- colnames(extractedSignatures)
    #}

    invisible(list(signature       = extractedSignatures,
                   exposure        = t(exposureCounts),
                   multi.chains    = multi.chains,
                   noise.individual.clusters = noise.individual.clusters,
                   noise.clusters.posterior.samples = nonselected.matched.samples
    ))


  }
