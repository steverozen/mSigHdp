#' Extract components and exposures from multiple posterior sample chains
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
#' @param min.sample A "component" (i.e. signature) must have at least
#'      this many samples; passed to \code{\link[hdpx]{hdp_merge_and_extract_components}}.
#'
#' @param categ.CI A number the range \eqn{[0,1]}. The level of the confidence
#'   interval used in step 4 of \code{\link{hdp_merge_and_extract_components}}.
#'   This governs when "averaged raw cluster" get assigned to component 0,
#'   i.e. if the the confidence interval overlaps 0. Lower values
#'   make it less likely that an averaged raw cluster will be assigned to
#'   component 0. The CI in question is for the number of mutations in
#'   a given mutation class (e.g. ACA > AAA, internally called a
#'   "category"). If, for every mutation class, this CI overlaps 0,
#'   then the averaged raw cluster goes to component 0.
#'
#' @param exposure.CI A number in the range \eqn{[0,1]}. The level of
#'   the confidence interval used in step 5 of hdp_merge_and_extract_components.
#'   The CI in question here for the total number of
#'   mutations assigned to an averaged raw cluster.
#'
#' @param diagnostic.folder If provided, diagnostic plots for hdp.0 components are provided
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
#' \item{sum_raw_clusters_after_cos_merge}{A matrix containing aggregated spectra of raw clusters after cosine
#'       similarity merge step in \code{\link[hdpx]{hdp_merge_and_extract_components}}.}
#'
#' \item{sum_raw_clusters_after_nonzero_categ}{A matrix containing aggregated spectra of raw clusters after non-zero category selecting
#'       step in \code{\link[hdpx]{hdp_merge_and_extract_components}}.}
#'
#' \item{clust_hdp0_ccc4}{A matrix containing aggregated spectra of raw clusters moving to hdp.0 after non-zero category selection step in              \code{\link[hdpx]{hdp_merge_and_extract_components}}.}
#'
#' \item{clust_hdp0_ccc5}{A matrix containing aggregated spectra of raw clusters moving to hdp.0 after non-zero observation selection step in           \code{\link[hdpx]{hdp_merge_and_extract_components}}.}
#'
#' }
#'
#'
#' @export
#'
CombinePosteriorChains <-
  function(clean.chlist,
           input.catalog,
           multi.types,
           verbose             = TRUE,
           cos.merge           = 0.9,
           categ.CI            = 0.95,
           exposure.CI         = 0.95,
           min.sample          = 1,
           diagnostic.folder   = NULL
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
        hdpx::hdp_merge_and_extract_components(multi.chains,
                                               exposure.CI    = exposure.CI,
                                               categ.CI       = categ.CI,
                                               cos.merge      = cos.merge,
                                               min.sample     = min.sample,
                                               diagnostic.folder = diagnostic.folder)
    )

    if (verbose) {
      message("hdp_extract_components time: ")
      for (xn in names(extract.time)) {
        message(" ", xn, " ", extract.time[[xn]])
      }
    }
    if (verbose) message("calling hdpx::comp_categ_distn ", Sys.time())
    extractedSignatures <- t(hdpx::comp_categ_distn(multi.chains)$mean)



    rownames(extractedSignatures) <- rownames(input.catalog)
    # Set signature names to "hdp.0","hdp.1","hdp.2", ...
    colnames(extractedSignatures) <-
      paste("hdp", colnames(extractedSignatures), sep = ".")

    ## Calculate the exposure probability of each signature (component) for each
    ## tumor sample (posterior sample corresponding to a Dirichlet process node).
    ## This is the probability distribution of signatures (components) for all
    ## tumor samples (DP nodes).

    if (verbose) message("Calling hdpx::comp_dp_distn ", Sys.time())
    exposureProbs <- hdpx::comp_dp_distn(multi.chains)$mean

    # Remove columns corresponding to parent or grandparent nodes
    # (leaving only columns corresponding to samples.
    # Transpose so it conforms to SynSigEval format
    exposureProbs <- t(exposureProbs[-c(1:(nrow(exposureProbs)-ncol(input.catalog))), ])

    # Now rows are signatures, columns are samples
    # Calculate exposure counts from exposure probabilities and total mutation
    # counts
    exposureCounts <- exposureProbs %*% diag(rowSums(convSpectra))

    colnames(exposureCounts) <- colnames(input.catalog)

    if(!(any(grepl("N",colnames(extractedSignatures)))||any(grepl("P",colnames(extractedSignatures))))){
      rownames(exposureCounts) <- colnames(extractedSignatures)
    }

    invisible(list(signature       = extractedSignatures,
                   exposure        = exposureCounts,
                   multi.chains    = multi.chains,

                   sum_raw_clusters_after_cos_merge  =
                     hdpx::comp_categ_distn(multi.chains)$raw_clusts_after_cos_merge,
                   sum_raw_clusters_after_nonzero_categ  =
                     hdpx::comp_categ_distn(multi.chains)$raw_clusts_after_cos_merge,
                   clust_hdp0_ccc4  = hdpx::comp_categ_distn(multi.chains)$clust_hdp0_ccc4,
                   clust_hdp0_ccc5  = hdpx::comp_categ_distn(multi.chains)$clust_hdp0_ccc5))


  }
