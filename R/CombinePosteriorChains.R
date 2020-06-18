#' Extract components and exposures from multiple posterior sample chains
#'
#' @param clean.chlist It collects the output of multiple independent hdp_posterior calls.
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
#'      this many samples; passed to \code{\link[hdpx]{hdp_extract_components}}.
#'
#' @return Invisibly, a list with the following elements:\describe{
#' \item{signature}{The extracted signature profiles as a matrix;
#'             rows are mutation types, columns are
#'             samples (e.g. tumors).}
#' \item{exposure}{The inferred exposures as a matrix of mutation counts;
#'            rows are signatures, columns are samples (e.g. tumors).}
#'
#' \item{multi.chains}{A \code{\link[hdpx]{hdpSampleMulti-class}} object.
#'     This object has the method \code{\link[hdpx]{chains}} which returns
#'     a list of \code{\link[hdpx]{hdpSampleChain-class}} objects. Each of these
#'     sample chains objects has a method \code{\link[hdpx]{final_hdpState}}
#'     (actually the methods seems to be just \code{hdp})
#'     that returns the \code{hdpState} from which it was generated.}
#'}
#' @export
#'
CombinePosteriorChains <-
  function(clean.chlist,
           input.catalog,
           multi.types,
           one.parent.hack,
           verbose             = TRUE,
           cos.merge           = 0.9,
           min.sample          = 1
  ) { # 6 arguments
    if (mode(input.catalog) == "character") {
      if (verbose) message("Reading input catalog file ", input.catalog)
      input.catalog <- ICAMS::ReadCatalog(input.catalog, strict = FALSE)
    } else {
      input.catalog <- input.catalog
    }
    # hdp gets confused if the class of its input is not matrix.
    convSpectra <- t(input.catalog)
    # class(convSpectra) <- "matrix"
    # convSpectra <- t(convSpectra)
    number.channels <- nrow(input.catalog)
    number.samples  <- ncol(input.catalog)

    #this function to generate num.tumor.type
    ppindex <- Generateppindex(multi.types = multi.types,
                               one.parent.hack = one.parent.hack,
                               input.catalog = input.catalog) ##clean up code

    multi.chains <- hdpx::hdp_multi_chain(clean.chlist)
    if (verbose) message("calling hdp_extract_components ", Sys.time())
    # Group raw "clusters" into "components" (i.e. signatures).
    extract.time <- system.time(
      multi.chains <-
        hdpx::hdp_extract_components(multi.chains,
                                     cos.merge  = cos.merge,
                                     min.sample = min.sample)
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
    ## tumor samples (DP nodes); exposureProbs is the normalized
    ## signature exposure all tumor samples # TODO what is this?

    if (verbose) message("Calling hdpx::comp_dp_distn ", Sys.time())
    exposureProbs <- hdpx::comp_dp_distn(multi.chains)$mean
    # Remove columns corresponding to parent or grandparent nodes
    # (leaving only columns corresponding to samples.
    # Transpose so it conforms to SynSigEval format
    exposureProbs <- t(exposureProbs[-(1:(ppindex$num.tumor.types + 1)), ])
    # Now rows are signatures, columns are samples
    # Calculate exposure counts from exposure probabilities and total mutation
    # counts
    exposureCounts <- exposureProbs %*% diag(rowSums(convSpectra))
    colnames(exposureCounts) <- colnames(input.catalog)
    rownames(exposureCounts) <- colnames(extractedSignatures)

    invisible(list(signature       = extractedSignatures,
                   exposure        = exposureCounts,
                   multi.chains    = multi.chains))
  }
