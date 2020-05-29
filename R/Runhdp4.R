#' Run hdp extraction and attribution on a spectra catalog file using hdpx
#'
#' @inheritParams RunhdpInternal4
#'
#' @param input.catalog Either a character string, in which case this
#'   is the path to a file containing a spectra catalog
#'   in \code{\link[ICAMS]{ICAMS}} format, or an \code{\link[ICAMS]{ICAMS}} catalog.
#'
#' @param out.dir Directory that will be created for the output;
#'   if \code{overwrite} is \code{FALSE} then
#'   abort if \code{out.dir} already exits.
#'
#' @param remove.noise Deprecated; ignored
#'
#' @param test.only If > 0, only analyze the first \code{test.only} columns
#'  in \code{input.catalog}.
#'
#' @param overwrite If \code{TRUE} overwrite \code{out.dir} if it exists, otherwise
#'  raise an error.
#'
#' @param plot.extracted.sig If \code{TRUE} then plot the extracted signatures.
#'
#' @return The same list as returned by \code{\link{RunhdpInternal4}}.
#'
#' @details Creates several files in \code{out.dir}. These are:
#'  call.and.session.info.txt, hdp.diagnostics.pdf, Runhdp4.retval.Rdata,
#'  extracted.signatures.csv, extracted.signature.pdf (optional),
#'  inferred.exposures.csv.
#'
#' @export

Runhdp4 <-
  function(input.catalog,
           out.dir,
           CPU.cores           = 1,
           seedNumber          = 1,
           K.guess,
           multi.types         = FALSE,
           remove.noise        = FALSE,
           test.only           = 0,
           overwrite           = FALSE,
           verbose             = TRUE,
           num.posterior       = 4,
           post.burnin         = 4000,
           post.n              = 50,
           post.space          = 50,
           post.cpiter         = 3,
           post.verbosity      = 0,
           cos.merge           = 0.9,
           min.sample          = 1,
           checkpoint.aft.post = NULL,
           plot.extracted.sig  = FALSE) {

    if (mode(input.catalog) == "character") {
      if (verbose) message("Reading input catalog file ", input.catalog)
      spectra <- ICAMS::ReadCatalog(input.catalog, strict = FALSE)
    } else {
      spectra <- input.catalog
    }
    if (test.only > 0) spectra <- spectra[ , 1:test.only]

    ## Create output directory
    if (dir.exists(out.dir)) {
      if (!overwrite) stop(out.dir, " already exits")
      if (verbose) message("Using existing out.dir ", out.dir)
    } else {
      dir.create(out.dir, recursive = T)
      if (verbose) message("Created new out.dir", out.dir)
    }

    utils::capture.output(date(), cat("\n\n"),
                          match.call(), cat("\n\n"),
                          utils::sessionInfo(),
                          file = file.path(out.dir, "call.and.session.info.txt"))

    retval <- mSigHdp::RunhdpInternal4(
      input.catalog       = spectra,
      CPU.cores           = CPU.cores,
      seedNumber          = seedNumber,
      K.guess             = K.guess,
      multi.types         = multi.types,
      num.posterior       = num.posterior,
      verbose             = verbose,
      post.burnin         = post.burnin,
      post.n              = post.n,
      post.space          = post.space,
      post.cpiter         = post.cpiter,
      post.verbosity      = post.verbosity,
      cos.merge           = cos.merge,
      min.sample          = min.sample,
      checkpoint.aft.post = checkpoint.aft.post
    ) # 14 Arguments

    save(retval, file = file.path(out.dir, "Runhdp4.retval.Rdata"))

    multi <- retval[["multi.chains"]] # class hdpSampleMulti
    chains <- hdpx::chains(multi)      # list of hdpSampleChain

    # Plot the diagnostics of sampling chains.
    if (verbose) message("Writing HDP diagnostics")
    graphics::par(mfrow=c(2,2), mar=c(4, 4, 2, 1))
    grDevices::pdf(file = file.path(out.dir,"diagnostics.likelihood.pdf"))
    lapply(chains, hdpx::plot_lik, bty = "L")
    grDevices::dev.off()

    grDevices::pdf(file = file.path(out.dir,"diagnostics.numcluster.pdf"))
    # This is the number of raw clusters sampled along each chain
    lapply(chains, hdpx::plot_numcluster, bty = "L")
    grDevices::dev.off()

    grDevices::pdf(file = file.path(out.dir,"diagnostics.data.assigned.pdf"))
    # This is the number of mutations assigned as a function of
    # the number of raw clusters
    lapply(chains, hdpx::plot_data_assigned, bty = "L")
    grDevices::dev.off()

    grDevices::pdf(file = file.path(out.dir,"diagnostics.comp.size.pdf"))
    graphics::par(mfrow=c(1,1), mar=c(5, 4, 4, 2))
    # Components were already extracted, so this call will work
    hdpx::plot_comp_size(multi, bty="L")
    grDevices::dev.off()

    grDevices::pdf(file = file.path(out.dir,"diagnostics.signatures.pdf"))
    graphics::par(mfrow=c(8, 1), mar = c(1, 1, 1, 1))
    # This plots the component (signature) profiles with
    # 95% credibility intervals
    hdpx::plot_comp_distn(multi)
    grDevices::dev.off()

    # TODO, need argument dpindices and col_comp;
    # Need to return the hdp object (perhaps) from RunhdpInternal
    # to get the required values.
    if (FALSE) { # Not finished
    num.dpindices <- length(chains[[1]]@hdp@ppindex)
    hdpx::plot_dp_comp_exposure(
      multi, dpindices = 3:num.dpindices,
      col_comp = myCol[1:ncol(retval$signature)],
      dpnames = colnames(retval$exposure))
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

    # Probably not needed; easily computed by caller:
    # SynSigGen::WriteExposure(retval$exposure.p,
    #              paste0(out.dir,"/exposure.probs.csv"))
    SynSigGen::WriteExposure(retval$exposure,
                  file.path(out.dir,"inferred.exposures.csv"))

    invisible(retval)
  }

