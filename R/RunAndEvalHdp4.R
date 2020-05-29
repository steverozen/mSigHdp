#' Run and evaluate hdp
#'
#' @inheritParams Runhdp4
#'
#' @param ground.truth.exp Ground truth exposure matrix or
#'   path to file with ground truth exposures.
#'   If \code{NULL} skip checks that need this information.
#'
#' @param ground.truth.sig Use this one. Either a string with the
#'   path to file with ground truth signatures or and
#'   \code{\link[ICAMS]{ICAMS}} catalog with the ground truth signatures.
#'   These are the signatures used to construct the ground truth spectra.
#'
#' @export
#'
#' @return See the return value for \code{\link{RunhdpInternal4}}.

RunAndEvalHdp4 <- function(
  input.catalog,
  ground.truth.exp         = NULL,
  ground.truth.sig,
  out.dir,
  CPU.cores                  = 1,
  seedNumber                 = 1,
  K.guess,
  multi.types                = FALSE,
  remove.noise               = FALSE,
  test.only                  = 0,
  overwrite                  = FALSE,
  verbose                    = TRUE,
  num.posterior              = 4,
  post.burnin                = 4000,
  post.n                     = 50,
  post.space                 = 50,
  post.cpiter                = 3,
  post.verbosity             = 0,
  cos.merge                  = 0.9,
  min.sample          = 1,
  checkpoint.aft.post = NULL) {

  if (mode(ground.truth.sig) == "character") {
    if (verbose) {
      message("Reading ground truth signatures from ",
              ground.truth.sig)
    }
    ground.truth.sig <- ICAMS::ReadCatalog(ground.truth.sig, strict = FALSE)
  }

  stopifnot(is.matrix(ground.truth.sig))

  # Do this early to catch any possible error before we do a lot
  # of computation
  if (!is.null(ground.truth.exp) &&
      class(ground.truth.exp) == "character") {
    ground.truth.exp <- SynSigGen::ReadExposure(ground.truth.exp)
  }

  retval <- Runhdp4(
    input.catalog         = input.catalog,
    out.dir               = out.dir,
    CPU.cores             = CPU.cores,
    seedNumber            = seedNumber,
    K.guess               = K.guess,
    multi.types           = multi.types,
    test.only             = test.only,
    overwrite             = overwrite,
    verbose               = verbose,
    num.posterior         = num.posterior,
    post.burnin           = post.burnin,
    post.n                = post.n,
    post.space            = post.space,
    post.cpiter           = post.cpiter,
    post.verbosity        = post.verbosity,
    cos.merge             = cos.merge,
    min.sample            = min.sample,
    checkpoint.aft.post   = checkpoint.aft.post)

  sigAnalysis0 <- SynSigEval::MatchSigsAndRelabel(
    ex.sigs  = retval$signature,
    gt.sigs  = ground.truth.sig,
    exposure = ground.truth.exp)

  ICAMS::PlotCatalogToPdf(ICAMS::as.catalog(sigAnalysis0$gt.sigs,
                                            catalog.type = "counts.signature"), # Need to fix this
                          file.path(out.dir, "ground.truth.sigs.w0.pdf"))

  ICAMS::PlotCatalogToPdf(ICAMS::as.catalog(sigAnalysis0$ex.sigs,
                                            catalog.type = "counts.signature"),
                          file.path(out.dir, "extracted.sigs.w0.pdf"))

  # Writes bi-directional matching and cos.sim calculation
  utils::write.csv(sigAnalysis0$match1, file = file.path(out.dir, "match1.w0.csv"))
  utils::write.csv(sigAnalysis0$match2, file = file.path(out.dir, "match2.w0.csv"))

  utils::capture.output(
    cat("Call\n"),
    match.call(),
    cat("\nAverage cosine similarity\n"),
    sigAnalysis0$averCosSim,
    cat("\nAverage cosine similarity to each ground-truth signature\n"),
    sigAnalysis0$cosSim,
    cat("\nNumber of ground-truth signatures\n"),
    ncol(sigAnalysis0$gt.sigs),
    cat("\nNumber of extracted signatures\n"),
    ncol(sigAnalysis0$ex.sigs),
    cat("\nsigAnalysis0$extracted.with.no.best.match\n"),
    sigAnalysis0$extracted.with.no.best.match,
    cat("\nsigAnalysis0$ground.truth.with.no.best.match\n"),
    sigAnalysis0$ground.truth.with.no.best.match,
    file = file.path(out.dir,"other.results.txt"))

  # Put analysis and plotting of exposures here

  invisible(retval)

}
