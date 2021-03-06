% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AnalyzeAndPlotretval.R
\name{AnalyzeAndPlotretval}
\alias{AnalyzeAndPlotretval}
\title{Evaluate and plot retval from \code{CombinePosteriorChains} or \code{CombineChainsAndExtractSigs}
This function now works for both NR's pipeline and Mo's pipeline}
\usage{
AnalyzeAndPlotretval(
  retval,
  input.catalog,
  out.dir = NULL,
  ground.truth.sig = NULL,
  ground.truth.exp = NULL,
  verbose = TRUE,
  overwrite = TRUE,
  diagnostic.plot = TRUE
)
}
\arguments{
\item{retval}{the output from function  \code{CombinePosteriorChains} or \code{CombineChainsAndExtractSigs}}

\item{input.catalog}{input catalog matrix or
path to file with input catalog}

\item{out.dir}{Directory that will be created for the output;
if \code{overwrite} is \code{FALSE} then
abort if \code{out.dir} already exits.}

\item{ground.truth.sig}{Optional. Either a string with the
path to file with ground truth signatures or and
\code{\link[ICAMS]{ICAMS}} catalog with the ground truth signatures.
These are the signatures used to construct the ground truth spectra.}

\item{ground.truth.exp}{Optional. Ground truth exposure matrix or
path to file with ground truth exposures.
If \code{NULL} skip checks that need this information.}

\item{verbose}{If \code{TRUE} then \code{message} progress information.}

\item{overwrite}{If \code{TRUE} overwrite \code{out.dir} if it exists, otherwise raise an error.}

\item{diagnostic.plot}{If \code{TRUE} plot diagnostic plot.
This is optional because there are cases having error}
}
\description{
Evaluate and plot retval from \code{CombinePosteriorChains} or \code{CombineChainsAndExtractSigs}
This function now works for both NR's pipeline and Mo's pipeline
}
