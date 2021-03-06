% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Generateppindex.R
\name{Generateppindex}
\alias{Generateppindex}
\title{Generate index for a HDP structure and num.tumor.types for other functions}
\usage{
Generateppindex(multi.types, input.catalog)
}
\arguments{
\item{multi.types}{A logical scalar or
a character vector.
If \code{FALSE}, The HDP analysis
will regard all input spectra as one tumor type.
HDP structure as one parent node for all tumors

If \code{TRUE}, the HDP analysis
will infer tumor types based on the string before "::" in their names.
e.g. tumor type for "Ovary-AdenoCA::S.500" would be "Ovary-AdenoCA"
HDP structure as a grandparent node for whole data and
one parent node for each tumor type

If \code{multi.types} is a character vector, then it should be of the same length
as the number of columns in \code{input.catalog}, and each value is the
name of the tumor type of the corresponding column in \code{input.catalog}.

If not \code{FALSE}, HDP will give a parent node for each tumor type and a grandparent node
for the whole dataset.

e.g. \code{c("Ovary-AdenoCA", "Kidney-RCC")}.}

\item{input.catalog}{Input spectra catalog as a matrix or
in \code{\link[ICAMS]{ICAMS}} format.}
}
\description{
Generate index for a HDP structure and num.tumor.types for other functions
}
