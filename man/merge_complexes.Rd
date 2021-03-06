% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bowtie_functions.R
\name{merge_complexes}
\alias{merge_complexes}
\title{Merge overlapping submatrices (protein complexes) produced by find.complexes}
\usage{
merge_complexes(
  prot_complex_ranges,
  merge_if_overlap = 1,
  min_complex_size = 8
)
}
\arguments{
\item{prot_complex_ranges}{A list of start and end ranges that indicate areas of protein complexes.}

\item{merge_if_overlap}{Complexes are merged if they have this number of overlapping proteins (default: 1). A value of 1 will produce no overlapping protein complexes.}

\item{min_complex_size}{The minimum size of a protein complex (default: 8).}
}
\value{
Returns a list of vectors that specifiy the start and end range of each submatrix
}
\description{
Merge submatrices given by find.complexes if they overlap and apply a minimum size threshold.
}
\examples{

}
\keyword{submatrices}
