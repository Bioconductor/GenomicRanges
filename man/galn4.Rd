\name{galn4}
\docType{data}

\alias{galn4}

\title{Small subset of aligned reads from the Pasilla experiment}

\description{
  \link{GappedAlignments} object containing a small subset of aligned reads
  from the Pasilla experiment.
}

\usage{
data(galn4)
}

\details{
  See the pasilla data package for details about the Pasilla experiment
  and data (RNA-seq, Fly).
  The \code{galn4} data set was obtained by extracting records from
  the untreated1.bam file that are on chr4.
}

\seealso{
  \link{GappedAlignments-class},
  \code{\link{readIsCompatibleWithSplicing}}
}

\examples{
data(galn4)
galn4  # 204355 alignments
}

\keyword{data}
\keyword{datasets}
