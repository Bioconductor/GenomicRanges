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
  See the "Overview of summarizeOverlaps" vignette in the GenomicRanges
  package for details about the BAM files located at
  \url{http://www.embl.de/~reyes/Graveley/bam}.
  Files treated1.bam, untreated1.bam and untreated2.bam contain reads
  aligned to reference genome \emph{BDGP Release 5} (aka the \emph{dm3}
  genome on the UCSC Genome Browser).

  The \code{galn4} data set was obtained by extracting the records from
  untreated1.bam that are located on chr4.
}

\seealso{
  \link{GappedAlignments-class},
  \code{\link{isCompatibleWithSplicing}}
}

\examples{
data(galn4)
galn4  # 204355 alignments
}

\keyword{data}
\keyword{datasets}
