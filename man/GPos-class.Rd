\name{GPos-class}
\docType{class}

\alias{class:GPos}
\alias{GPos-class}
\alias{GPos}

\alias{length,GPos-method}
\alias{names<-,GPos-method}
\alias{seqnames,GPos-method}
\alias{pos}
\alias{pos,GPos-method}
\alias{start,GPos-method}
\alias{end,GPos-method}
\alias{width,GPos-method}
\alias{ranges,GPos-method}
\alias{strand,GPos-method}
\alias{seqinfo,GPos-method}

\title{GPos objects}

\description{
  The GPos class is a container for storing a set of genomic \emph{positions},
  that is, genomic ranges of length 1. It's more memory-efficient than
  the \link{GRanges} container when the object contains long runs of adjacent
  positions.
}

\usage{
GPos(pos_runs)
}

\arguments{
  \item{pos_runs}{
    A \link{GRanges} object (or any \link{GenomicRanges} object) where each
    range is interpreted as a run of adjacent genomic positions.
  }
}

\details{
}

\author{Hervé Pagès}

\seealso{
  \itemize{
    \item \link{GRanges} objects.
  }
}

\examples{
pos_runs <- GRanges(c("chr1", "chr1", "chr2"),
                    IRanges(c(44, 5, 2), width=c(10, 6, 4)))
GPos(pos_runs)
}
\keyword{methods}
\keyword{classes}
