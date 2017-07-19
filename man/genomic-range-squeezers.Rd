\name{genomic-range-squeezers}

\alias{genomic-range-squeezers}

\alias{granges}
\alias{grglist}

\alias{grglist,Pairs-method}

\title{Squeeze the genomic ranges out of a range-based object}

\description{
  S4 generic functions for squeezing the genomic ranges out of a range-based
  object.

  These are analog to range squeezers \code{\link[IRanges]{ranges}} and
  \code{\link[IRanges]{rglist}} defined in the \pkg{IRanges} package, except
  that \code{granges} returns the ranges in a \link{GRanges} object (instead
  of an \link[IRanges]{IRanges} object for \code{\link[IRanges]{ranges}}),
  and \code{grglist} returns them in a \link{GRangesList} object (instead of
  an \link[IRanges]{IRangesList} object for \code{\link[IRanges]{rglist}}).
}

\usage{
granges(x, use.names=TRUE, use.mcols=FALSE, ...)
grglist(x, use.names=TRUE, use.mcols=FALSE, ...)
}

\arguments{
  \item{x}{
    An object containing genomic ranges e.g. a
    \link[GenomicRanges]{GenomicRanges},
    \link[SummarizedExperiment]{RangedSummarizedExperiment},
    \link[GenomicAlignments]{GAlignments},
    \link[GenomicAlignments]{GAlignmentPairs}, or
    \link[GenomicAlignments]{GAlignmentsList} object,
    or a \link[S4Vectors]{Pairs} object containing genomic ranges.
  }
  \item{use.names, use.mcols, ...}{
    See \code{\link[IRanges]{ranges}} in the \pkg{IRanges} package for
    a description of these arguments.
  }
}

\details{
  See \code{\link[IRanges]{ranges}} in the \pkg{IRanges} package for
  some details.

  For some objects (e.g. \link[GenomicAlignments]{GAlignments} and
  \link[GenomicAlignments]{GAlignmentPairs} objects defined in the
  \pkg{GenomicAlignments} package), \code{as(x, "GRanges")} and
  \code{as(x, "GRangesList")}, are equivalent to
  \code{granges(x, use.names=TRUE, use.mcols=TRUE)} and
  \code{grglist(x, use.names=TRUE, use.mcols=TRUE)}, respectively.
}

\value{
  A \link{GRanges} object for \code{granges}.

  A \link{GRangesList} object for \code{grglist}.

  If \code{x} is a vector-like object (e.g.
  \link[GenomicAlignments]{GAlignments}), the returned object is expected
  to be \emph{parallel} to \code{x}, that is, the i-th element in the output
  corresponds to the i-th element in the input.

  If \code{use.names} is TRUE, then the names on \code{x}
  (if any) are propagated to the returned object.
  If \code{use.mcols} is TRUE, then the metadata columns on \code{x}
  (if any) are propagated to the returned object.
}

\author{H. Pagès}

\seealso{
  \itemize{
    \item \link{GRanges} and \link{GRangesList} objects.

    \item \link[SummarizedExperiment]{RangedSummarizedExperiment} objects
          in the \pkg{SummarizedExperiment} packages.

    \item \link[GenomicAlignments]{GAlignments},
          \link[GenomicAlignments]{GAlignmentPairs},
          and \link[GenomicAlignments]{GAlignmentsList} objects in the
          \pkg{GenomicAlignments} package.
  }
}

\examples{
## See ?GAlignments in the GenomicAlignments package for examples of
## "ranges" and "rglist" methods.
}

\keyword{methods}
