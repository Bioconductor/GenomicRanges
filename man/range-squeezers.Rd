\name{range-squeezers}

\alias{range-squeezers}

\alias{granges}
\alias{grglist}
\alias{rglist}


\title{Squeeze the ranges out of a range-based object}

\description{
  S4 generic functions for squeezing the ranges out of a range-based object.

  \code{granges} returns them as a \link{GRanges} object, \code{grglist}
  as a \link{GRangesList} object, and \code{rglist} as a
  \link[IRanges]{RangesList} object.
}

\usage{
granges(x, ...)
grglist(x, ...)
rglist(x, ...)
}

\arguments{
  \item{x}{
    A range-based object e.g. a \link{SummarizedExperiment},
    \link[GenomicAlignments]{GAlignments},
    \link[GenomicAlignments]{GAlignmentPairs},
    or \link[GenomicAlignments]{GAlignmentsList} object.
  }
  \item{...}{
    Additional arguments, for use in specific methods.
  }
}

\details{
  The \pkg{GenomicRanges} and \pkg{GenomicAlignments} packages define and
  document methods for various types of range-based objects (e.g. for
  \link{SummarizedExperiment}, \link[GenomicAlignments]{GAlignments},
  \link[GenomicAlignments]{GAlignmentPairs}, and
  \link[GenomicAlignments]{GAlignmentsList} objects).
  Other Bioconductor packages might as well.

  Note that these functions can be seen as a specific kind of \emph{object
  getters} as well as functions performing coercion.
  For some objects (e.g. \link[GenomicAlignments]{GAlignments}),
  \code{granges(x)}, \code{grglist(x)}, and \code{rglist(x)} are actually
  equivalent to \code{as(x, "GRanges")}, \code{as(x, "GRangesList")},
  and \code{as(x, "RangesList")}, respectively.
}

\value{
  A \link{GRanges} object for \code{granges}.

  A \link{GRangesList} object for \code{grglist}.

  A \link[IRanges]{RangesList} object for \code{rglist}.

  If \code{x} is a vector-like object (e.g.
  \link[GenomicAlignments]{GAlignments}), the returned object is expected
  to be \emph{parallel} to \code{x}, that is, the i-th element in the output
  corresponds to the i-th element in the input.
}

\author{H. Pages}

\seealso{
  \itemize{
    \item \link{GRanges} and \link{GRangesList} objects.

    \item \link[IRanges]{RangesList} objects in the \pkg{IRanges} package.

    \item \link{SummarizedExperiment} objects.

    \item \link[GenomicAlignments]{GAlignments},
          \link[GenomicAlignments]{GAlignmentPairs},
          and \link[GenomicAlignments]{GAlignmentsList} objects in the
          \pkg{GenomicAlignments} package.
  }
}

\examples{
## See ?GAlignments in the GenomicAlignments package for some
## examples.
}

\keyword{methods}
