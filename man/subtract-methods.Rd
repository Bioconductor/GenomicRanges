\name{subtract-methods}
\alias{subtract-methods}

\alias{subtract}
\alias{subtract,GenomicRanges,GenomicRanges-method}

\title{Subtract a set of genomic ranges from a GRanges object}

\description{
  Similar to bedtools subtract.
}

\usage{
subtract(x, y, minoverlap=1L, ...)

\S4method{subtract}{GenomicRanges,GenomicRanges}(x, y, minoverlap=1L, ignore.strand=FALSE)
}

\arguments{
  \item{x, y}{
    Two \link{GRanges} objects, typically, but any \link{GenomicRanges}
    derivative should be supported.
    Note that \code{y} gets immediately replaced with:
    \preformatted{  reduce(y, ignore.strand=ignore.strand)
    } internally.
  }
  \item{minoverlap}{
    Minimum overlap (in number of genomic positions) between a range in
    \code{x} and a range in \code{reduce(y, ignore.strand=ignore.strand)}
    for the 2 ranges to be considered overlapping, and for their overlapping
    portion to be removed from the range in \code{x}.
  }
  \item{ignore.strand}{
    If set to TRUE, the strand information is ignored in the computation
    and the strand of \code{x} is propagated to the result.
  }
  \item{...}{
    Further arguments to be passed to specific methods.
  }
}

\details{
  \code{subtract()} first replaces its second argument \code{y} with:
  \preformatted{  reduce(y, ignore.strand=ignore.strand)
  }

  Then it searches for genomic ranges in \code{y} that overlap genomic
  ranges in \code{x} by at least the number of base pairs specified via
  the \code{minoverlap} argument. If an overlapping range is found in
  \code{y}, the overlapping portion is removed from any range in \code{x}
  involved in the overlap.

  Note that by default \code{subtract(x, y)} is equivalent to:
  \preformatted{  psetdiff(x, rep(GRangesList(y), length(x)))
  } but will typically be hundred times more efficient.
}

\value{
  A \link{GRangesList} object \emph{parallel} to \code{x}, that is, with
  one list element per range in \code{x}.
  The names and metadata columns on \code{x} are propagated to the result.
}

\author{H. Pagès}

\seealso{
  \itemize{
    \item bedtools subtract at
          \url{https://bedtools.readthedocs.io/en/latest/content/tools/subtract.html}

    \item \link[GenomicRanges]{setops-methods} for set operations on
          \link{GRanges} objects.

    \item \link[GenomicRanges]{findOverlaps-methods} for finding/counting
          overlapping genomic ranges.

    \item \link[GenomicRanges]{intra-range-methods} and
          \link[GenomicRanges]{inter-range-methods} for \emph{intra range}
          and \emph{inter range} transformations of a \link{GRanges} object.

    \item \link{GRanges} and \link{GRangesList} objects.
  }
}

\examples{
x <- GRanges(c(A="chr1:1-50", B="chr1:40-110", C="chrX:1-500"))
y <- GRanges(c("chr1:21-25", "chr1:38-150"))
z <- subtract(x, y)
z

unlist(z)
}

\keyword{methods}
\keyword{utilities}
