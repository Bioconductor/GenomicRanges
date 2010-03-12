\name{GRanges-setops}

\alias{pintersect,GRanges,GRanges-method}
\alias{pintersect,GRanges,GRangesList-method}
\alias{pintersect,GRangesList,GRanges-method}


\title{GRanges and GRangesList Set Operations}

\description{
  Finds interval overlaps between a GRanges/GRangesList object
  and a GRanges/GRangesList object.
}

\usage{
  \S4method{pintersect}{GRanges,GRanges}(x, y, resolve.empty=c("none", "max.start", "start.x"), ...)
}

\arguments{
  \item{x, y}{
    \link{GRanges} or \link{GRangesList} objects of equal length
    (i.e. \code{length(x) == length(y)}).
    For \code{pintersect}, \code{x} and \code{y} cannot both be
    \link{GRangesList} objects.
  }
  \item{resolve.empty}{
    One of \code{"none"}, \code{"max.start"}, or \code{"start.x"} denoting
    how to handle ambiguous empty ranges formed by intersections.
    \code{"none"} - throw an error if an ambiguous empty range is formed,
    \code{"max.start"} - associate the maximum start value with any
    ambiguous empty range, and \code{"start.x"} - associate the start value
    of \code{x} with any ambiguous empty range. (See
    \code{\link[IRanges]{pintersect}} for the definition of an ambiguous
    range.)
  }
  \item{...}{
    Further arguments to be passed to or from other methods.
  }
}

\details{
  The \code{pintersect} methods involving \link{GRanges} and
  \link{GRangesList} objects use the triplet (sequence name, range,
  strand) to determine the element by element intersection of features
  (see paragraph below for the definition of feature), where a strand
  value of \code{"*"} is treated as occurring on both the \code{"+"}
  and \code{"-"} strand.
}

\value{
  For \code{pintersect} either a \link{GRanges} object when both \code{x}
  and \code{y} are \link{GRanges} objects or a \link{GRangesList} object
  when one of the arguments is a \link{GRangesList} object.
}

\author{P. Aboyoun}

\seealso{
  \code{\link[IRanges]{pintersect}},
  \link{GRanges},
  \link{GRangesList}
}

\examples{
  ## GRanges object
  gr <-
    GRanges(seqnames =
            Rle(c("chr2", "chr1", "chr1")),
            ranges = IRanges(1:3, width = 12),
            strand = Rle(strand(c("-", "*", "-"))))
  gr

  ## GRangesList object
  gr1 <-
    GRanges(seqnames = "chr2", ranges = IRanges(3, 6),
            strand = "+", score = 5L, GC = 0.45)
  gr2 <-
    GRanges(seqnames = c("chr1", "chr1"),
            ranges = IRanges(c(7,13), width = 3),
            strand = c("+", "-"), score = 3:4, GC = c(0.3, 0.5))
  gr3 <-
    GRanges(seqnames = c("chr1", "chr2"),
            ranges = IRanges(c(1, 4), c(3, 9)),
            strand = c("-", "-"), score = c(6L, 2L), GC = c(0.4, 0.1))
  grlist <- GRangesList("gr1" = gr1, "gr2" = gr2, "gr3" = gr3)

  ## Parallel intersection of two GRanges objects
  pintersect(gr2, shift(gr2, 3))

  ## Parallel intersection of a GRanges and a GRangesList object
  pintersect(gr, grlist)
  pintersect(grlist, gr)
}

\keyword{methods}
\keyword{utilities}
