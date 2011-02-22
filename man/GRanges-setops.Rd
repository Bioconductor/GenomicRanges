\name{GRanges-setops}
\alias{GRanges-setops}

\alias{union,GRanges,GRanges-method}
\alias{intersect,GRanges,GRanges-method}
\alias{setdiff,GRanges,GRanges-method}

\alias{punion,GRanges,GRanges-method}
\alias{punion,GRanges,GRangesList-method}
\alias{punion,GRangesList,GRanges-method}

\alias{pintersect,GRanges,GRanges-method}
\alias{pintersect,GRanges,GRangesList-method}
\alias{pintersect,GRangesList,GRanges-method}

\alias{psetdiff,GRanges,GRanges-method}
\alias{psetdiff,GRanges,GRangesList-method}
\alias{psetdiff,GRangesList,GRangesList-method}


\title{GRanges and GRangesList Set and Parallel Set Operations}

\description{
  Set and parallel set operations for GRanges/GRangesList objects.
}

\usage{
## Set operations
\S4method{union}{GRanges,GRanges}(x, y, ignore.strand = FALSE,...)
\S4method{intersect}{GRanges,GRanges}(x, y, ignore.strand = FALSE,...)
\S4method{setdiff}{GRanges,GRanges}(x, y, ignore.strand = FALSE,...)

## Parallel set operations
\S4method{punion}{GRanges,GRanges}(x, y, fill.gap = FALSE, ignore.strand = FALSE, ...)
\S4method{pintersect}{GRanges,GRanges}(x, y, resolve.empty = c("none", "max.start", "start.x"), ignore.strand = FALSE, ...)
\S4method{psetdiff}{GRanges,GRanges}(x, y, ignore.strand = FALSE, ...)
}

\arguments{
  \item{x, y}{
    \link{GRanges} or \link{GRangesList} objects of equal length
    (i.e. \code{length(x) == length(y)}).
    For \code{union}, \code{intersect}, \code{setdiff}, \code{x} and
    \code{y} must both be \link{GRanges} objects.
    For \code{pintersect}, \code{x} and \code{y} cannot both be
    \link{GRangesList} objects.
  }
  \item{fill.gap}{
    Logical indicating whether or not to force a union by using the rule
    \code{start = min(start(x), start(y)), end = max(end(x), end(y))}.
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
  \item{ignore.strand}{
    If set to TRUE, the strand information is ignored in the computation and the 
    result has the strand information of x.
  }
  \item{...}{
    Further arguments to be passed to or from other methods.
  }
}

\details{
  The \code{pintersect} methods involving \link{GRanges} and
  \link{GRangesList} objects use the triplet (sequence name, range,
  strand) to determine the element by element intersection of features,
  where a strand value of \code{"*"} is treated as occurring on both the
  \code{"+"} and \code{"-"} strand.

  The \code{psetdiff} methods involving \link{GRanges} and
  \link{GRangesList} objects use the triplet (sequence name, range,
  strand) to determine the element by element set difference of features,
  where a strand value of \code{"*"} is treated as occurring on both the
  \code{"+"} and \code{"-"} strand.
}

\value{
  For \code{union}, \code{intersect}, and \code{setdiff} a \link{GRanges}.

  For \code{punion} and \code{pintersect} either a \link{GRanges} object
  when both \code{x} and \code{y} are \link{GRanges} objects or a
  \link{GRangesList} object when one of the arguments is a
  \link{GRangesList} object.

  For \code{psetdiff} either a \link{GRanges} object when both \code{x}
  and \code{y} are \link{GRanges} objects or a \link{GRangesList} object
  when \code{y} is a \link{GRangesList} object.
}

\author{P. Aboyoun}

\seealso{
  \link{GRanges},
  \link{GRangesList},
  \link{findOverlaps,GenomicRanges,GenomicRanges-method},
  \link{pintersect,GRanges,GRanges-method}
}

\examples{
  ## GRanges object
  gr <-
    GRanges(seqnames = c("chr2", "chr1", "chr1"),
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


  ## Union, intersection, and set difference of two GRanges objects
  union(gr2, gr3)
  intersect(gr2, gr3)
  setdiff(gr2, gr3)


  ## Parallel intersection of two GRanges objects
  pintersect(gr2, shift(gr2, 3))

  ## Parallel intersection of a GRanges and a GRangesList object
  pintersect(gr, grlist)
  pintersect(grlist, gr)


  ## Parallel set difference of two GRanges objects
  psetdiff(gr2, shift(gr2, 3))

  ## Parallel set difference of a GRanges and a GRangesList object
  psetdiff(gr, grlist)

  ## Parallel set difference of two GRangesList objects
  psetdiff(grlist, shift(grlist, 3))
}

\keyword{methods}
\keyword{utilities}
