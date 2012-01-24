\name{setops-methods}
\alias{setops-methods}

\alias{union,GRanges,GRanges-method}
\alias{intersect,GRanges,GRanges-method}
\alias{setdiff,GRanges,GRanges-method}

\alias{punion,GRanges,GRanges-method}
\alias{punion,GRanges,GRangesList-method}
\alias{punion,GRangesList,GRanges-method}

\alias{pintersect,GRanges,GRanges-method}
\alias{pintersect,GRanges,GRangesList-method}
\alias{pintersect,GRangesList,GRanges-method}
\alias{pintersect,GRangesList,GRangesList-method}
\alias{pintersect,GappedAlignments,GRanges-method}
\alias{pintersect,GRanges,GappedAlignments-method}

\alias{psetdiff,GRanges,GRanges-method}
\alias{psetdiff,GRanges,GRangesList-method}
\alias{psetdiff,GRangesList,GRangesList-method}

\alias{pgap,GRanges,GRanges-method}

\title{Set operations on GRanges/GRangesList/GappedAlignments objects}

\description{
  Performs set operations on GRanges/GRangesList/GappedAlignments objects.
}

\usage{
## Set operations
\S4method{union}{GRanges,GRanges}(x, y, ignore.strand=FALSE, ...)
\S4method{intersect}{GRanges,GRanges}(x, y, ignore.strand=FALSE, ...)
\S4method{setdiff}{GRanges,GRanges}(x, y, ignore.strand=FALSE, ...)

## Parallel set operations
\S4method{punion}{GRanges,GRanges}(x, y, fill.gap=FALSE, ignore.strand=FALSE, ...)
\S4method{pintersect}{GRanges,GRanges}(x, y, resolve.empty=c("none", "max.start", "start.x"), ignore.strand=FALSE, ...)
\S4method{pintersect}{GappedAlignments,GRanges}(x, y, ...)
\S4method{psetdiff}{GRanges,GRanges}(x, y, ignore.strand=FALSE, ...)
}

\arguments{
  \item{x, y}{
    For \code{union}, \code{intersect}, \code{setdiff}, \code{pgap}:
    \code{x} and \code{y} must both be \link{GRanges} objects.

    For \code{punion}: one of \code{x} or \code{y} must be a \link{GRanges}
    object, the other one can be a \link{GRanges} or \link{GRangesList} object.

    For \code{pintersect}: one of \code{x} or \code{y} must be a \link{GRanges}
    object, the other one can be a \link{GRanges}, \link{GRangesList} or
    \link{GappedAlignments} object.

    For \code{psetdiff}: \code{x} and \code{y} can be any combination of
    \link{GRanges} and/or \link{GRangesList} objects, with the exception
    that if \code{x} is a \link{GRangesList} object then \code{y} must be
    a \link{GRangesList} too.

    In addition, for the "parallel" operations, \code{x} and \code{y} must
    be of equal length (i.e. \code{length(x) == length(y)}).
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
    For set operations: If set to TRUE, then the strand of \code{x} and
    \code{y} is set to \code{"*"} prior to any computation.

    For parallel set operations: If set to TRUE, the strand information is
    ignored in the computation and the result has the strand information of
    \code{x}.
  }
  \item{...}{
    Further arguments to be passed to or from other methods.
  }
}

\details{
  The \code{pintersect} methods involving \link{GRanges}, \link{GRangesList}
  and/or \link{GappedAlignments} objects use the triplet (sequence name, range,
  strand) to determine the element by element intersection of features,
  where a strand value of \code{"*"} is treated as occurring on both the
  \code{"+"} and \code{"-"} strand.

  The \code{psetdiff} methods involving \link{GRanges} and/or
  \link{GRangesList} objects use the triplet (sequence name, range,
  strand) to determine the element by element set difference of features,
  where a strand value of \code{"*"} is treated as occurring on both the
  \code{"+"} and \code{"-"} strand.
}

\value{
  For \code{union}, \code{intersect}, \code{setdiff}, and \code{pgap}: a
  \link{GRanges}.

  For \code{punion} and \code{pintersect}: when \code{x} or \code{y} is
  not a \link{GRanges} object, an object of the same class as this
  non-\link{GRanges} object. Otherwise, a \link{GRanges} object.

  For \code{psetdiff}: either a \link{GRanges} object when both \code{x}
  and \code{y} are \link{GRanges} objects, or a \link{GRangesList} object
  when \code{y} is a \link{GRangesList} object.
}

\author{P. Aboyoun}

\seealso{
  \link[IRanges]{IRanges-setops},
  \link{GRanges-class},
  \link{GRangesList-class},
  \link{GappedAlignments-class},
  \link{findOverlaps-methods}
}

\examples{
## ---------------------------------------------------------------------
## A. SET OPERATIONS
## ---------------------------------------------------------------------

x <- GRanges("chr1", IRanges(c(2, 9) , c(7, 19)), strand=c("+", "-"))
y <- GRanges("chr1", IRanges(5, 10), strand="-") 

union(x, y)
union(x, y, ignore.strand=TRUE)

intersect(x, y)
intersect(x, y, ignore.strand=TRUE)

setdiff(x, y)
setdiff(x, y, ignore.strand=TRUE)

## ---------------------------------------------------------------------
## B. PARALLEL SET OPERATIONS
## ---------------------------------------------------------------------

\dontrun{
punion(x, shift(x, 7))  # will fail
}
punion(x, shift(x, 7), fill.gap=TRUE)

pintersect(x, shift(x, 6))
\dontrun{
pintersect(x, shift(x, 7))  # will fail
}
pintersect(x, shift(x, 7), resolve.empty="max.start")

psetdiff(x, shift(x, 7))

## ---------------------------------------------------------------------
## C. MORE EXAMPLES
## ---------------------------------------------------------------------

## GRanges object:
gr <- GRanges(seqnames=c("chr2", "chr1", "chr1"),
              ranges=IRanges(1:3, width = 12),
              strand=Rle(strand(c("-", "*", "-"))))

## GRangesList object
gr1 <- GRanges(seqnames="chr2",
               ranges=IRanges(3, 6))
gr2 <- GRanges(seqnames=c("chr1", "chr1"),
               ranges=IRanges(c(7,13), width = 3),
               strand=c("+", "-"))
gr3 <- GRanges(seqnames=c("chr1", "chr2"),
               ranges=IRanges(c(1, 4), c(3, 9)),
               strand=c("-", "-"))
grlist <- GRangesList(gr1=gr1, gr2=gr2, gr3=gr3)

## Parallel intersection of a GRanges and a GRangesList object
pintersect(gr, grlist)
pintersect(grlist, gr)

## Parallel intersection of a GappedAlignments and a GRanges object
library(Rsamtools)  # because file ex1.bam is in this package
galn_file <- system.file("extdata", "ex1.bam", package="Rsamtools")
galn <- readGappedAlignments(galn_file)
pintersect(galn, shift(as(galn, "GRanges"), 6L))

## Parallel set difference of a GRanges and a GRangesList object
psetdiff(gr, grlist)

## Parallel set difference of two GRangesList objects
psetdiff(grlist, shift(grlist, 3))
}

\keyword{methods}
\keyword{utilities}
