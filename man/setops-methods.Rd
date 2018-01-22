\name{setops-methods}
\alias{setops-methods}

\alias{union}
\alias{union,GenomicRanges,GenomicRanges-method}
\alias{union,GRangesList,GRangesList-method}
\alias{union,GenomicRanges,Vector-method}
\alias{union,Vector,GenomicRanges-method}

\alias{intersect}
\alias{intersect,GenomicRanges,GenomicRanges-method}
\alias{intersect,GRangesList,GRangesList-method}
\alias{intersect,GenomicRanges,Vector-method}
\alias{intersect,Vector,GenomicRanges-method}

\alias{setdiff}
\alias{setdiff,GenomicRanges,GenomicRanges-method}
\alias{setdiff,GRangesList,GRangesList-method}
\alias{setdiff,GenomicRanges,Vector-method}
\alias{setdiff,Vector,GenomicRanges-method}

\alias{punion}
\alias{punion,GRanges,GRanges-method}
\alias{punion,GRanges,GRangesList-method}
\alias{punion,GRangesList,GRanges-method}

\alias{pintersect}
\alias{pintersect,GRanges,GRanges-method}
\alias{pintersect,GRanges,GRangesList-method}
\alias{pintersect,GRangesList,GRanges-method}

\alias{psetdiff}
\alias{psetdiff,GRanges,GRanges-method}
\alias{psetdiff,GRanges,GRangesList-method}

\alias{pgap}
\alias{pgap,GRanges,GRanges-method}


\title{Set operations on genomic ranges}

\description{
  Performs set operations on \link{GRanges} and \link{GRangesList} objects.

  NOTE: The \code{\link[IRanges]{punion}}, \code{\link[IRanges]{pintersect}},
  \code{\link[IRanges]{psetdiff}}, and \code{\link[IRanges]{pgap}} generic
  functions and methods for \link[IRanges]{IntegerRanges} objects are defined
  and documented in the \pkg{IRanges} package.
}

\usage{
## Vector-wise set operations
## --------------------------

\S4method{union}{GenomicRanges,GenomicRanges}(x, y, ignore.strand=FALSE)

\S4method{intersect}{GenomicRanges,GenomicRanges}(x, y, ignore.strand=FALSE)

\S4method{setdiff}{GenomicRanges,GenomicRanges}(x, y, ignore.strand=FALSE)

## Element-wise (aka "parallel") set operations
## --------------------------------------------

\S4method{punion}{GRanges,GRanges}(x, y, fill.gap=FALSE, ignore.strand=FALSE)

\S4method{pintersect}{GRanges,GRanges}(x, y, drop.nohit.ranges=FALSE,
           ignore.strand=FALSE, strict.strand=FALSE)

\S4method{psetdiff}{GRanges,GRanges}(x, y, ignore.strand=FALSE)
}

\arguments{
  \item{x, y}{
    For \code{union}, \code{intersect}, and \code{setdiff}:
    2 \link{GenomicRanges} objects or 2 \link{GRangesList} objects.

    For \code{punion} and \code{pintersect}: 2 \link{GRanges} objects, or
    1 \link{GRanges} object and 1 \link{GRangesList} object.

    For \code{psetdiff}: \code{x} must be a \link{GRanges} object and
    \code{y} can be a \link{GRanges} or \link{GRangesList} object.

    For \code{pgap}: 2 \link{GRanges} objects.

    In addition, for the \emph{parallel} operations, \code{x} and \code{y}
    must be of equal length (i.e. \code{length(x) == length(y)}).
  }
  \item{fill.gap}{
    Logical indicating whether or not to force a union by using the rule
    \code{start = min(start(x), start(y)), end = max(end(x), end(y))}.
  }
  \item{ignore.strand}{
    For set operations: If set to TRUE, then the strand of \code{x} and
    \code{y} is set to \code{"*"} prior to any computation.

    For parallel set operations: If set to TRUE, the strand information is
    ignored in the computation and the result has the strand information of
    \code{x}.
  }
  \item{drop.nohit.ranges}{
    If TRUE then elements in \code{x} that don't intersect with their
    corresponding element in \code{y} are removed from the result (so
    the returned object is no more parallel to the input).

    If FALSE (the default) then nothing is removed and a \code{hit} metadata
    column is added to the returned object to indicate elements in \code{x}
    that intersect with the corresponding element in \code{y}. For those
    that don't, the reported intersection is a zero-width range that has
    the same start as \code{x}.
  }
  \item{strict.strand}{
    If set to FALSE (the default), features on the \code{"*"} strand are
    treated as occurring on both the \code{"+"} and \code{"-"} strand.
    If set to TRUE, the strand of intersecting elements must be strictly
    the same.
  }
}

\details{
  The \code{pintersect} methods involving \link{GRanges} and/or
  \link{GRangesList} objects use the triplet (sequence name, range, strand)
  to determine the element by element intersection of features, where a
  strand value of \code{"*"} is treated as occurring on both the \code{"+"}
  and \code{"-"} strand (unless \code{strict.strand} is set to TRUE, in
  which case the strand of intersecting elements must be strictly the same).

  The \code{psetdiff} methods involving \link{GRanges} and/or
  \link{GRangesList} objects use the triplet (sequence name, range,
  strand) to determine the element by element set difference of features,
  where a strand value of \code{"*"} is treated as occurring on both the
  \code{"+"} and \code{"-"} strand.
}

\value{
  For \code{union}, \code{intersect}, and \code{setdiff}: a \link{GRanges}
  object if \code{x} and \code{y} are \link{GenomicRanges} objects,
  and a \link{GRangesList} object if they are \link{GRangesList} objects.

  For \code{punion} and \code{pintersect}: when \code{x} or \code{y} is
  not a \link{GRanges} object, an object of the same class as this
  non-\link{GRanges} object. Otherwise, a \link{GRanges} object.

  For \code{psetdiff}: either a \link{GRanges} object when both \code{x}
  and \code{y} are \link{GRanges} objects, or a \link{GRangesList} object
  when \code{y} is a \link{GRangesList} object.

  For \code{pgap}: a \link{GRanges} object.
}

\author{P. Aboyoun and H. Pagès}

\seealso{
  \itemize{
    \item \link[IRanges]{setops-methods} in the \pkg{IRanges} package for set
          operations on \link[IRanges]{IntegerRanges} and
          \link[IRanges]{IntegerRangesList} objects.

    \item \link[GenomicRanges]{findOverlaps-methods} for finding/counting
          overlapping genomic ranges.

    \item \link[GenomicRanges]{intra-range-methods} and
          \link[GenomicRanges]{inter-range-methods} for \emph{intra range}
          and \emph{inter range} transformations of a \link{GRanges} object.

    \item \link{GRanges} and \link{GRangesList} objects.

    \item \code{\link[S4Vectors]{mendoapply}} in the \pkg{S4Vectors} package.
  }
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

## With 2 GRangesList objects:
gr1 <- GRanges(seqnames="chr2",
               ranges=IRanges(3, 6))
gr2 <- GRanges(seqnames=c("chr1", "chr1"),
               ranges=IRanges(c(7,13), width = 3),
               strand=c("+", "-"))
gr3 <- GRanges(seqnames=c("chr1", "chr2"),
               ranges=IRanges(c(1, 4), c(3, 9)),
               strand=c("-", "-"))
grlist <- GRangesList(gr1=gr1, gr2=gr2, gr3=gr3)

union(grlist, shift(grlist, 3))
intersect(grlist, shift(grlist, 3))
setdiff(grlist, shift(grlist, 3))

## Sanity checks:
grlist2 <- shift(grlist, 3)
stopifnot(identical(
    union(grlist, grlist2),
    mendoapply(union, grlist, grlist2)
))
stopifnot(identical(
    intersect(grlist, grlist2),
    mendoapply(intersect, grlist, grlist2)
))
stopifnot(identical(
    setdiff(grlist, grlist2),
    mendoapply(setdiff, grlist, grlist2)
))

## ---------------------------------------------------------------------
## B. PARALLEL SET OPERATIONS
## ---------------------------------------------------------------------

punion(x, shift(x, 6))
\dontrun{
punion(x, shift(x, 7))  # will fail
}
punion(x, shift(x, 7), fill.gap=TRUE)

pintersect(x, shift(x, 6))
pintersect(x, shift(x, 7))

psetdiff(x, shift(x, 7))

## ---------------------------------------------------------------------
## C. MORE EXAMPLES
## ---------------------------------------------------------------------

## GRanges object:
gr <- GRanges(seqnames=c("chr2", "chr1", "chr1"),
              ranges=IRanges(1:3, width = 12),
              strand=Rle(strand(c("-", "*", "-"))))

## Parallel intersection of a GRanges and a GRangesList object
pintersect(gr, grlist)
pintersect(grlist, gr)

## For a fast 'mendoapply(intersect, grlist, as(gr, "GRangesList"))'
## call pintersect() with 'strict.strand=TRUE' and call reduce() on
## the result with 'drop.empty.ranges=TRUE':
reduce(pintersect(grlist, gr, strict.strand=TRUE),
       drop.empty.ranges=TRUE)

## Parallel set difference of a GRanges and a GRangesList object
psetdiff(gr, grlist)
}

\keyword{methods}
\keyword{utilities}
