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

\alias{psetdiff,GRanges,GRanges-method}
\alias{psetdiff,GRanges,GRangesList-method}
\alias{psetdiff,GRangesList,GRangesList-method}

\alias{pgap,GRanges,GRanges-method}


\title{Set operations on GRanges and GRangesList objects}

\description{
  Performs set operations on \link{GRanges} and \link{GRangesList} objects.
}

\usage{
## Set operations
\S4method{union}{GRanges,GRanges}(x, y, ignore.strand=FALSE, ...)
\S4method{intersect}{GRanges,GRanges}(x, y, ignore.strand=FALSE, ...)
\S4method{setdiff}{GRanges,GRanges}(x, y, ignore.strand=FALSE, ...)

## Parallel set operations
\S4method{punion}{GRanges,GRanges}(x, y, fill.gap=FALSE, ignore.strand=FALSE, ...)
\S4method{pintersect}{GRanges,GRanges}(x, y, drop.nohit.ranges=FALSE,
           ignore.strand=FALSE, strict.strand=FALSE)
\S4method{psetdiff}{GRanges,GRanges}(x, y, ignore.strand=FALSE, ...)
}

\arguments{
  \item{x, y}{
    For \code{union}, \code{intersect}, \code{setdiff}, \code{pgap}:
    \code{x} and \code{y} must both be \link{GRanges} objects.

    For \code{punion}: one of \code{x} or \code{y} must be a \link{GRanges}
    object, the other one can be a \link{GRanges} or \link{GRangesList} object.

    For \code{pintersect}: \code{x} and \code{y} can be any combination of
    \link{GRanges} and/or \link{GRangesList} objects.

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
  \item{...}{
    Further arguments to be passed to or from other methods.
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
  For \code{union}, \code{intersect}, \code{setdiff}, and \code{pgap}: a
  \link{GRanges} object.

  For \code{punion} and \code{pintersect}: when \code{x} or \code{y} is
  not a \link{GRanges} object, an object of the same class as this
  non-\link{GRanges} object. Otherwise, a \link{GRanges} object.

  For \code{psetdiff}: either a \link{GRanges} object when both \code{x}
  and \code{y} are \link{GRanges} objects, or a \link{GRangesList} object
  when \code{y} is a \link{GRangesList} object.
}

\author{P. Aboyoun and H. Pages}

\seealso{
  \link[IRanges]{setops-methods},
  \link{GRanges-class},
  \link{GRangesList-class},
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

## For a fast 'mendoapply(intersect, grlist, as(gr, "GRangesList"))'
## call pintersect() with 'strict.strand=TRUE' and call reduce() on
## the result with 'drop.empty.ranges=TRUE':
reduce(pintersect(grlist, gr, strict.strand=TRUE),
       drop.empty.ranges=TRUE)

## Parallel set difference of a GRanges and a GRangesList object
psetdiff(gr, grlist)

## Parallel set difference of two GRangesList objects
psetdiff(grlist, shift(grlist, 3))
}

\keyword{methods}
\keyword{utilities}
