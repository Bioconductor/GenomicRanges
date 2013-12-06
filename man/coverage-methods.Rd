\name{coverage-methods}
\alias{coverage-methods}

\alias{coverage,GenomicRanges-method}
\alias{coverage,GRangesList-method}


\title{Coverage of a GRanges or GRangesList object}

\description{
  \code{\link[IRanges]{coverage}} methods for \link{GRanges} and
  \link{GRangesList} objects.

  NOTE: The \code{\link[IRanges]{coverage}} generic function and methods
  for \link[IRanges]{Ranges} and \link[IRanges]{RangesList} objects
  are defined and documented in the \pkg{IRanges} package.
  Methods for \link[GenomicAlignments]{GAlignments} and
  \link[GenomicAlignments]{GAlignmentPairs} objects are defined and
  documented in the \pkg{GenomicAlignments} package.
}

\usage{
\S4method{coverage}{GenomicRanges}(x, shift=0L, width=NULL, weight=1L,
            method=c("auto", "sort", "hash"))

\S4method{coverage}{GRangesList}(x, shift=0L, width=NULL, weight=1L,
            method=c("auto", "sort", "hash"))
}

\arguments{
  \item{x}{
    A \link{GRanges} or \link{GRangesList} object.
  }
  \item{shift}{
    A numeric vector or a list-like object. If numeric, it must be parallel
    to \code{x} (recycled if necessary). If a list-like object, it must
    have 1 list element per seqlevel in \code{x}, and its names must be
    exactly \code{seqlevels(x)}.

    Alternatively, \code{shift} can also be specified as a single string
    naming a metadata column in \code{x} (i.e. a column in \code{mcols(x)})
    to be used as the \code{shift} vector.

    See \code{?\link[IRanges]{coverage}} in the \pkg{IRanges} package
    for more information about this argument.
  }
  \item{width}{
    Either \code{NULL} (the default), or an integer vector.
    If \code{NULL}, it is replaced with \code{seqlengths(x)}.
    Otherwise, the vector must have the length and names of
    \code{seqlengths(x)} and contain NAs or non-negative integers.

    See \code{?\link[IRanges]{coverage}} in the \pkg{IRanges} package
    for more information about this argument.
  }
  \item{weight}{
    A numeric vector or a list-like object. If numeric, it must be parallel
    to \code{x} (recycled if necessary). If a list-like object, it must
    have 1 list element per seqlevel in \code{x}, and its names must be
    exactly \code{seqlevels(x)}.

    Alternatively, \code{weight} can also be specified as a single string
    naming a metadata column in \code{x} (i.e. a column in \code{mcols(x)})
    to be used as the \code{weight} vector.

    See \code{?\link[IRanges]{coverage}} in the \pkg{IRanges} package
    for more information about this argument.
  }
  \item{method}{
    See \code{?\link[IRanges]{coverage}} in the \pkg{IRanges} package
    for a description of this argument.
  }
}

\details{
  When \code{x} is a \link{GRangesList} object, \code{coverage(x, ...)}
  is equivalent to \code{coverage(unlist(x), ...)}.
}

\value{
  A named \link[IRanges]{RleList} object with one coverage vector per
  seqlevel in \code{x}.
}

\author{H. Pages and P. Aboyoun}

\seealso{
  \itemize{
    \item \code{\link[IRanges]{coverage}} in the \pkg{IRanges} package.

    \item \link[GenomicAlignments]{coverage-methods} in the
          \pkg{GenomicAlignments} package.

    \item \link[IRanges]{RleList} objects in the \pkg{IRanges} package.

    \item \link{GRanges} and \link{GRangesList} objects.
  }
}

\examples{
## Coverage of a GRanges object:
gr <- GRanges(
        seqnames=Rle(c("chr1", "chr2", "chr1", "chr3"), c(1, 3, 2, 4)),
        ranges=IRanges(1:10, end=10),
        strand=Rle(strand(c("-", "+", "*", "+", "-")), c(1, 2, 2, 3, 2)),
        seqlengths=c(chr1=11, chr2=12, chr3=13))
cvg <- coverage(gr)
pcvg <- coverage(gr[strand(gr) == "+"])
mcvg <- coverage(gr[strand(gr) == "-"])
scvg <- coverage(gr[strand(gr) == "*"])
stopifnot(identical(pcvg + mcvg + scvg, cvg))

## Coverage of a GRangesList object:
gr1 <- GRanges(seqnames="chr2",
               ranges=IRanges(3, 6),
               strand = "+")
gr2 <- GRanges(seqnames=c("chr1", "chr1"),
               ranges=IRanges(c(7,13), width=3),
               strand=c("+", "-"))
gr3 <- GRanges(seqnames=c("chr1", "chr2"),
               ranges=IRanges(c(1, 4), c(3, 9)),
               strand=c("-", "-"))
grl <- GRangesList(gr1=gr1, gr2=gr2, gr3=gr3)
stopifnot(identical(coverage(grl), coverage(unlist(grl))))
}

\keyword{methods}
\keyword{utilities}
