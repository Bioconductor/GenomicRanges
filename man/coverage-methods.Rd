\name{coverage-methods}
\alias{coverage-methods}

\alias{coverage,GenomicRanges-method}
\alias{coverage,GRangesList-method}
\alias{coverage,GAlignments-method}
\alias{coverage,GAlignmentPairs-method}


\title{Coverage of a GRanges, GRangesList, GAlignments, or
       GAlignmentPairs object}

\description{
  \code{\link[IRanges]{coverage}} methods for \link{GRanges},
  \link{GRangesList}, \link{GAlignments}, and \link{GAlignmentPairs}
  objects.

  NOTE: The \code{\link[IRanges]{coverage}} generic function and methods
  for \link[IRanges]{Ranges} and \link[IRanges]{RangesList} objects
  are defined and documented in the \pkg{IRanges} package.
}

\usage{
\S4method{coverage}{GenomicRanges}(x, shift=0L, width=NULL, weight=1L,
            method=c("auto", "sort", "hash"))

\S4method{coverage}{GRangesList}(x, shift=0L, width=NULL, weight=1L,
            method=c("auto", "sort", "hash"))

\S4method{coverage}{GAlignments}(x, shift=0L, width=NULL, weight=1L,
            method=c("auto", "sort", "hash"), drop.D.ranges=FALSE)

\S4method{coverage}{GAlignmentPairs}(x, shift=0L, width=NULL, weight=1L,
            method=c("auto", "sort", "hash"), drop.D.ranges=FALSE)
}

\arguments{
  \item{x}{
    A \link{GRanges}, \link{GRangesList}, \link{GAlignments},
    or \link{GAlignmentPairs} object.
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
  \item{drop.D.ranges}{
    Whether the coverage calculation should ignore ranges corresponding
    to D (deletion) in the CIGAR string.
  }
}

\details{
  When \code{x} is a \link{GRangesList} object, \code{coverage(x, ...)}
  is equivalent to \code{coverage(unlist(x), ...)}.

  When \code{x} is a \link{GAlignments} or \link{GAlignmentPairs}
  object, \code{coverage(x, ...)} is equivalent to
  \code{coverage(grglist(x), ...)}.
}

\value{
  A named \link[IRanges]{RleList} object with one coverage vector per
  seqlevel in \code{x}.
}

\author{H. Pages and P. Aboyoun}

\seealso{
  \itemize{
    \item \code{\link[IRanges]{coverage}}.
    \item \link[IRanges]{RleList-class}.
    \item \link{GRanges-class}.
    \item \link{GRangesList-class}.
    \item \link{GAlignments-class}.
    \item \link{GAlignmentPairs-class}.
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

## Coverage of a GAlignments or GAlignmentPairs object:
library(Rsamtools)  # because file ex1.bam is in this package
ex1_file <- system.file("extdata", "ex1.bam", package="Rsamtools")
galn <- readGAlignments(ex1_file)
stopifnot(identical(coverage(galn), coverage(as(galn, "GRangesList"))))
galp <- readGAlignmentPairs(ex1_file)
stopifnot(identical(coverage(galp), coverage(as(galp, "GRangesList"))))
}

\keyword{methods}
\keyword{utilities}
