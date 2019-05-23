\name{inter-range-methods}

\alias{inter-range-methods}

\alias{range}
\alias{range,GenomicRanges-method}
\alias{range,StitchedGPos-method}
\alias{range,GRangesList-method}

\alias{reduce}
\alias{reduce,GenomicRanges-method}
\alias{reduce,GRangesList-method}

\alias{gaps}
\alias{gaps,GenomicRanges-method}

\alias{disjoin}
\alias{disjoin,GenomicRanges-method}
\alias{disjoin,GRangesList-method}

\alias{isDisjoint}
\alias{isDisjoint,GenomicRanges-method}
\alias{isDisjoint,StitchedGPos-method}
\alias{isDisjoint,GRangesList-method}

\alias{disjointBins}
\alias{disjointBins,GenomicRanges-method}

\title{Inter range transformations of a GRanges or GRangesList object}

\description{
  This man page documents \emph{inter range transformations} of a
  \link{GenomicRanges} object (i.e. of an object that belongs to the
  \link{GenomicRanges} class or one of its subclasses, this includes
  for example \link{GRanges} objects), or a \link{GRangesList} object.

  See \code{?`\link[IRanges]{intra-range-methods}`} and
  \code{?`\link[IRanges]{inter-range-methods}`} in the \pkg{IRanges}
  package for a quick introduction to \emph{intra range} and \emph{inter
  range transformations}.

  See \code{?`\link[GenomicRanges]{intra-range-methods}`} for
  \emph{intra range transformations} of a \link{GenomicRanges} object or
  \link{GRangesList} object.
}

\usage{
\S4method{range}{GenomicRanges}(x, ..., with.revmap=FALSE, ignore.strand=FALSE, na.rm=FALSE)
\S4method{range}{GRangesList}(x, ..., with.revmap=FALSE, ignore.strand=FALSE, na.rm=FALSE)

\S4method{reduce}{GenomicRanges}(x, drop.empty.ranges=FALSE, min.gapwidth=1L, with.revmap=FALSE,
          with.inframe.attrib=FALSE, ignore.strand=FALSE)
\S4method{reduce}{GRangesList}(x, drop.empty.ranges=FALSE, min.gapwidth=1L, with.revmap=FALSE,
          with.inframe.attrib=FALSE, ignore.strand=FALSE)

\S4method{gaps}{GenomicRanges}(x, start=1L, end=seqlengths(x))

\S4method{disjoin}{GenomicRanges}(x, with.revmap=FALSE, ignore.strand=FALSE)
\S4method{disjoin}{GRangesList}(x, with.revmap=FALSE, ignore.strand=FALSE)

\S4method{isDisjoint}{GenomicRanges}(x, ignore.strand=FALSE)
\S4method{isDisjoint}{GRangesList}(x, ignore.strand=FALSE)

\S4method{disjointBins}{GenomicRanges}(x, ignore.strand=FALSE)
}

\arguments{
  \item{x}{A \link{GenomicRanges} or \link{GenomicRangesList} object.}
  \item{drop.empty.ranges, min.gapwidth, with.revmap, with.inframe.attrib, start, end}{
    See \code{?`\link[IRanges]{inter-range-methods}`} in the \pkg{IRanges}
    package.
  } 
  \item{ignore.strand}{
    \code{TRUE} or \code{FALSE}. Whether the strand of the input ranges
    should be ignored or not. See details below.
  }
  \item{...}{
    For \code{range}, additional \code{GenomicRanges} objects to consider.
    Ignored otherwise.
  }
  \item{na.rm}{Ignored.}
}

\details{
  \subsection{On a GRanges object}{
    \code{range} returns an object of the same type as \code{x}
    containing range bounds for each distinct (seqname, strand) pairing.
    The names (\code{names(x)}) and the metadata columns in \code{x} are
    dropped.

    \code{reduce} returns an object of the same type as \code{x}
    containing reduced ranges for each distinct (seqname, strand) pairing.
    The names (\code{names(x)}) and the metadata columns in \code{x} are
    dropped.
    See \code{?\link[IRanges]{reduce}} for more information about range
    reduction and for a description of the optional arguments.

    \code{gaps} returns an object of the same type as \code{x}
    containing complemented ranges for each distinct (seqname, strand) pairing.
    The names (\code{names(x)}) and the columns in \code{x} are dropped.
    For the start and end arguments of this gaps method, it is expected that
    the user will supply a named integer vector (where the names correspond to
    the appropriate seqlevels). See \code{?\link[IRanges]{gaps}} for more
    information about range complements and for a description of the optional
    arguments.

    \code{disjoin} returns an object of the same type as \code{x}
    containing disjoint ranges for each distinct (seqname, strand) pairing.
    The names (\code{names(x)}) and the metadata columns in \code{x} are
    dropped. If \code{with.revmap=TRUE}, a metadata column that maps the
    ouput ranges to the input ranges is added to the returned object.
    See \code{?\link[IRanges]{disjoin}} for more information.

    \code{isDisjoint} returns a logical value indicating whether the ranges
    in \code{x} are disjoint (i.e. non-overlapping).

    \code{disjointBins} returns bin indexes for the ranges in \code{x}, such
    that ranges in the same bin do not overlap. If \code{ignore.strand=FALSE},
    the two features cannot overlap if they are on different strands.
  }
  \subsection{On a GRangesList/GenomicRangesList object}{
    When they are supported on GRangesList object \code{x}, the above inter
    range transformations will apply the transformation to each of the list
    elements in \code{x} and return a list-like object \emph{parallel} to
    \code{x} (i.e. with 1 list element per list element in \code{x}).
    If \code{x} has names on it, they're propagated to the returned object.
  }
}

\author{H. Pagès and P. Aboyoun}

\seealso{
  \itemize{
    \item The \link{GenomicRanges} and \link{GRanges} classes.

    \item The \link[IRanges]{IntegerRanges} class in the \pkg{IRanges} package.

    \item The \link[IRanges]{inter-range-methods} man page in the
          \pkg{IRanges} package.

    \item \link{GenomicRanges-comparison} for comparing and ordering
          genomic ranges.

    \item \code{\link[S4Vectors]{endoapply}} in the \pkg{S4Vectors} package.
  }
}

\examples{
gr <- GRanges(
        seqnames=Rle(paste("chr", c(1, 2, 1, 3), sep=""), c(1, 3, 2, 4)),
        ranges=IRanges(1:10, width=10:1, names=letters[1:10]),
        strand=Rle(strand(c("-", "+", "*", "+", "-")), c(1, 2, 2, 3, 2)),
        score=1:10,
        GC=seq(1, 0, length=10)
      )
gr

gr1 <- GRanges(seqnames="chr2", ranges=IRanges(3, 6),
               strand="+", score=5L, GC=0.45)
gr2 <- GRanges(seqnames="chr1",
               ranges=IRanges(c(10, 7, 19), width=5),
               strand=c("+", "-", "+"), score=3:5, GC=c(0.3, 0.5, 0.66))
gr3 <- GRanges(seqnames=c("chr1", "chr2"),
               ranges=IRanges(c(1, 4), c(3, 9)),
               strand=c("-", "-"), score=c(6L, 2L), GC=c(0.4, 0.1))
grl <- GRangesList(gr1=gr1, gr2=gr2, gr3=gr3)
grl

## ---------------------------------------------------------------------
## range()
## ---------------------------------------------------------------------

## On a GRanges object:
range(gr)
range(gr, with.revmap=TRUE)

## On a GRangesList object:
range(grl)
range(grl, ignore.strand=TRUE)
range(grl, with.revmap=TRUE, ignore.strand=TRUE)

# ---------------------------------------------------------------------
## reduce()
## ---------------------------------------------------------------------
reduce(gr)

gr2 <- reduce(gr, with.revmap=TRUE)
revmap <- mcols(gr2)$revmap  # an IntegerList

## Use the mapping from reduced to original ranges to group the original
## ranges by reduced range:
relist(gr[unlist(revmap)], revmap)

## Or use it to split the DataFrame of original metadata columns by
## reduced range:
relist(mcols(gr)[unlist(revmap), ], revmap)  # a SplitDataFrameList

## [For advanced users] Use this reverse mapping to compare the reduced
## ranges with the ranges they originate from:
expanded_gr2 <- rep(gr2, elementNROWS(revmap))
reordered_gr <- gr[unlist(revmap)]
codes <- pcompare(expanded_gr2, reordered_gr)
## All the codes should translate to "d", "e", "g", or "h" (the 4 letters
## indicating that the range on the left contains the range on the right):
alphacodes <- rangeComparisonCodeToLetter(pcompare(expanded_gr2, reordered_gr))
stopifnot(all(alphacodes \%in\% c("d", "e", "g", "h")))

## On a big GRanges object with a lot of seqlevels:
mcols(gr) <- NULL
biggr <- c(gr, GRanges("chr1", IRanges(c(4, 1), c(5, 2)), strand="+"))
seqlevels(biggr) <- paste0("chr", 1:2000)
biggr <- rep(biggr, 25000)
set.seed(33)
seqnames(biggr) <- sample(factor(seqlevels(biggr), levels=seqlevels(biggr)),
                          length(biggr), replace=TRUE)

biggr2 <- reduce(biggr, with.revmap=TRUE)
revmap <- mcols(biggr2)$revmap
expanded_biggr2 <- rep(biggr2, elementNROWS(revmap))
reordered_biggr <- biggr[unlist(revmap)]
codes <- pcompare(expanded_biggr2, reordered_biggr)
alphacodes <- rangeComparisonCodeToLetter(pcompare(expanded_biggr2,
                                                  reordered_biggr))
stopifnot(all(alphacodes \%in\% c("d", "e", "g", "h")))
table(alphacodes)

## On a GRangesList object:
reduce(grl)  # Doesn't really reduce anything but note the reordering
             # of the inner elements in the 2nd and 3rd list elements:
             # the ranges are reordered by sequence name first (which
             # should appear in the same order as in 'seqlevels(grl)'),
             # and then by strand.
reduce(grl, ignore.strand=TRUE)  # 2nd list element got reduced

## ---------------------------------------------------------------------
## gaps()
## ---------------------------------------------------------------------
gaps(gr, start=1, end=10)

## ---------------------------------------------------------------------
## disjoin(), isDisjoint(), disjointBins()
## ---------------------------------------------------------------------
disjoin(gr)
disjoin(gr, with.revmap=TRUE)
disjoin(gr, with.revmap=TRUE, ignore.strand=TRUE)
isDisjoint(gr)
stopifnot(isDisjoint(disjoin(gr)))
disjointBins(gr)
stopifnot(all(sapply(split(gr, disjointBins(gr)), isDisjoint)))

## On a GRangesList object:
disjoin(grl)  # doesn't really disjoin anything but note the reordering
disjoin(grl, with.revmap=TRUE)
}
 
\keyword{utilities}
