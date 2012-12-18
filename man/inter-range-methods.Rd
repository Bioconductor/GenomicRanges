\name{inter-range-methods}

\alias{inter-range-methods}

\alias{reduce,GenomicRanges-method}
\alias{range,GenomicRanges-method}
\alias{gaps,GenomicRanges-method}
\alias{disjoin,GenomicRanges-method}
\alias{isDisjoint,GenomicRanges-method}
\alias{disjointBins,GenomicRanges-method}


\title{Inter range transformations of a GenomicRanges object}

\description{
  See \code{?`\link[IRanges]{intra-range-methods}`} and
  \code{?`\link[IRanges]{inter-range-methods}`} in the IRanges package
  for a quick introduction to intra range and inter range transformations.

  This man page documents inter range transformations of a
  \link{GenomicRanges} object (i.e. of an object that belongs to the
  \link{GenomicRanges} class or one of its subclasses, this includes
  for example \link{GRanges} objects).

  See \code{?`\link[GenomicRanges]{intra-range-methods}`} for
  intra range transformations of a GenomicRanges object.
}

\usage{
\S4method{range}{GenomicRanges}(x, ..., ignore.strand=FALSE, na.rm=FALSE)

\S4method{reduce}{GenomicRanges}(x, drop.empty.ranges=FALSE, min.gapwidth=1L,
        with.inframe.attrib=FALSE, ignore.strand=FALSE)

\S4method{gaps}{GenomicRanges}(x, start=1L, end=seqlengths(x))

\S4method{disjoin}{GenomicRanges}(x, ignore.strand=FALSE)

\S4method{isDisjoint}{GenomicRanges}(x, ignore.strand=FALSE)

\S4method{disjointBins}{GenomicRanges}(x, ignore.strand=FALSE)
}

\arguments{
  \item{x}{A \link{GenomicRanges} object.}
  \item{start, end, drop.empty.ranges, min.gapwidth, with.inframe.attrib}{
    See \code{?`\link[IRanges]{inter-range-methods}`}.
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
\code{range} returns an object of the same type as \code{x}
containing range bounds for each distinct (seqname, strand) pairing.
The names (\code{names(x)}) and the metadata columns in \code{x} are dropped.

\code{reduce} returns an object of the same type as \code{x}
containing reduced ranges for each distinct (seqname, strand) pairing.
The names (\code{names(x)}) and the metadata columns in \code{x} are dropped.
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
The names (\code{names(x)}) and the metadata columns in \code{x} are dropped.

\code{isDisjoint} returns a logical value indicating whether the ranges
in \code{x} are disjoint (i.e. non-overlapping).

\code{disjointBins} returns bin indexes for the ranges in \code{x}, such
that ranges in the same bin do not overlap. If \code{ignore.strand=FALSE},
the two features cannot overlap if they are on different strands.
}

\author{H. Pages and P. Aboyoun}

\seealso{
  \itemize{
    \item The \link{GenomicRanges} and \link{GRanges} classes.
    \item The \link[IRanges]{Ranges} class in the IRanges package.
    \item The \link[IRanges]{inter-range-methods} man page in the
          IRanges package.
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

range(gr)
reduce(gr)
gaps(gr, start = 1, end = 10)

## disjoin(), isDisjoint(), disjointBins():
disjoin(gr)
isDisjoint(gr)
stopifnot(isDisjoint(disjoin(gr)))
disjointBins(gr)
stopifnot(all(sapply(split(gr, disjointBins(gr)), isDisjoint)))
}
 
\keyword{utilities}
