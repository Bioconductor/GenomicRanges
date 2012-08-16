\name{GenomicRanges-utils}

\alias{GenomicRanges-utils}

\alias{shift,GenomicRanges-method}
\alias{narrow,GenomicRanges-method}
\alias{flank,GenomicRanges-method}
\alias{resize,GenomicRanges-method}
\alias{Ops,GenomicRanges,numeric-method}
\alias{restrict,GenomicRanges-method}
\alias{gaps,GenomicRanges-method}
\alias{range,GenomicRanges-method}
\alias{reduce,GenomicRanges-method}
\alias{isDisjoint,GenomicRanges-method}
\alias{disjoin,GenomicRanges-method}
\alias{disjointBins,GenomicRanges-method}
\alias{precede,GenomicRanges,GenomicRanges-method}
\alias{precede,GenomicRanges,missing-method}
\alias{follow,GenomicRanges,GenomicRanges-method}
\alias{follow,GenomicRanges,missing-method}
\alias{nearest,GenomicRanges,GenomicRanges-method}
\alias{nearest,GenomicRanges,missing-method}
\alias{distance,GenomicRanges,GenomicRanges-method}
\alias{distanceToNearest,GenomicRanges,GenomicRanges-method}
\alias{distanceToNearest,GenomicRanges,missing-method}

\title{Basic ranges operations on GenomicRanges objects}

\description{
  Utility functions for performing basic ranges operations on
  \link{GenomicRanges} objects (i.e. on objects that belong to the
  \link{GenomicRanges} class or one of its subclasses, this includes
  for example \link{GRanges} objects).
}

\usage{
\S4method{shift}{GenomicRanges}(x, shift=0L, use.names=TRUE)

\S4method{narrow}{GenomicRanges}(x, start=NA, end=NA, width=NA, use.names=TRUE)

\S4method{flank}{GenomicRanges}(x, width, start=TRUE, both=FALSE,
        use.names=TRUE, ignore.strand=FALSE)

\S4method{resize}{GenomicRanges}(x, width, fix="start", use.names=TRUE,
        ignore.strand=FALSE)

\S4method{restrict}{GenomicRanges}(x, start=NA, end=NA, keep.all.ranges=FALSE,
        use.names=TRUE)

\S4method{gaps}{GenomicRanges}(x, start=1L, end=seqlengths(x))

\S4method{range}{GenomicRanges}(x, ..., ignore.strand=FALSE, na.rm=FALSE)

\S4method{reduce}{GenomicRanges}(x, drop.empty.ranges=FALSE, min.gapwidth=1L,
        with.inframe.attrib=FALSE, ignore.strand=FALSE)

\S4method{isDisjoint}{GenomicRanges}(x, ignore.strand=FALSE)

\S4method{disjoin}{GenomicRanges}(x, ignore.strand=FALSE)

\S4method{disjointBins}{GenomicRanges}(x, ignore.strand=FALSE)

\S4method{precede}{GenomicRanges,GenomicRanges}(x, subject,
        select=c("arbitrary", "all"), ignore.strand=FALSE, ...)

\S4method{precede}{GenomicRanges,missing}(x, subject,
        select=c("arbitrary", "all"), ignore.strand=FALSE, ...)

\S4method{follow}{GenomicRanges,GenomicRanges}(x, subject,
        select=c("arbitrary", "all"), ignore.strand=FALSE, ...)

\S4method{follow}{GenomicRanges,missing}(x, subject,
        select=c("arbitrary", "all"), ignore.strand=FALSE, ...)

\S4method{distance}{GenomicRanges,GenomicRanges}(x, y, ignore.strand=FALSE, ...)

\S4method{distanceToNearest}{GenomicRanges,GenomicRanges}(x, subject,
        ignore.strand=FALSE, ...)
\S4method{distanceToNearest}{GenomicRanges,missing}(x, subject,
        ignore.strand=FALSE, ...)
}

\arguments{
  \item{x, subject, y}{\link{GenomicRanges} objects.}
  \item{shift, use.names, start, end, width, both, fix, keep.all.ranges, drop.empty.ranges, min.gapwidth, with.inframe.attrib, select}{
    See \code{?`\link[IRanges]{Ranges-utils}`}.
  } 
  \item{ignore.strand}{
    \code{TRUE} or \code{FALSE}. Whether the strand of the input ranges
    should be ignored or not. See details below.
  }
  \item{na.rm}{Ignored.}
  \item{...}{
    For \code{range}, additional \code{GenomicRanges} objects to consider.
    Ignored otherwise.
  }
}

\details{
\code{shift} behaves like the \code{shift} method for \link[IRanges]{Ranges}
objects. See \code{?`\link[IRanges]{Ranges-utils}`} for the details.

\code{narrow} behaves like the \code{narrow} method for \link[IRanges]{Ranges}
objects. See \code{?`\link[IRanges]{Ranges-utils}`} for the details.

\code{flank} returns an object of the same type and length as \code{x}
containing intervals of width \code{width} that flank the intervals in
\code{x}. The \code{start} argument takes a logical indicating whether
\code{x} should be flanked at the "start" (\code{TRUE}) or the "end"
(\code{FALSE}), which for \code{strand(x) != "-"} is \code{start(x)} and
\code{end(x)} respectively and for \code{strand(x) == "-"} is 
\code{end(x)} and \code{start(x)} respectively.
The \code{both} argument takes a single logical value indicating whether
the flanking region \code{width} positions extends \emph{into} the range.
If \code{both=TRUE}, the resulting range thus straddles the end point,
with \code{width} positions on either side.

\code{resize} returns an object of the same type and length as \code{x}
containing intervals that have been resized to width \code{width} based
on the \code{strand(x)} values.  Elements where \code{strand(x) == "+"}
or \code{strand(x) == "*"} are anchored at \code{start(x)} and elements
where \code{strand(x) == "-"} are anchored at the \code{end(x)}. The
\code{use.names} argument determines whether or not to keep the names on
the ranges.

\code{restrict} returns an object of the same type and length as \code{x}
containing restricted ranges for distinct seqnames. The \code{start} and
\code{end} arguments can be a named numeric vector of seqnames for the
ranges to be resticted or a numeric vector or length 1 if the restriction
operation is to be applied to all the sequences in \code{x}.
See \code{?\link[IRanges:RangesList-utils]{restrict}} for more information
about range restriction and for a description of the optional arguments.

\code{gaps} returns an object of the same type as \code{x}
containing complemented ranges for each distinct (seqname, strand) pairing.
The names (\code{names(x)}) and the columns in \code{x} are dropped.
For the start and end arguments of this gaps method, it is expected that
the user will supply a named integer vector (where the names correspond to
the appropriate seqlevels). See \code{?\link[IRanges:RangesList-utils]{gaps}}
for more information about range complements and for a description of the
optional arguments.

\code{range} returns an object of the same type as \code{x}
containing range bounds for each distinct (seqname, strand) pairing.
The names (\code{names(x)}) and the metadata columns in \code{x} are dropped.

\code{reduce} returns an object of the same type as \code{x}
containing reduced ranges for each distinct (seqname, strand) pairing.
The names (\code{names(x)}) and the metadata columns in \code{x} are dropped.
See \code{?\link[IRanges:RangesList-utils]{reduce}} for more information
about range reduction and for a description of the optional arguments.

\code{isDisjoint} returns a logical value indicating whether the ranges
in \code{x} are disjoint (i.e. non-overlapping).

\code{disjoin} returns an object of the same type as \code{x}
containing disjoint ranges for each distinct (seqname, strand) pairing.
The names (\code{names(x)}) and the metadata columns in \code{x} are dropped.

\code{disjointBins} returns bin indexes for the ranges in \code{x}, such
that ranges in the same bin do not overlap. If \code{ignore.strand=FALSE},
the two features cannot overlap if they are on different strands.

\code{precede} identifies which subject(s) the query precedes. Returns the
index of the range in \code{subject} that is directly preceded by the range
in \code{x}. When \code{ignore.strand=TRUE}, the strand for both \code{x}
and \code{subject} are set to \code{"+"}. When \code{select="arbitrary"}
an integer vector is returned with a single match per \code{x}. If no match
is found an \code{NA} is returned. When \code{select="all"} a
\link[IRanges]{Hits} object is returned with all matches for \code{x}.
If \code{x} does not have a match in \code{subject} the \code{x} is not
included in the \code{Hits} object. Overlapping ranges are excluded.
Matching by strand :
\itemize{
  \item{}{x on + strand can match to ranges on both + and 
    * strands. In the case of a tie the first range by order 
    is chosen.
  } 
  \item{}{x on - strand can match to ranges on both - and 
    * strands. In the case of a tie the first range by order 
    is chosen.
  } 
  \item{}{x on * strand can match to ranges on any of +, 
    - or * strands. In the case of a tie the first range by
    order is chosen.
  }
}

\code{follow} identifies which subject(s) the query follows. Returns the
index of the interval in \code{subject} that is directly followed by the
range in \code{x}. When \code{ignore.strand=TRUE}, both \code{x} and
\code{subject} strand are set to \code{"+"}.  When \code{select="arbitrary"}
an integer vector is returned with a single match per \code{x}. If no match
is found an \code{NA} is returned. When \code{select="all"} a
\link[IRanges]{Hits} object is returned with all matches for \code{x}.
If \code{x} does not have a match in \code{subject} the \code{x} is not
included in the \code{Hits} object. Overlapping ranges are excluded. See
the \code{precede} above for details of strand matching.

\code{nearest} does conventional nearest neighbor finding. It returns an
integer vector containing the index of the nearest neightbor range in
\code{subject} for each range in \code{x}. If there is no nearest neighbor
\code{NA} is returned. See the \code{\link[IRanges]{nearest}} man page in
the \code{IRanges} package for a description of the algorithm used.
If \code{ignore.strand=TRUE}, both \code{x} and \code{subject} strand
are set to \code{"+"}. \code{nearest} calls \code{precede} and \code{follow}
and thus the strand matching for \code{*} follows the conventions documented
under those functions. 

\code{distance} calculates the number of positions separating two features.
The value is zero if the features overlap and \code{NA} if the features
are on different sequences, or different strands (if \code{ignore.strand} is
\code{FALSE}).

\code{distanceToNearest} returns the distance for each range in \code{x}
to its nearest neighbor in the \code{subject}. If \code{ignore.strand=TRUE},
both \code{x} and \code{subject} strand are set to \code{"+"}. 
}

\author{P. Aboyoun and V. Obenchain <vobencha@fhcrc.org>}

\seealso{
  \itemize{
    \item The \link{GenomicRanges} and \link{GRanges} classes.
    \item The \link[IRanges]{Ranges} and \link[IRanges]{Hits} classes in the
          IRanges package.
    \item The \link[IRanges]{Ranges-utils} man page in the IRanges package.
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

## Intra-interval operations:
shift(gr, 1)
narrow(gr[-10], start=2, end=-2)
flank(gr, 10)
resize(gr, 10)
restrict(gr, start=3, end=7)

## Inter-interval operations:
gaps(gr, start = 1, end = 10)
range(gr)
reduce(gr)

## isDisjoint(), disjoin(), disjointBins():
isDisjoint(gr)
disjoin(gr)
disjointBins(gr)
stopifnot(all(sapply(split(gr, disjointBins(gr)), isDisjoint)))

## precede() and follow():

## query on "+" strand
query <- GRanges("A", IRanges(c(5, 20), width=1), strand="+")
subject <- GRanges("A", IRanges(rep(c(10, 15), 2), width=1),
                        strand=c("+", "+", "-", "-"))
precede(query, subject)
follow(query, subject)

## query on "+" strand
strand(query) <- "-"
precede(query, subject)
follow(query, subject)

## ties choose first in order
query <- GRanges("A", IRanges(10, width=1), c("+", "-", "*"))
subject <- GRanges("A", IRanges(c(5, 5, 5, 15, 15, 15), width=1),
                        rep(c("+", "-", "*"), 2))
precede(query, subject)
precede(query, rev(subject))
}

\keyword{utilities}
