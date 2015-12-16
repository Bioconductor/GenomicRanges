\name{findOverlaps-methods}
\alias{findOverlaps-methods}

\alias{findOverlaps}
\alias{findOverlaps,GenomicRanges,GenomicRanges-method}
\alias{findOverlaps,GRangesList,GenomicRanges-method}
\alias{findOverlaps,GenomicRanges,GRangesList-method}
\alias{findOverlaps,GRangesList,GRangesList-method}
\alias{findOverlaps,RangesList,GenomicRanges-method}
\alias{findOverlaps,RangesList,GRangesList-method}
\alias{findOverlaps,GenomicRanges,RangesList-method}
\alias{findOverlaps,GRangesList,RangesList-method}
\alias{findOverlaps,RangedData,GenomicRanges-method}
\alias{findOverlaps,RangedData,GRangesList-method}
\alias{findOverlaps,GenomicRanges,RangedData-method}
\alias{findOverlaps,GRangesList,RangedData-method}

\alias{countOverlaps}
\alias{countOverlaps,GenomicRanges,GenomicRanges-method}
\alias{countOverlaps,GenomicRanges,Vector-method}
\alias{countOverlaps,Vector,GenomicRanges-method}
\alias{countOverlaps,GRangesList,Vector-method}
\alias{countOverlaps,Vector,GRangesList-method}
\alias{countOverlaps,GRangesList,GRangesList-method}
\alias{countOverlaps,GRanges,GRangesList-method}
\alias{countOverlaps,GRangesList,GRanges-method}

\alias{overlapsAny}
\alias{overlapsAny,GenomicRanges,GenomicRanges-method}
\alias{overlapsAny,GRangesList,GenomicRanges-method}
\alias{overlapsAny,GenomicRanges,GRangesList-method}
\alias{overlapsAny,GRangesList,GRangesList-method}
\alias{overlapsAny,RangesList,GenomicRanges-method}
\alias{overlapsAny,RangesList,GRangesList-method}
\alias{overlapsAny,GenomicRanges,RangesList-method}
\alias{overlapsAny,GRangesList,RangesList-method}
\alias{overlapsAny,RangedData,GenomicRanges-method}
\alias{overlapsAny,RangedData,GRangesList-method}
\alias{overlapsAny,GenomicRanges,RangedData-method}
\alias{overlapsAny,GRangesList,RangedData-method}

\alias{subsetByOverlaps}
\alias{subsetByOverlaps,GenomicRanges,GenomicRanges-method}
\alias{subsetByOverlaps,GRangesList,GenomicRanges-method}
\alias{subsetByOverlaps,GenomicRanges,GRangesList-method}
\alias{subsetByOverlaps,GRangesList,GRangesList-method}
\alias{subsetByOverlaps,RangesList,GenomicRanges-method}
\alias{subsetByOverlaps,RangesList,GRangesList-method}
\alias{subsetByOverlaps,GenomicRanges,RangesList-method}
\alias{subsetByOverlaps,GRangesList,RangesList-method}
\alias{subsetByOverlaps,RangedData,GenomicRanges-method}
\alias{subsetByOverlaps,RangedData,GRangesList-method}
\alias{subsetByOverlaps,GenomicRanges,RangedData-method}
\alias{subsetByOverlaps,GRangesList,RangedData-method}


\title{Finding overlapping genomic ranges}

\description{
  Various methods for finding/counting overlaps between objects containing
  genomic ranges. This man page describes the methods that operate on
  \link{GenomicRanges} or \link{GRangesList} objects.

  NOTE: The \code{\link[IRanges]{findOverlaps}} generic function and
  methods for \link[IRanges]{Ranges} and \link[IRanges]{RangesList}
  objects are defined and documented in the \pkg{IRanges} package.
  The methods for \link[GenomicAlignments]{GAlignments},
  \link[GenomicAlignments]{GAlignmentPairs}, and
  \link[GenomicAlignments]{GAlignmentsList} objects are defined and
  documented in the \pkg{GenomicAlignments} package.
}

\usage{
\S4method{findOverlaps}{GenomicRanges,GenomicRanges}(query, subject,
    maxgap=0L, minoverlap=1L,
    type=c("any", "start", "end", "within", "equal"),
    select=c("all", "first", "last", "arbitrary"),
    ignore.strand=FALSE)

\S4method{countOverlaps}{GenomicRanges,GenomicRanges}(query, subject,
    maxgap=0L, minoverlap=1L,
    type=c("any", "start", "end", "within", "equal"),
    ignore.strand=FALSE)

\S4method{overlapsAny}{GenomicRanges,GenomicRanges}(query, subject,
    maxgap=0L, minoverlap=1L,
    type=c("any", "start", "end", "within", "equal"),
    ignore.strand=FALSE)

\S4method{subsetByOverlaps}{GenomicRanges,GenomicRanges}(query, subject,
    maxgap=0L, minoverlap=1L,
    type=c("any", "start", "end", "within", "equal"),
    ignore.strand=FALSE)
}

\arguments{
  \item{query, subject}{
    A \link{GRanges} or \link{GRangesList} object.
    \link[IRanges]{RangesList} and \link[IRanges]{RangedData} are 
    also accepted for one of \code{query} or \code{subject}.
  }
  \item{maxgap, minoverlap, type}{
    See \code{\link[IRanges]{findOverlaps}} in the \pkg{IRanges} package
    for a description of these arguments.
  }
  \item{select}{
    When \code{select} is \code{"all"} (the default), the results are
    returned as a \link{Hits} object.
    Otherwise the returned value is an integer vector parallel to \code{query}
    (i.e. same length) containing the first, last, or arbitrary overlapping
    interval in \code{subject}, with \code{NA} indicating intervals that did
    not overlap any intervals in \code{subject}.
  }
  \item{ignore.strand}{
    When set to \code{TRUE}, the strand information is ignored in the
    overlap calculations.
  }
}

\details{
  When the query and the subject are \link{GRanges} or
  \link{GRangesList} objects, \code{findOverlaps} uses the triplet
  (sequence name, range, strand) to determine which features (see
  paragraph below for the definition of feature) from the \code{query}
  overlap which features in the \code{subject}, where a strand value
  of \code{"*"} is treated as occurring on both the \code{"+"} and
  \code{"-"} strand.
  An overlap is recorded when a feature in the \code{query} and a feature
  in the \code{subject} have the same sequence name, have a compatible
  pairing of strands (e.g. \code{"+"}/\code{"+"}, \code{"-"}/\code{"-"},
  \code{"*"}/\code{"+"}, \code{"*"}/\code{"-"}, etc.), and satisfy the
  interval overlap requirements. Strand is taken as \code{"*"}
  for \code{RangedData} and \code{RangesList}.

  In the context of \code{findOverlaps}, a feature is a collection of
  ranges that are treated as a single entity. For \link{GRanges} objects,
  a feature is a single range; while for \link{GRangesList} objects,
  a feature is a list element containing a set of ranges. In the results,
  the features are referred to by number, which run from 1 to
  \code{length(query)}/\code{length(subject)}.
}

\value{
  For \code{findOverlaps} either a \link[S4Vectors]{Hits} object when
  \code{select="all"} or an integer vector otherwise. 

  For \code{countOverlaps} an integer vector containing the tabulated
  query overlap hits.

  For \code{overlapsAny} a logical vector of length equal to the number of
  ranges in \code{query} indicating those that overlap any of the ranges
  in \code{subject}.

  For \code{subsetByOverlaps} an object of the same class as \code{query}
  containing the subset that overlapped at least one entity in \code{subject}.

  For \code{RangedData} and \code{RangesList}, with the exception of
  \code{subsetByOverlaps}, the results align to the unlisted
  form of the object. This turns out to be fairly convenient for
  \code{RangedData} (not so much for \code{RangesList}, but something
  has to give).
}

\author{P. Aboyoun, S. Falcon, M. Lawrence, and H. Pagès}

\seealso{
  \itemize{
    \item The \link{Hits} class for representing a set of hits between 2
          vector-like objects.

    \item The \code{\link[IRanges]{findOverlaps}} generic function defined
          in the \pkg{IRanges} package.

    \item The \link{GNCList} constructor and class for preprocessing and
          representing a \link{GenomicRanges} or object as a data structure
          based on Nested Containment Lists.

    \item The \link{GRanges} and \link{GRangesList} classes.
  }
}

\examples{
## ---------------------------------------------------------------------
## BASIC EXAMPLES
## ---------------------------------------------------------------------

## GRanges object:
gr <- GRanges(
        seqnames=Rle(c("chr1", "chr2", "chr1", "chr3"), c(1, 3, 2, 4)),
        ranges=IRanges(1:10, width=10:1, names=head(letters,10)),
        strand=Rle(strand(c("-", "+", "*", "+", "-")), c(1, 2, 2, 3, 2)),
        score=1:10,
        GC=seq(1, 0, length=10)
      )
gr

## GRangesList object:
gr1 <- GRanges(seqnames="chr2", ranges=IRanges(4:3, 6),
               strand="+", score=5:4, GC=0.45)
gr2 <- GRanges(seqnames=c("chr1", "chr1"),
               ranges=IRanges(c(7,13), width=3),
               strand=c("+", "-"), score=3:4, GC=c(0.3, 0.5))
gr3 <- GRanges(seqnames=c("chr1", "chr2"),
               ranges=IRanges(c(1, 4), c(3, 9)),
               strand=c("-", "-"), score=c(6L, 2L), GC=c(0.4, 0.1))
grl <- GRangesList("gr1"=gr1, "gr2"=gr2, "gr3"=gr3)

## Overlapping two GRanges objects:
table(!is.na(findOverlaps(gr, gr1, select="arbitrary")))
countOverlaps(gr, gr1)
findOverlaps(gr, gr1)
subsetByOverlaps(gr, gr1)

countOverlaps(gr, gr1, type="start")
findOverlaps(gr, gr1, type="start")
subsetByOverlaps(gr, gr1, type="start")

findOverlaps(gr, gr1, select="first")
findOverlaps(gr, gr1, select="last")

findOverlaps(gr1, gr)
findOverlaps(gr1, gr, type="start")
findOverlaps(gr1, gr, type="within")
findOverlaps(gr1, gr, type="equal")

## ---------------------------------------------------------------------
## MORE EXAMPLES
## ---------------------------------------------------------------------

table(!is.na(findOverlaps(gr, gr1, select="arbitrary")))
countOverlaps(gr, gr1)
findOverlaps(gr, gr1)
subsetByOverlaps(gr, gr1)

## Overlaps between a GRanges and a GRangesList object:

table(!is.na(findOverlaps(grl, gr, select="first")))
countOverlaps(grl, gr)
findOverlaps(grl, gr)
subsetByOverlaps(grl, gr)
countOverlaps(grl, gr, type="start")
findOverlaps(grl, gr, type="start")
subsetByOverlaps(grl, gr, type="start")
findOverlaps(grl, gr, select="first")

table(!is.na(findOverlaps(grl, gr1, select="first")))
countOverlaps(grl, gr1)
findOverlaps(grl, gr1)
subsetByOverlaps(grl, gr1)
countOverlaps(grl, gr1, type="start")
findOverlaps(grl, gr1, type="start")
subsetByOverlaps(grl, gr1, type="start")
findOverlaps(grl, gr1, select="first")

## Overlaps between two GRangesList objects:
countOverlaps(grl, rev(grl))
findOverlaps(grl, rev(grl))
subsetByOverlaps(grl, rev(grl))
}

\keyword{methods}
\keyword{utilities}
