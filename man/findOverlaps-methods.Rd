\name{findOverlaps-methods}
\alias{findOverlaps-methods}

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
\alias{findOverlaps,GappedAlignments,ANY-method}
\alias{findOverlaps,ANY,GappedAlignments-method}
\alias{findOverlaps,GappedAlignments,GappedAlignments-method}

\alias{countOverlaps,GenomicRanges,GenomicRanges-method}
\alias{countOverlaps,GRangesList,GenomicRanges-method}
\alias{countOverlaps,GenomicRanges,GRangesList-method}
\alias{countOverlaps,GRangesList,GRangesList-method}
\alias{countOverlaps,RangesList,GenomicRanges-method}
\alias{countOverlaps,RangesList,GRangesList-method}
\alias{countOverlaps,GenomicRanges,RangesList-method}
\alias{countOverlaps,GRangesList,RangesList-method}
\alias{countOverlaps,RangedData,GenomicRanges-method}
\alias{countOverlaps,RangedData,GRangesList-method}
\alias{countOverlaps,GenomicRanges,RangedData-method}
\alias{countOverlaps,GRangesList,RangedData-method}
\alias{countOverlaps,GappedAlignments,ANY-method}
\alias{countOverlaps,ANY,GappedAlignments-method}
\alias{countOverlaps,GappedAlignments,GappedAlignments-method}

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
\alias{subsetByOverlaps,GappedAlignments,ANY-method}
\alias{subsetByOverlaps,ANY,GappedAlignments-method}
\alias{subsetByOverlaps,GappedAlignments,GappedAlignments-method}

\alias{match,GenomicRanges,GenomicRanges-method}
\alias{match,GRangesList,GenomicRanges-method}
\alias{match,GenomicRanges,GRangesList-method}
\alias{match,GRangesList,GRangesList-method}
\alias{match,RangesList,GenomicRanges-method}
\alias{match,RangesList,GRangesList-method}
\alias{match,GenomicRanges,RangesList-method}
\alias{match,GRangesList,RangesList-method}
\alias{match,RangedData,GenomicRanges-method}
\alias{match,RangedData,GRangesList-method}
\alias{match,GenomicRanges,RangedData-method}
\alias{match,GRangesList,RangedData-method}
\alias{match,GappedAlignments,ANY-method}
\alias{match,ANY,GappedAlignments-method}
\alias{match,GappedAlignments,GappedAlignments-method}

\alias{\%in\%,GenomicRanges,GenomicRanges-method}
\alias{\%in\%,GRangesList,GenomicRanges-method}
\alias{\%in\%,GenomicRanges,GRangesList-method}
\alias{\%in\%,GRangesList,GRangesList-method}
\alias{\%in\%,RangesList,GenomicRanges-method}
\alias{\%in\%,RangesList,GRangesList-method}
\alias{\%in\%,GenomicRanges,RangesList-method}
\alias{\%in\%,GRangesList,RangesList-method}
\alias{\%in\%,RangedData,GenomicRanges-method}
\alias{\%in\%,RangedData,GRangesList-method}
\alias{\%in\%,GenomicRanges,RangedData-method}
\alias{\%in\%,GRangesList,RangedData-method}
\alias{\%in\%,GappedAlignments,ANY-method}
\alias{\%in\%,ANY,GappedAlignments-method}
\alias{\%in\%,GappedAlignments,GappedAlignments-method}


\title{GRanges, GRangesList and GappedAlignments Interval Overlaps}

\description{
  Finds interval overlaps between a GRanges, GRangesList or GappedAlignments
  object and another object containing ranges.
}

\usage{
\S4method{findOverlaps}{GenomicRanges,GenomicRanges}(query, subject,
    maxgap = 0L, minoverlap = 1L,
    type = c("any", "start", "end", "within", "equal"),
    select = c("all", "first"), ignore.strand = FALSE)
\S4method{countOverlaps}{GenomicRanges,GenomicRanges}(query, subject,
    maxgap = 0L, minoverlap = 1L,
    type = c("any", "start", "end", "within", "equal"), 
    ignore.strand = FALSE)
\S4method{subsetByOverlaps}{GenomicRanges,GenomicRanges}(query, subject,
    maxgap = 0L, minoverlap = 1L,
    type = c("any", "start", "end", "within", "equal"),
    ignore.strand = FALSE)
\S4method{match}{GenomicRanges,GenomicRanges}(x, table,
    nomatch = NA_integer_, incomparables = NULL)
# Also:  x \%in\% table
}

\arguments{
  \item{query, subject, x, table}{
    A \link{GRanges}, \link{GRangesList} or \link{GappedAlignments}
    object. \link[IRanges]{RangesList} and \link[IRanges]{RangedData}
    are also accepted for one of \code{query} or \code{subject}
    (\code{x} or \code{table} for \code{match}).
  }
  \item{maxgap, minoverlap, type, select}{
    See \code{\link[IRanges]{findOverlaps}} in the IRanges package for
    a description of these arguments.
  }
  \item{ignore.strand}{
    When set to \code{TRUE}, the strand information is ignored in the
    overlap calculations.
  }
  \item{nomatch}{
    The integer value to be returned in the case when no match is found.
  }
  \item{incomparables}{
    This value is ignored.
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

  When the query or the subject (or both) is a \link{GappedAlignments}
  object, it is first turned into a \link{GRangesList} object (with
  \code{as( , "GRangesList")}) and then the rules described previously
  apply.
}

\value{
  For \code{findOverlaps} either a RangesMatching object when
  \code{select = "all"} or an integer vector when \code{select = "first"}. 

  For \code{countOverlaps} an integer vector containing the tabulated
  query overlap hits.

  For \code{subsetByOverlaps} an object of the same class as \code{query}
  containing the subset that overlapped at least one entity in \code{subject}.

  For \code{match} same as \code{findOverlaps} when \code{select = "first"}.

  For \code{\%in\%} the logical vector produced by
  \code{!is.na(match(x, table))}.

  For \code{RangedData} and \code{RangesList}, with the exception of
  \code{subsetByOverlaps}, the results align to the unlisted
  form of the object. This turns out to be fairly convenient for
  \code{RangedData} (not so much for \code{RangesList}, but something
  has to give).
}

\author{P. Aboyoun, S. Falcon, M. Lawrence, N. Gopalakrishnan and H. Pages}

\seealso{
  \code{\link[IRanges]{findOverlaps}},
  \link[IRanges]{RangesMatching-class},
  \link{GRanges},
  \link{GRangesList},
  \link{GappedAlignments}
}

\examples{
## ---------------------------------------------------------------------
## WITH GRanges AND/OR GRangesList OBJECTS
## ---------------------------------------------------------------------

## GRanges object:
gr <-
  GRanges(seqnames =
          Rle(c("chr1", "chr2", "chr1", "chr3"), c(1, 3, 2, 4)),
          ranges =
          IRanges(1:10, width = 10:1, names = head(letters,10)),
          strand =
          Rle(strand(c("-", "+", "*", "+", "-")),
              c(1, 2, 2, 3, 2)),
          score = 1:10,
          GC = seq(1, 0, length=10))
gr

## GRangesList object:
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

## Overlapping two GRanges objects:
table(gr \%in\% gr1)
countOverlaps(gr, gr1)
findOverlaps(gr, gr1)
subsetByOverlaps(gr, gr1)

countOverlaps(gr, gr1, type = "start")
findOverlaps(gr, gr1, type = "start")
subsetByOverlaps(gr, gr1, type = "start")

findOverlaps(gr, gr1, select = "first")

findOverlaps(gr1, gr)
findOverlaps(gr1, gr, type = "start")
findOverlaps(gr1, gr, type = "within")
findOverlaps(gr1, gr, type = "equal")

## Overlapping a GRanges and a GRangesList object:
table(grlist \%in\% gr)
countOverlaps(grlist, gr)
findOverlaps(grlist, gr)
subsetByOverlaps(grlist, gr)
countOverlaps(grlist, gr, type = "start")
findOverlaps(grlist, gr, type = "start")
subsetByOverlaps(grlist, gr, type = "start")
findOverlaps(grlist, gr, select = "first")

## Overlapping two GRangesList objects:
countOverlaps(grlist, rev(grlist))
findOverlaps(grlist, rev(grlist))
subsetByOverlaps(grlist, rev(grlist))

## ---------------------------------------------------------------------
## WITH A GappedAlignments OBJECT
## ---------------------------------------------------------------------
library(Rsamtools)  # because file ex1.bam is in this package
galn_file <- system.file("extdata", "ex1.bam", package="Rsamtools")
galn <- readGappedAlignments(galn_file)

subject <- granges(galn)[1]

## Note the absence of query no. 9 (i.e. 'galn[9]') in this result:
as.matrix(findOverlaps(galn, subject))

## This is because, by default, findOverlaps()/countOverlaps() are
## strand specific:
galn[8:10]
countOverlaps(galn[8:10], subject)
countOverlaps(galn[8:10], subject, ignore.strand=TRUE)

## Advanced examples:
subsetByOverlaps(galn, subject)
table(match(galn, subject), useNA = "ifany")
table(galn \%in\% subject)
}

\keyword{methods}
\keyword{utilities}
