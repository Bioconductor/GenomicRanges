\name{GenomicRanges-comparison}

\alias{GenomicRanges-comparison}

\alias{compare,GenomicRanges,GenomicRanges-method}

\alias{==,GenomicRanges,GenomicRanges-method}
\alias{<=,GenomicRanges,GenomicRanges-method}

\alias{duplicated,GenomicRanges-method}
\alias{duplicated.GenomicRanges}

\alias{match,GenomicRanges,GenomicRanges-method}
\alias{\%in\%,GenomicRanges,GenomicRanges-method}

\alias{findMatches,GenomicRanges,GenomicRanges-method}
\alias{countMatches,GenomicRanges,GenomicRanges-method}

\alias{order,GenomicRanges-method}
\alias{rank,GenomicRanges-method}


\title{Comparing and ordering genomic ranges}

\description{
  Methods for comparing and ordering the elements in one or more
  \link{GenomicRanges} objects.
}

\usage{
## Element-wise (aka "parallel") comparison of 2 GenomicRanges objects
## -------------------------------------------------------------------

\S4method{==}{GenomicRanges,GenomicRanges}(e1, e2)

\S4method{<=}{GenomicRanges,GenomicRanges}(e1, e2)

## duplicated()
## ------------

\S4method{duplicated}{GenomicRanges}(x, incomparables=FALSE, fromLast=FALSE,
           method=c("auto", "quick", "hash"))

## match()
## -------

\S4method{match}{GenomicRanges,GenomicRanges}(x, table, nomatch=NA_integer_, incomparables=NULL,
      method=c("auto", "quick", "hash"), ignore.strand=FALSE, match.if.overlap=FALSE)

## order() and related methods
## ----------------------------

\S4method{order}{GenomicRanges}(..., na.last=TRUE, decreasing=FALSE)

\S4method{rank}{GenomicRanges}(x, na.last=TRUE,
     ties.method=c("average", "first", "random", "max", "min"))

## Generalized element-wise (aka "parallel") comparison of 2 GenomicRanges
## objects
## ------------------------------------------------------------------------

\S4method{compare}{GenomicRanges,GenomicRanges}(x, y)
}

\arguments{
  \item{e1, e2, x, table, y}{
    \link{GenomicRanges} objects.
  }
  \item{incomparables}{
    Not supported.
  }
  \item{fromLast, method, nomatch}{
    See \code{?`\link[IRanges]{Ranges-comparison}`} in the IRanges
    package for a description of these arguments.
  }
  \item{ignore.strand}{
    Whether or not the strand should be ignored when comparing 2 genomic
    ranges.
  }
  \item{match.if.overlap}{
    This argument is deprecated in BioC 2.13 and won't be supported
    anymore in BioC 2.14. Please use
    \code{findOverlaps(x, table, select="first", ignore.strand=ignore.strand)}
    instead of
    \code{match(x, table, ignore.strand=ignore.strand, match.if.overlap=TRUE)}.
  }
  \item{...}{
    Additional \link{GenomicRanges} objects used for breaking ties.
  }
  \item{na.last}{
    Ignored.
  }
  \item{decreasing}{
    \code{TRUE} or \code{FALSE}.
  }
  \item{ties.method}{
    A character string specifying how ties are treated. Only \code{"first"}
    is supported for now.
  }
}

\details{
  Two elements of a \link{GenomicRanges} object (i.e. two genomic ranges) are
  considered equal iff they are on the same underlying sequence and strand,
  and have the same start and width. \code{duplicated()} and \code{unique()}
  on a \link{GenomicRanges} object are conforming to this.

  The "natural order" for the elements of a \link{GenomicRanges} object is to
  order them (a) first by sequence level, (b) then by strand, (c) then by
  start, (d) and finally by width.
  This way, the space of genomic ranges is totally ordered.
  Note that the \code{reduce} method for \link{GenomicRanges} uses this
  "natural order" implicitly. Also, note that, because we already do (c)
  and (d) for regular ranges (see \code{?`\link[IRanges]{Ranges-comparison}`}),
  genomic ranges that belong to the same underlying sequence and strand are
  ordered like regular ranges.
  \code{order()}, \code{sort()}, and \code{rank()} on a \link{GenomicRanges}
  object are using this "natural order".

  Also the \code{==}, \code{!=}, \code{<=}, \code{>=}, \code{<} and \code{>}
  operators between 2 \link{GenomicRanges} objects are using this "natural
  order".
}

\author{H. Pages}

\seealso{
  \itemize{
    \item The \link{GenomicRanges} class.
    \item \link[IRanges]{Ranges-comparison} in the IRanges
          package for comparing and ordering genomic ranges.
    \item \link[GenomicRanges]{intra-range-methods} and
          \link[GenomicRanges]{inter-range-methods} for
          intra and inter range transformations.
    \item \link[GenomicRanges]{setops-methods} for set operations on
          \link{GenomicRanges} objects.
    \item \link[GenomicRanges]{findOverlaps-methods} for finding
          overlapping genomic ranges.
  }
}

\examples{
gr0 <- GRanges(
        seqnames=Rle(c("chr1", "chr2", "chr1", "chr3"), c(1, 3, 2, 4)),
        ranges=IRanges(c(1:9,7L), end=10),
        strand=Rle(strand(c("-", "+", "*", "+", "-")), c(1, 2, 2, 3, 2)),
        seqlengths=c(chr1=11, chr2=12, chr3=13))
gr <- c(gr0, gr0[7:3])

## ---------------------------------------------------------------------
## A. ELEMENT-WISE (AKA "PARALLEL") COMPARISON OF 2 GenomicRanges OBJECTS
## ---------------------------------------------------------------------
gr[2] == gr[2]  # TRUE
gr[2] == gr[5]  # FALSE
gr == gr[4]
gr >= gr[3]

## ---------------------------------------------------------------------
## B. duplicated(), unique()
## ---------------------------------------------------------------------
duplicated(gr)
unique(gr)

## ---------------------------------------------------------------------
## C. match(), %in%
## ---------------------------------------------------------------------
table <- gr[1:7]
match(gr, table)
match(gr, table, ignore.strand=TRUE)

gr \%in\% table  # Warning! The warning will be removed in BioC 2.14.
## In the meantime, use suppressWarnings() to suppress the warning:
suppressWarnings(gr \%in\% table)

## ---------------------------------------------------------------------
## D. findMatches(), countMatches()
## ---------------------------------------------------------------------
findMatches(gr, table)
countMatches(gr, table)

findMatches(gr, table, ignore.strand=TRUE)
countMatches(gr, table, ignore.strand=TRUE)

gr_levels <- unique(gr)
countMatches(gr_levels, gr)

## ---------------------------------------------------------------------
## E. order() AND RELATED METHODS
## ---------------------------------------------------------------------
order(gr)
sort(gr)
rank(gr)

## ---------------------------------------------------------------------
## F. GENERALIZED ELEMENT-WISE COMPARISON OF 2 GenomicRanges OBJECTS
## ---------------------------------------------------------------------
gr2 <- GRanges(c(rep("chr1", 12), "chr2"), IRanges(c(1:11, 6:7), width=3))
strand(gr2)[12] <- "+"
gr3 <- GRanges("chr1", IRanges(5, 9))

compare(gr2, gr3)
rangeComparisonCodeToLetter(compare(gr2, gr3))
}

\keyword{methods}
