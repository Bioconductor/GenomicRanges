\name{GRanges-findOverlaps}

\alias{findOverlaps,GRanges,GRanges-method}
\alias{findOverlaps,GRanges,GRangesList-method}
\alias{findOverlaps,GRangesList,GRanges-method}
\alias{findOverlaps,GRangesList,GRangesList-method}

\alias{countOverlaps,GRanges,GRanges-method}
\alias{countOverlaps,GRanges,GRangesList-method}
\alias{countOverlaps,GRangesList,GRanges-method}
\alias{countOverlaps,GRangesList,GRangesList-method}

\alias{match,GRanges,GRanges-method}
\alias{match,GRanges,GRangesList-method}
\alias{match,GRangesList,GRanges-method}
\alias{match,GRangesList,GRangesList-method}

\alias{\%in\%,GRanges,GRanges-method}
\alias{\%in\%,GRanges,GRangesList-method}
\alias{\%in\%,GRangesList,GRanges-method}
\alias{\%in\%,GRangesList,GRangesList-method}


\title{GRanges and GRangesList Interval Overlaps}

\description{
  Finds interval overlaps between a GRanges/GRangesList object
  and a GRanges/GRangesList object.
}

\usage{
  \S4method{findOverlaps}{GRanges,GRanges}(query, subject, maxgap = 0L,
             type = c("any", "start", "end"), select = c("all", "first"))
  \S4method{countOverlaps}{GRanges,GRanges}(query, subject, maxgap = 0L,
              type = c("any", "start", "end"))
  \S4method{match}{GRanges,GRanges}(x, table, nomatch = NA_integer_, incomparables = NULL)
  # Also:  x \%in\% table
}

\arguments{
  \item{query, subject, x, table}{
    A \link{GRanges} or \link{GRangesList} object.
  }
  \item{maxgap}{
    A non-negative integer representing the maximum distance between
    a query interval and a subject interval.
  }
  \item{type}{
    The type of acceptable overlap:
    \code{"any"} - any overlap within \code{maxgap},
    \code{"start"} - the start of the \code{query} overlaps the start
    of the \code{subject} within \code{maxgap}, and
    \code{"end"} - the end of the \code{query} overlaps the end of the
    \code{subject} within \code{maxgap}.
  }
  \item{select}{
    Overlaps to return:
    \code{"all"} - select all overlaps, and
    \code{"first"} - select the first overlap.
  }
  \item{nomatch}{
    The integer value to be returned in the case when no match is found.
  }
  \item{incomparables}{
    This value is ignored.
  }
}

\details{
  The \code{findOverlaps} methods involving \link{GRanges} and
  \link{GRangesList} objects use the triplet (sequence name, range,
  strand) to determine which features (see paragraph below for the
  definition of feature) from the \code{query} overlap which features
  in the \code{subject}, where a strand value of \code{"*"} is treated
  as occurring on both the \code{"+"} and \code{"-"} strand. An overlap
  is recorded when a feature in the \code{query} and a feature in the
  \code{subject} have the same sequence name, have a compatible pairing
  of strands (e.g. \code{"+"}/\code{"+"}, \code{"-"}/\code{"-"},
  \code{"*"}/\code{"+"}, \code{"*"}/\code{"-"}, etc.), and satisfy the
  interval overlap requirements.

  In the context of \code{findOverlaps}, a feature is a collection of
  ranges that are treated as a single entity. For \link{GRanges} objects,
  a feature is a single range; while for \link{GRangesList} objects,
  a feature and a list element containing a set of ranges. In the results,
  the features are referred to by number, which run from 1 to
  \code{length(query)}/\code{length(subject)}.
}

\value{
  For \code{findOverlaps} either a
  \link[IRanges:RangesMatching-class]{RangesMatching} object when
  \code{select = "all"} or an integer vector when \code{select = "first"}.

  For \code{countOverlaps} an integer vector containing the tabulated
  query overlap hits.

  For \code{match} same as \code{findOverlaps} when \code{select = "first"}.

  For \code{\%in\%} the logical vector produced by
  \code{!is.na(match(x, table))}.
}

\author{P. Aboyoun and S. Falcon}

\seealso{
  \code{\link[IRanges]{findOverlaps}},
  \link{GRanges},
  \link{GRangesList}
}

\examples{
  ## GRanges object
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

  ## GRangesList object
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

  ## Overlapping two GRanges objects
  table(gr \%in\% gr1)
  countOverlaps(gr, gr1)
  findOverlaps(gr, gr1)
  countOverlaps(gr, gr1, type = "start")
  findOverlaps(gr, gr1, type = "start")
  findOverlaps(gr, gr1, select = "first")

  ## Overlapping a GRanges and a GRangesList object
  table(grlist \%in\% gr)
  countOverlaps(grlist, gr)
  findOverlaps(grlist, gr)
  countOverlaps(grlist, gr, type = "start")
  findOverlaps(grlist, gr, type = "start")
  findOverlaps(grlist, gr, select = "first")

  ## Overlapping two GRangesList objects
  countOverlaps(grlist, rev(grlist))
  findOverlaps(grlist, rev(grlist))
}

\keyword{methods}
\keyword{utilities}
