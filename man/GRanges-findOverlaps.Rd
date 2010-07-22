\name{GRanges-findOverlaps}

\alias{findOverlaps,GRanges,GRanges-method}
\alias{findOverlaps,GRanges,GRangesList-method}
\alias{findOverlaps,GRangesList,GRanges-method}
\alias{findOverlaps,GRangesList,GRangesList-method}
\alias{findOverlaps,RangesList,GRanges-method}
\alias{findOverlaps,GRanges,RangesList-method}
\alias{findOverlaps,GRangesList,RangesList-method}
\alias{findOverlaps,RangesList,GRangesList-method}
\alias{findOverlaps,RangedData,GRanges-method}
\alias{findOverlaps,GRanges,RangedData-method}
\alias{findOverlaps,GRangesList,RangedData-method}
\alias{findOverlaps,RangedData,GRangesList-method}

\alias{countOverlaps,GRanges,GRanges-method}
\alias{countOverlaps,GRanges,GRangesList-method}
\alias{countOverlaps,GRangesList,GRanges-method}
\alias{countOverlaps,GRangesList,GRangesList-method}
\alias{countOverlaps,RangesList,GRanges-method}
\alias{countOverlaps,GRanges,RangesList-method}
\alias{countOverlaps,GRangesList,RangesList-method}
\alias{countOverlaps,RangesList,GRangesList-method}
\alias{countOverlaps,RangedData,GRanges-method}
\alias{countOverlaps,GRanges,RangedData-method}
\alias{countOverlaps,GRangesList,RangedData-method}
\alias{countOverlaps,RangedData,GRangesList-method}

\alias{subsetByOverlaps,GRanges,GRanges-method}
\alias{subsetByOverlaps,GRanges,GRangesList-method}
\alias{subsetByOverlaps,GRangesList,GRanges-method}
\alias{subsetByOverlaps,GRangesList,GRangesList-method}
\alias{subsetByOverlaps,RangesList,GRanges-method}
\alias{subsetByOverlaps,GRanges,RangesList-method}
\alias{subsetByOverlaps,GRangesList,RangesList-method}
\alias{subsetByOverlaps,RangesList,GRangesList-method}
\alias{subsetByOverlaps,RangedData,GRanges-method}
\alias{subsetByOverlaps,GRanges,RangedData-method}
\alias{subsetByOverlaps,GRangesList,RangedData-method}
\alias{subsetByOverlaps,RangedData,GRangesList-method}

\alias{match,GRanges,GRanges-method}
\alias{match,GRanges,GRangesList-method}
\alias{match,GRangesList,GRanges-method}
\alias{match,GRangesList,GRangesList-method}
\alias{match,RangesList,GRanges-method}
\alias{match,GRanges,RangesList-method}
\alias{match,GRangesList,RangesList-method}
\alias{match,RangesList,GRangesList-method}
\alias{match,RangedData,GRanges-method}
\alias{match,GRanges,RangedData-method}
\alias{match,GRangesList,RangedData-method}
\alias{match,RangedData,GRangesList-method}

\alias{\%in\%,GRanges,GRanges-method}
\alias{\%in\%,GRanges,GRangesList-method}
\alias{\%in\%,GRangesList,GRanges-method}
\alias{\%in\%,GRangesList,GRangesList-method}
\alias{\%in\%,RangesList,GRanges-method}
\alias{\%in\%,GRanges,RangesList-method}
\alias{\%in\%,GRangesList,RangesList-method}
\alias{\%in\%,RangesList,GRangesList-method}
\alias{\%in\%,RangedData,GRanges-method}
\alias{\%in\%,GRanges,RangedData-method}
\alias{\%in\%,GRangesList,RangedData-method}
\alias{\%in\%,RangedData,GRangesList-method}


\title{GRanges and GRangesList Interval Overlaps}

\description{
  Finds interval overlaps between a GRanges/GRangesList object
  and a GRanges/GRangesList object.
}

\usage{
  \S4method{findOverlaps}{GRanges,GRanges}(query, subject, maxgap = 0L, minoverlap = 1L,
             type = c("any", "start", "end"), select = c("all", "first"))
  \S4method{countOverlaps}{GRanges,GRanges}(query, subject, maxgap = 0L, minoverlap = 1L,
              type = c("any", "start", "end"))
  \S4method{subsetByOverlaps}{GRanges,GRanges}(query, subject, maxgap = 0L, minoverlap = 1L,
                 type = c("any", "start", "end"))
  \S4method{match}{GRanges,GRanges}(x, table, nomatch = NA_integer_, incomparables = NULL)
  # Also:  x \%in\% table
}

\arguments{
  \item{query, subject, x, table}{
    A \link{GRanges} or \link{GRangesList}
    object. \link[IRanges:RangesList-class]{RangesList} and
    \link[IRanges:RangedData-class]{RangedData} are also accepted
    for one of \code{query} or \code{subject} (\code{x} or \code{table}
    for \code{match}).
  }
  \item{maxgap}{
    A non-negative integer representing the maximum distance between
    a query interval and a subject interval.
  }
  \item{minoverlap}{
    Ignored.
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
  For \code{findOverlaps} either a
  \link[IRanges:RangesMatching-class]{RangesMatching} object when
  \code{select = "all"} or an integer vector when \code{select =
    "first"}. 

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

\author{P. Aboyoun and S. Falcon and M. Lawrence}

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
  subsetByOverlaps(gr, gr1)
  countOverlaps(gr, gr1, type = "start")
  findOverlaps(gr, gr1, type = "start")
  subsetByOverlaps(gr, gr1, type = "start")
  findOverlaps(gr, gr1, select = "first")

  ## Overlapping a GRanges and a GRangesList object
  table(grlist \%in\% gr)
  countOverlaps(grlist, gr)
  findOverlaps(grlist, gr)
  subsetByOverlaps(grlist, gr)
  countOverlaps(grlist, gr, type = "start")
  findOverlaps(grlist, gr, type = "start")
  subsetByOverlaps(grlist, gr, type = "start")
  findOverlaps(grlist, gr, select = "first")

  ## Overlapping two GRangesList objects
  countOverlaps(grlist, rev(grlist))
  findOverlaps(grlist, rev(grlist))
  subsetByOverlaps(grlist, rev(grlist))
}

\keyword{methods}
\keyword{utilities}
