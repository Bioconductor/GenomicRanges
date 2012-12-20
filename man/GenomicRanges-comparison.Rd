\name{GenomicRanges-comparison}

\alias{GenomicRanges-comparison}

\alias{compare,GenomicRanges,GenomicRanges-method}

\alias{==,GenomicRanges,GenomicRanges-method}
\alias{<=,GenomicRanges,GenomicRanges-method}

\alias{duplicated,GenomicRanges-method}

\alias{match,GenomicRanges,GenomicRanges-method}
\alias{\%in\%,GenomicRanges,GenomicRanges-method}

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
      method=c("auto", "quick", "hash"), match.if.overlap=FALSE)

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
  \item{match.if.overlap}{
    You temporarily need to explicitly provide this argument otherwise
    \code{match()} will issue a warning.

    Starting with BioC 2.12, the default behavior of \code{match()} on
    \link{GenomicRanges} objects has changed to use \emph{equality} instead
    of \emph{overlap} for comparing elements between \link{GenomicRanges}
    objects \code{x} and \code{table}. Now \code{x[i]} and \code{table[j]} are
    considered to match when they are equal (i.e. \code{x[i] == table[j]}),
    instead of when they overlap. This new behavior is consistent with
    \code{base::\link[base]{match}()}.

    If you need the new behavior, use \code{match.if.overlap=FALSE}.

    If you need the old behavior, you can either do:
    \code{findOverlaps(x, table, select="first")}
    which is the recommended way. Alternatively, you can call \code{match()}
    with \code{match.if.overlap=TRUE}.
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
  }
}

\examples{
gr <- GRanges(
        seqnames=Rle(c("chr1", "chr2", "chr1", "chr3"), c(1, 3, 2, 4)),
        ranges=IRanges(1:10, end=10),
        strand=Rle(strand(c("-", "+", "*", "+", "-")), c(1, 2, 2, 3, 2)),
        seqlengths=c(chr1=11, chr2=12, chr3=13))

duplicated(gr)
duplicated(c(gr[4], gr))
unique(gr)
unique(c(gr[4], gr))
order(gr)
sort(gr)
rank(gr)

gr[2] == gr[2]  # TRUE
gr[2] == gr[5]  # FALSE
gr == gr[4]
gr >= gr[3]
compare(gr, gr[3])
rangeComparisonCodeToLetter(compare(gr, gr[3]))
}

\keyword{methods}
