\name{GenomicRanges-comparison}

\alias{GenomicRanges-comparison}

\alias{==,GenomicRanges,GenomicRanges-method}
\alias{duplicated,GenomicRanges-method}
\alias{unique,GenomicRanges-method}

\alias{<=,GenomicRanges,GenomicRanges-method}
\alias{>=,GenomicRanges,GenomicRanges-method}
\alias{<,GenomicRanges,GenomicRanges-method}
\alias{>,GenomicRanges,GenomicRanges-method}
\alias{order,GenomicRanges-method}
\alias{sort,GenomicRanges-method}
\alias{rank,GenomicRanges-method}

\alias{compare,GenomicRanges,GenomicRanges-method}


\title{Comparing and ordering genomic ranges}

\description{
  Methods for comparing and ordering the elements in one or more
  \link{GenomicRanges} objects.
}

\details{
  Two elements of a \link{GenomicRanges} object (i.e. two genomic ranges) are
  considered equal iff they are on the same underlying sequence and strand,
  and have the same start and width.
  The \code{duplicated} and \code{unique} methods for GenomicRanges objects
  are using this equality.

  The "natural order" for the elements of a \link{GenomicRanges} object is to
  order them (a) first by sequence level, (b) then by strand, (c) then by
  start, (d) and finally by width.
  This way, the space of genomic ranges is totally ordered.
  Note that the \code{reduce} method for \link{GenomicRanges} uses this
  "natural order" implicitly. Also, note that, because we already do (c)
  and (d) for regular ranges (see \code{?`\link[IRanges]{Ranges-comparison}`}),
  genomic ranges that belong to the same underlying sequence and strand are
  ordered like regular ranges.
  The \code{order}, \code{sort} and \code{rank} methods for \link{GenomicRanges}
  objects are using this "natural order".

  Also the \code{==}, \code{!=}, \code{<=}, \code{>=}, \code{<} and \code{>}
  operators between 2 \link{GenomicRanges} objects are using this "natural
  order".
}

\seealso{
  \link{GenomicRanges-class},
  \link[IRanges]{Ranges-comparison}
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
