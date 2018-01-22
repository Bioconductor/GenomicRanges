\name{GNCList-class}
\docType{class}

% GNCList objects:
\alias{class:GNCList}
\alias{GNCList-class}
\alias{GNCList}

\alias{granges,GNCList-method}
\alias{length,GNCList-method}
\alias{names,GNCList-method}
\alias{seqnames,GNCList-method}
\alias{start,GNCList-method}
\alias{end,GNCList-method}
\alias{width,GNCList-method}
\alias{ranges,GNCList-method}
\alias{strand,GNCList-method}
\alias{seqinfo,GNCList-method}
\alias{coerce,GNCList,GRanges-method}
\alias{coerce,GenomicRanges,GNCList-method}


\title{GNCList objects}

\description{
  The GNCList class is a container for storing the Nested Containment List
  representation of a vector of genomic ranges (typically represented as
  a \link{GRanges} object).
  To preprocess a \link{GRanges} object, simply call the \code{GNCList}
  constructor function on it. The resulting GNCList object can then be used
  for efficient overlap-based operations on the genomic ranges.
}

\usage{
GNCList(x)
}

\arguments{
  \item{x}{
    The \link{GRanges} (or more generally \link{GenomicRanges}) object to
    preprocess.
  }
}

\details{
  The \pkg{IRanges} package also defines the \code{\link[IRanges]{NCList}}
  and \code{\link[IRanges]{NCLists}} constructors and classes for
  preprocessing and representing a \link[IRanges]{IntegerRanges} or
  \link[IRanges]{IntegerRangesList} object as a data structure based on
  Nested Containment Lists.

  Note that GNCList objects (introduced in BioC 3.1) are replacements for
  GIntervalTree objects (BioC < 3.1).

  See \code{?\link[IRanges]{NCList}} in the \pkg{IRanges} package for
  some important differences between the new algorithm based on Nested
  Containment Lists and the old algorithm based on interval trees.
  In particular, the new algorithm supports preprocessing of a
  \link{GenomicRanges} object with ranges defined on circular sequences
  (e.g. on the mitochnodrial chromosome). See below for some examples.
}

\value{
  A GNCList object.
}

\author{H. Pagès}

\references{
  Alexander V. Alekseyenko and Christopher J. Lee --
  Nested Containment List (NCList): a new algorithm for accelerating interval
  query of genome alignment and interval databases.
  Bioinformatics (2007) 23 (11): 1386-1393.
  doi: 10.1093/bioinformatics/btl647
}

\seealso{
  \itemize{
    \item The \code{\link[IRanges]{NCList}} and \code{\link[IRanges]{NCLists}}
          constructors and classs defined in the \pkg{IRanges} package.

    \item \code{\link{findOverlaps}} for finding/counting interval overlaps
          between two \emph{range-based} objects.

    \item \link{GRanges} objects.
  }
}

\examples{
## The examples below are for illustration purpose only and do NOT
## reflect typical usage. This is because, for a one time use, it is
## NOT advised to explicitely preprocess the input for findOverlaps()
## or countOverlaps(). These functions will take care of it and do a
## better job at it (by preprocessing only what's needed when it's
## needed, and release memory as they go).

## ---------------------------------------------------------------------
## PREPROCESS QUERY OR SUBJECT
## ---------------------------------------------------------------------

query <- GRanges(Rle(c("chrM", "chr1", "chrM", "chr1"), 4:1),
                 IRanges(1:10, width=5), strand=rep(c("+", "-"), 5))
subject <- GRanges(Rle(c("chr1", "chr2", "chrM"), 3:1),
                   IRanges(6:1, width=5), strand="+")

## Either the query or the subject of findOverlaps() can be preprocessed:

ppsubject <- GNCList(subject)
hits1a <- findOverlaps(query, ppsubject)
hits1a
hits1b <- findOverlaps(query, ppsubject, ignore.strand=TRUE)
hits1b

ppquery <- GNCList(query)
hits2a <- findOverlaps(ppquery, subject)
hits2a
hits2b <- findOverlaps(ppquery, subject, ignore.strand=TRUE)
hits2b

## Note that 'hits1a' and 'hits2a' contain the same hits but not
## necessarily in the same order.
stopifnot(identical(sort(hits1a), sort(hits2a)))
## Same for 'hits1b' and 'hits2b'.
stopifnot(identical(sort(hits1b), sort(hits2b)))

## ---------------------------------------------------------------------
## WITH CIRCULAR SEQUENCES
## ---------------------------------------------------------------------

seqinfo <- Seqinfo(c("chr1", "chr2", "chrM"),
                   seqlengths=c(100, 50, 10),
                   isCircular=c(FALSE, FALSE, TRUE))
seqinfo(query) <- seqinfo[seqlevels(query)]
seqinfo(subject) <- seqinfo[seqlevels(subject)]

ppsubject <- GNCList(subject)
hits3 <- findOverlaps(query, ppsubject)
hits3

## Circularity introduces more hits:

stopifnot(all(hits1a \%in\% hits3))
new_hits <- setdiff(hits3, hits1a)
new_hits  # 1 new hit
query[queryHits(new_hits)]
subject[subjectHits(new_hits)]  # positions 11:13 on chrM are the same
                                # as positions 1:3

## Sanity checks:
stopifnot(identical(new_hits, Hits(9, 6, 10, 6, sort.by.query=TRUE)))
ppquery <- GNCList(query)
hits4 <- findOverlaps(ppquery, subject)
stopifnot(identical(sort(hits3), sort(hits4)))
}

\keyword{classes}
\keyword{methods}
