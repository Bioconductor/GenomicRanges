\name{intra-range-methods}

\alias{intra-range-methods}

\alias{update_ranges,GenomicRanges-method}

\alias{shift,GenomicRanges-method}

\alias{narrow,GenomicRanges-method}

\alias{resize}
\alias{resize,GenomicRanges-method}

\alias{flank}
\alias{flank,GenomicRanges-method}

\alias{promoters}
\alias{promoters,GenomicRanges-method}

\alias{restrict}
\alias{restrict,GenomicRanges-method}

\alias{trim}
\alias{trim,GenomicRanges-method}
\alias{trim,GRangesList-method}

\title{Intra range transformations of a GRanges or GRangesList object}

\description{
  This man page documents \emph{intra range transformations} of a
  \link{GenomicRanges} object (i.e. of an object that belongs to the
  \link{GenomicRanges} class or one of its subclasses, this includes
  for example \link{GRanges} objects), or a \link{GRangesList} object.

  See \code{?`\link[IRanges]{intra-range-methods}`} and
  \code{?`\link[IRanges]{inter-range-methods}`} in the \pkg{IRanges}
  package for a quick introduction to \emph{intra range} and \emph{inter
  range transformations}.

  \emph{Intra range} methods for \link[GenomicAlignments]{GAlignments}
  and \link[GenomicAlignments]{GAlignmentsList} objects are defined and
  documented in the \pkg{GenomicAlignments} package.

  See \code{?`\link[GenomicRanges]{inter-range-methods}`} for
  \emph{inter range transformations} of a \link{GenomicRanges} or
  \link{GRangesList} object.
}

\usage{
\S4method{shift}{GenomicRanges}(x, shift=0L, use.names=TRUE)

\S4method{narrow}{GenomicRanges}(x, start=NA, end=NA, width=NA, use.names=TRUE)

\S4method{resize}{GenomicRanges}(x, width, fix="start", use.names=TRUE, ignore.strand=FALSE)

\S4method{flank}{GenomicRanges}(x, width, start=TRUE, both=FALSE, use.names=TRUE,
      ignore.strand=FALSE)

\S4method{promoters}{GenomicRanges}(x, upstream=2000, downstream=200, use.names=TRUE)

\S4method{restrict}{GenomicRanges}(x, start=NA, end=NA, keep.all.ranges=FALSE, use.names=TRUE)

\S4method{trim}{GenomicRanges}(x, use.names=TRUE)
}

\arguments{
  \item{x}{A \link{GenomicRanges} object.}
  \item{shift, use.names, start, end, width, both, fix, keep.all.ranges,
        upstream, downstream}{
    See \code{?`\link[IRanges]{intra-range-methods}`}.
  }
  \item{ignore.strand}{
    \code{TRUE} or \code{FALSE}. Whether the strand of the input ranges
    should be ignored or not. See details below.
  }
}

\details{
  \itemize{
    \item{}{\code{shift} behaves like the \code{shift} method for
      \link[IRanges]{IntegerRanges} objects. See
      \code{?`\link[IRanges]{intra-range-methods}`} for the details.
    }
    \item{}{\code{narrow} on a \link{GenomicRanges} object behaves
      like on an \link[IRanges]{IntegerRanges} object. See
      \code{?`\link[IRanges]{intra-range-methods}`} for the details.

      A major difference though is that it returns a \link{GenomicRanges}
      object instead of an \link[IRanges]{IntegerRanges} object.
      The returned object is \emph{parallel} (i.e. same length and names)
      to the original object \code{x}.
    }
    \item{}{\code{resize} returns an object of the same type and length as
      \code{x} containing intervals that have been resized to width
      \code{width} based on the \code{strand(x)} values.  Elements where
      \code{strand(x) == "+"} or \code{strand(x) == "*"} are anchored at
      \code{start(x)} and elements where \code{strand(x) == "-"} are anchored
      at the \code{end(x)}. The \code{use.names} argument determines whether
      or not to keep the names on the ranges.
    }
    \item{}{\code{flank} returns an object of the same type and length
      as \code{x} containing intervals of width \code{width} that flank
      the intervals in \code{x}. The \code{start} argument takes a
      logical indicating whether \code{x} should be flanked at the
      "start" (\code{TRUE}) or the "end" (\code{FALSE}), which for
      \code{strand(x) != "-"} is \code{start(x)} and \code{end(x)}
      respectively and for \code{strand(x) == "-"} is \code{end(x)} and
      \code{start(x)} respectively. The \code{both} argument takes a
      single logical value indicating whether the flanking region
      \code{width} positions extends \emph{into} the range. If
      \code{both=TRUE}, the resulting range thus straddles the end
      point, with \code{width} positions on either side.
    }
    \item{}{
      \code{promoters} returns an object of the same type and length
      as \code{x} containing promoter ranges. Promoter ranges extend
      around the transcription start site (TSS) which is defined as
      \code{start(x)} for ranges on the \code{+} or \code{*} strand
      and as \code{end(x)} for ranges on the \code{-} strand.
      The \code{upstream} and \code{downstream} arguments define the
      number of nucleotides in the 5' and 3' direction, respectively.
      More precisely, the output range is defined as
      \preformatted{    (start(x) - upstream) to (start(x) + downstream - 1)
      }
      for ranges on the \code{+} or \code{*} strand, and as
      \preformatted{    (end(x) - downstream + 1) to (end(x) + upstream)
      }
      for ranges on the \code{-} strand.

      Note that the returned object might contain \emph{out-of-bound} ranges
      i.e. ranges that start before the first nucleotide position and/or end
      after the last nucleotide position of the underlying sequence.
    }
    \item{}{\code{restrict} returns an object of the same type and length as
      \code{x} containing restricted ranges for distinct seqnames. The
      \code{start} and \code{end} arguments can be a named numeric vector of
      seqnames for the ranges to be resticted or a numeric vector or length 1
      if the restriction operation is to be applied to all the sequences in
      \code{x}. See \code{?`\link[IRanges]{intra-range-methods}`} for more
      information about range restriction and for a description of the optional
      arguments.
    }
    \item{}{\code{trim} trims out-of-bound ranges located on non-circular
      sequences whose length is not NA.
    }
  }
}

\author{P. Aboyoun and V. Obenchain}

\seealso{
  \itemize{
    \item \link{GenomicRanges}, \link{GRanges}, and \link{GRangesList} objects.

    \item The \link[IRanges]{intra-range-methods} man page in the
          \pkg{IRanges} package.

    \item The \link[IRanges]{IntegerRanges} class in the \pkg{IRanges} package.
  }
}

\examples{
## ---------------------------------------------------------------------
## A. ON A GRanges OBJECT
## ---------------------------------------------------------------------
gr <- GRanges(
        seqnames=Rle(paste("chr", c(1, 2, 1, 3), sep=""), c(1, 3, 2, 4)),
        ranges=IRanges(1:10, width=10:1, names=letters[1:10]),
        strand=Rle(strand(c("-", "+", "*", "+", "-")), c(1, 2, 2, 3, 2)),
        score=1:10,
        GC=seq(1, 0, length=10)
      )
gr

shift(gr, 1)
narrow(gr[-10], start=2, end=-2)
resize(gr, width=10)
flank(gr, width=10)
restrict(gr, start=3, end=7)

gr <- GRanges("chr1", IRanges(rep(10, 3), width=6), c("+", "-", "*"))
promoters(gr, 2, 2)

## ---------------------------------------------------------------------
## B. ON A GRangesList OBJECT
## ---------------------------------------------------------------------
gr1 <- GRanges("chr2", IRanges(3, 6))
gr2 <- GRanges(c("chr1", "chr1"), IRanges(c(7,13), width=3),
               strand=c("+", "-"))
gr3 <- GRanges(c("chr1", "chr2"), IRanges(c(1, 4), c(3, 9)),
               strand="-")
grl <- GRangesList(gr1= gr1, gr2=gr2, gr3=gr3)
grl

resize(grl, width=20)
flank(grl, width=20)
restrict(grl, start=3)
}

\keyword{methods}
\keyword{utilities}
