\name{intra-range-methods}

\alias{intra-range-methods}

\alias{shift,GenomicRanges-method}
\alias{narrow,GenomicRanges-method}
\alias{flank,GenomicRanges-method}
\alias{promoters,GRanges-method}
\alias{promoters,GRangesList-method}
\alias{resize,GenomicRanges-method}
\alias{Ops,GenomicRanges,numeric-method}
\alias{restrict,GenomicRanges-method}
\alias{trim,GenomicRanges-method}


\title{Intra range transformations of a GenomicRanges object}

\description{
  See \code{?`\link[IRanges]{intra-range-methods}`} and
  \code{?`\link[IRanges]{inter-range-methods}`} in the IRanges package
  for a quick introduction to intra range and inter range transformations.

  This man page documents intra range transformations of a
  \link{GenomicRanges} object (i.e. of an object that belongs to the
  \link{GenomicRanges} class or one of its subclasses, this includes
  for example \link{GRanges} objects). 

  See \code{?`\link[GenomicRanges]{inter-range-methods}`} for
  inter range transformations of a GenomicRanges object.
}

\usage{
\S4method{shift}{GenomicRanges}(x, shift=0L, use.names=TRUE)

\S4method{narrow}{GenomicRanges}(x, start=NA, end=NA, width=NA, use.names=TRUE)

\S4method{flank}{GenomicRanges}(x, width, start=TRUE, both=FALSE,
        use.names=TRUE, ignore.strand=FALSE)

\S4method{promoters}{GRanges}(x, upstream=2000, downstream=200, ...)
\S4method{promoters}{GRangesList}(x, upstream=2000, downstream=200, ...)

\S4method{resize}{GenomicRanges}(x, width, fix="start", use.names=TRUE,
        ignore.strand=FALSE)

\S4method{restrict}{GenomicRanges}(x, start=NA, end=NA, keep.all.ranges=FALSE,
        use.names=TRUE)

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
      \link[IRanges]{Ranges} objects. See 
      \code{?`\link[IRanges]{intra-range-methods}`} for the details.
    }
    \item(){\code{narrow} behaves like the \code{narrow} method for 
      \link[IRanges]{Ranges} objects. See 
      \code{?`\link[IRanges]{intra-range-methods}`} for the details.
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
      \code{promoters} generates promoter ranges for each range in 
      \code{x} relative to the transcription start site (TSS), where TSS is
      \code{start(x)}. The promoter range is expanded around the TSS 
      according to the \code{upsteam} and \code{downstream} arguments.
      \code{upstream} represents the number of nucleotides in the 5' 
      direction and \code{downstream} the number in the 3' direction. 
      The full range is defined as,

      (start(x) - upstream) to (start(x) + downstream - 1).

      Ranges on the \code{*} strand are treated the same as those on the
      \code{+} strand. When no seqlengths are present in \code{x}, it is 
      possible to have non-positive start values in the promoter ranges. This 
      occurs when (TSS - \code{upstream}) < 1. In the equal but opposite case, 
      the \code{end} values of the ranges may extend beyond the chromosome end 
      when (TSS + \code{downstream} + 1) > 'chromosome end'. When 
      \code{seqlengths} are not \code{NA} the promoter ranges are kept within 
      the bounds of the defined \code{seqlengths}.
    }
    \item{}{\code{resize} returns an object of the same type and length as 
      \code{x} containing intervals that have been resized to width 
      \code{width} based on the \code{strand(x)} values.  Elements where 
      \code{strand(x) == "+"} or \code{strand(x) == "*"} are anchored at 
      \code{start(x)} and elements where \code{strand(x) == "-"} are anchored 
      at the \code{end(x)}. The \code{use.names} argument determines whether 
      or not to keep the names on the ranges.
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
    \item{}{\code{trim} subsets the ranges in \code{x} to fall within the valid
      seqlengths. If no seqlenths are present, negative start values are set 
      to 1L and end values are not touched.
    }
  }
}

\author{P. Aboyoun and V. Obenchain <vobencha@fhcrc.org>}

\seealso{
  \itemize{
    \item The \link{GenomicRanges} and \link{GRanges} classes.
    \item The \link[IRanges]{Ranges} class in the IRanges package.
    \item The \link[IRanges]{intra-range-methods} man page in the
          IRanges package.
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

shift(gr, 1)
narrow(gr[-10], start=2, end=-2)
flank(gr, 10)
resize(gr, 10)
restrict(gr, start=3, end=7)

gr <- GRanges("chr1", IRanges(rep(10, 3), width=6), c("+", "-", "*"))
promoters(gr, 2, 2)
}
 
\keyword{utilities}
