\name{absoluteRanges}

\alias{absoluteRanges}
\alias{relativeRanges}

\title{Transform genomic ranges into "absolute" ranges}

\description{
  \code{absoluteRanges} transforms the genomic ranges in \code{x} into
  \emph{absolute} ranges i.e. into ranges counted from the beginning of
  the virtual sequence obtained by concatenating all the sequences in the
  underlying genome (in the order reported by \code{seqlevels(x)}).

  \code{relativeRanges} performs the reverse transformation.
}

\usage{
absoluteRanges(x)
relativeRanges(x, seqlengths)
}

\arguments{
  \item{x}{
    A \link{GenomicRanges} object with sequence lengths available for
    \code{absoluteRanges}.

    A \link{Ranges} object for \code{relativeRanges}.
  }
  \item{seqlengths}{
    Either a named numeric vector of chromosome lengths or a \link{Seqinfo}
    object. More precisely, if a named numeric vector, it must have a length
    >= 1, cannot contain NAs or negative values, and cannot have duplicated
    names. If a \link{Seqinfo} object, then it's first replaced with the
    vector of sequence lengths stored in the object (extracted from the object
    with the \code{\link{seqlengths}} getter), then the restrictions described
    previously apply to this vector.
  }
}

\value{
  An \link[IRanges]{IRanges} object for \code{absoluteRanges}.

  A \link{GRanges} object for \code{relativeRanges}.
}

\author{
  H. Pages
}

\seealso{
  \itemize{
    \item \link{GRanges} objects.

    \item \link[IRanges]{IRanges} objects in the \pkg{IRanges} package.

    \item \link{Seqinfo} objects and the \code{\link{seqlengths}} getter.

    \item The \code{\link{tileGenome}} function for putting tiles on a
          genome.
  }
}

\examples{
gr0 <- GRanges(Rle(c("chr2", "chr1", "chr3", "chr1"), 4:1),
               IRanges(1:10, width=5),
               seqinfo=Seqinfo(c("chr1", "chr2", "chr3"), c(100, 50, 20)))

ir0 <- absoluteRanges(gr0)
ir0

gr1 <- relativeRanges(ir0, seqlengths(gr0))
gr1

## Sanity check:
stopifnot(all(gr0 == gr1))
}

\keyword{manip}
