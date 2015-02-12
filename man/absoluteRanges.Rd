\name{absoluteRanges}

\alias{absoluteRanges}
\alias{relativeRanges}
\alias{isSmallGenome}

\title{Transform genomic ranges into "absolute" ranges}

\description{
  \code{absoluteRanges} transforms the genomic ranges in \code{x} into
  \emph{absolute} ranges i.e. into ranges counted from the beginning of
  the virtual sequence obtained by concatenating all the sequences in the
  underlying genome (in the order reported by \code{seqlevels(x)}).

  \code{relativeRanges} performs the reverse transformation.

  NOTE: These functions only work on \emph{small} genomes. See Details
  section below.
}

\usage{
absoluteRanges(x)
relativeRanges(x, seqlengths)

## Related utility:
isSmallGenome(seqlengths)
}

\arguments{
  \item{x}{
    For \code{absoluteRanges}: a \link{GenomicRanges} object with ranges
    defined on a \emph{small} genome (see Details section below).

    For \code{relativeRanges}: a \link{Ranges} object.
  }
  \item{seqlengths}{
    An object holding sequence lengths. This can be a named integer
    (or numeric) vector with no duplicated names as returned by
    \code{\link[GenomeInfoDb]{seqlengths}()}, or any object from
    which sequence lengths can be extracted with
    \code{\link[GenomeInfoDb]{seqlengths}()}.

    For \code{relativeRanges}, \code{seqlengths} must describe a \emph{small}
    genome (see Details section below).
  }
}

\details{
  Because ranges in Bioconductor are using integer vectors for their
  internal representation, \code{absoluteRanges} and \code{relativeRanges}
  cannot be used on genomes for which the total sequence length is >
  \code{.Machine$integer.max}. Utility function \code{isSmallGenome} is
  provided as a mean for the user to check upfront whether the size of a
  genome is <= \code{.Machine$integer.max}, and thus compatible with
  \code{absoluteRanges} and \code{relativeRanges}.
}

\value{
  An \link[IRanges]{IRanges} object for \code{absoluteRanges}.

  A \link{GRanges} object for \code{relativeRanges}.

  \code{isSmallGenome} returns TRUE if the total length of the underlying
  sequences is <= \code{.Machine$integer.max} (e.g. Fly genome),
  FALSE if not (e.g. Human genome), or NA if it cannot be computed (because
  some sequence lengths are NA).
}

\author{
  H. Pages
}

\seealso{
  \itemize{
    \item \link{GRanges} objects.

    \item \link[IRanges]{IRanges} objects in the \pkg{IRanges} package.

    \item \link{Seqinfo} objects and the \code{\link{seqlengths}} getter in
          the \pkg{GenomeInfoDb} package.

    \item The \code{\link{tileGenome}} function for putting tiles on a
          genome.
  }
}

\examples{
## ---------------------------------------------------------------------
## TOY EXAMPLE
## ---------------------------------------------------------------------

gr <- GRanges(Rle(c("chr2", "chr1", "chr3", "chr1"), 4:1),
              IRanges(1:10, width=5),
              seqinfo=Seqinfo(c("chr1", "chr2", "chr3"), c(100, 50, 20)))

ar <- absoluteRanges(gr)
ar

gr2 <- relativeRanges(ar, seqlengths(gr))
gr2

## Sanity check:
stopifnot(all(gr == gr2))

## ---------------------------------------------------------------------
## ON REAL DATA
## ---------------------------------------------------------------------

library(TxDb.Dmelanogaster.UCSC.dm3.ensGene)
txdb <- TxDb.Dmelanogaster.UCSC.dm3.ensGene
isSmallGenome(txdb)
ex <- exons(txdb)
ex

ar <- absoluteRanges(ex)
ar

ex2 <- relativeRanges(ar, seqlengths(ex))
ex2  # original strand is not restored
strand(ex2) <- strand(ex)  # set it back
stopifnot(all(ex == ex2))

library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
isSmallGenome(txdb)
if (interactive()) {
    ex <- exons(txdb)
    absoluteRanges(ex)  # error!
}
}

\keyword{manip}
