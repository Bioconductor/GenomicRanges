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

    For \code{relativeRanges}: an \link{IntegerRanges} object.
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
  Because \code{absoluteRanges} returns the \emph{absolute} ranges in an
  \link[IRanges]{IRanges} object, and because an \link[IRanges]{IRanges}
  object cannot hold ranges with an end > \code{.Machine$integer.max}
  (i.e. >= 2^31 on most platforms), \code{absoluteRanges} cannot be used
  if the size of the underlying genome (i.e. the total length of the
  sequences in it) is > \code{.Machine$integer.max}. Utility function
  \code{isSmallGenome} is provided as a mean for the user to check
  upfront whether the genome is \emph{small} (i.e. its size is <=
  \code{.Machine$integer.max}) or not, and thus compatible with
  \code{absoluteRanges} or not.

  \code{relativeRanges} applies the same restriction by looking at the
  \code{seqlengths} argument.
}

\value{
  An \link[IRanges]{IRanges} object for \code{absoluteRanges}.

  A \link{GRanges} object for \code{relativeRanges}.

  \code{absoluteRanges} and \code{relativeRanges} both return an object that
  is \emph{parallel} to the input object (i.e. same length and names).

  \code{isSmallGenome} returns TRUE if the total length of the underlying
  sequences is <= \code{.Machine$integer.max} (e.g. Fly genome),
  FALSE if not (e.g. Human genome), or NA if it cannot be computed (because
  some sequence lengths are NA).
}

\author{
  H. Pagès
}

\seealso{
  \itemize{
    \item \link{GRanges} objects.

    \item \link[IRanges]{IntegerRanges} objects in the \pkg{IRanges} package.

    \item \link{Seqinfo} objects and the \code{\link{seqlengths}} getter in
          the \pkg{GenomeInfoDb} package.

    \item \link{genomicvars} for manipulating genomic variables.

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

## With a "small" genome

library(TxDb.Dmelanogaster.UCSC.dm3.ensGene)
txdb <- TxDb.Dmelanogaster.UCSC.dm3.ensGene
ex <- exons(txdb)
ex

isSmallGenome(ex)

## Note that because isSmallGenome() can return NA (see Value section
## above), its result should always be wrapped inside isTRUE() when
## used in an if statement:
if (isTRUE(isSmallGenome(ex))) {
    ar <- absoluteRanges(ex)
    ar

    ex2 <- relativeRanges(ar, seqlengths(ex))
    ex2  # original strand is not restored

    ## Sanity check:
    strand(ex2) <- strand(ex)  # restore the strand
    stopifnot(all(ex == ex2))
}

## With a "big" genome (but we can reduce it)

library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
ex <- exons(txdb)
isSmallGenome(ex)
\dontrun{
    absoluteRanges(ex)  # error!
}

## However, if we are only interested in some chromosomes, we might
## still be able to use absoluteRanges():
seqlevels(ex, pruning.mode="coarse") <- paste0("chr", 1:10)
isSmallGenome(ex)  # TRUE!
ar <- absoluteRanges(ex)
ex2 <- relativeRanges(ar, seqlengths(ex))

## Sanity check:
strand(ex2) <- strand(ex) 
stopifnot(all(ex == ex2))
}

\keyword{manip}
