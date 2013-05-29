\name{utils}
\docType{methods}
\alias{keepSeqlevels}
\alias{dropSeqlevels}
\alias{renameSeqlevels}
\alias{restoreSeqlevels}

\title{
  seqlevels utility functions
}

\description{
  Keep, drop or rename seqlevels in objects with a \link{Seqinfo} class.
}

\usage{
keepSeqlevels(x, value, ...)
dropSeqlevels(x, value, ...)
renameSeqlevels(x, value, ...)
restoreSeqlevels(x, ...)
}

\arguments{
  \item{x}{
    Any object having a \link{Seqinfo} class in which the seqlevels will
    be kept, dropped or renamed.
  }
  \item{value}{
    A named or unnamed character vector.

    Names are ignored by \code{keepSeqlevels} and \code{dropSeqlevels}. 
    Only the values in the character vector dictate which seqlevels to keep
    or drop.

    In the case of \code{renameSeqlevels}, the names are used to map new 
    sequence levels to the old (names correspond to the old levels). When 
    \code{value} is unnamed, the replacement vector must the same length
    and in the same order as the original \code{seqlevels(x)}.
  }
  \item{...}{
    Arguments passed to other functions.
  }
}

\value{
  The \code{x} object with seqlevels removed or renamed. If \code{x} has
  no seqlevels (empty object) or no replacement values match the current
  seqlevels in \code{x} the unchanged \code{x} is returned.
}

\details{
  Matching and overlap operations on range objects often require that the
  seqlevels match before a comparison can be made (e.g., \code{findOverlaps}).
  \code{keepSeqlevels}, \code{dropSeqlevels} and \code{renameSeqlevels} are 
  high-level convenience functions that wrap the low-level \code{seqlevels}
  function.

  \code{keepSeqlevels}, \code{dropSeqlevels}: Subsetting operations 
  that modify the size of \code{x}. \code{keepSeqlevels} keeps only the
  seqlevels in \code{value} and removes all others. \code{dropSeqlevels}
  drops the levels in \code{value} and retains all others. If \code{value} 
  does not match any seqlevels in \code{x} an empty object is returned.

  \code{renameSeqlevels}: Rename the seqlevels in \code{x} to those in 
  \code{value}. If \code{value} is a named character vector, the names are used
  to map the new seqlevels to the old. When \code{value} is unnamed, the
  replacment vector must be the same length and in the same order as the
  original \code{seqlevels(x)}.

  \code{restoreSeqlevels}: Restore the seqlevels in \code{x} back to the
  original values. Applicable only when \code{x} is a TranscriptDb. The 
  function re-initializes the TranscriptDb which resets the seqlevels, 
  removes masks and any other previous modifications.
}

\seealso{
  \itemize{
    \item \link{seqinfo} ## Accessing sequence information
    \item \link{Seqinfo} ## The Seqinfo class
  }
}

\examples{

## ---------------------------------------------------------------------
## keepSeqlevels / dropSeqlevels 
## ---------------------------------------------------------------------

## GRanges / GAlignments:

gr <- GRanges(c("chr1", "chr1", "chr2", "chr3"), IRanges(1:4, width=3))
seqlevels(gr)
## Keep only 'chr1'
chr1 <- keepSeqlevels(gr, "chr1")
## Drop 'chr1'. Both 'chr2' and 'chr3' are kept.
chr2 <- dropSeqlevels(gr, "chr1")

library(Rsamtools)
fl <- system.file("extdata", "ex1.bam", package="Rsamtools")
gal <- readGAlignments(fl)
seqlevels(gal)
## If 'value' is named, the names are ignored.
seq2 <- keepSeqlevels(gal, c(foo="seq2"))
seqlevels(seq2)

## GRangesList / GAlignmentsList:

grl <- split(gr, as.character(seqnames(gr)))
dropSeqlevels(grl, c("chr1", "chr2"))
galist <- split(gal, as.character(seqnames(gal)))
keepSeqlevels(galist, "seq2")

## TranscriptDb:

## A TranscriptDb cannot be directly subset with 'keepSeqlevels' 
## and 'dropSeqlevels'.
library(TxDb.Dmelanogaster.UCSC.dm3.ensGene)
txdb <- TxDb.Dmelanogaster.UCSC.dm3.ensGene
seqlevels(txdb)
\dontrun{
keepSeqlevels(txdb, "chr2L") ## fails
}

## GRanges or GRangesLists extracted from the TranscriptDb can be subset.
txbygene <- transcriptsBy(txdb, "gene")
seqlevels(txbygene)
chr2L <- keepSeqlevels(txbygene, "chr2L")
seqlevels(chr2L)

## ---------------------------------------------------------------------
## renameSeqlevels 
## ---------------------------------------------------------------------

## GAlignments:

seqlevels(gal)
## Rename 'seq2' to 'chr2' with a named or unnamed vector.
gal2a <- renameSeqlevels(gal, c(seq2="chr2"))
gal2b <- renameSeqlevels(gal, c("seq1", "chr2"))
## Names that do not match existing seqlevels are ignored.
## This attempt at renaming does nothing.
gal3 <- renameSeqlevels(gal, c(foo="chr2"))
identical(seqlevels(gal), seqlevels(gal3))

## TranscriptDb:

seqlevels(txdb)
## When the seqlevels of a TranscriptDb are renamed, all future 
## extractions reflect the modified seqlevels.
renameSeqlevels(txdb, sub("chr", "CH", seqlevels(txdb)))
renameSeqlevels(txdb, c(CHM="M"))
seqlevels(txdb)

transcripts <- transcripts(txdb)
identical(seqlevels(txdb), seqlevels(transcripts))

## ---------------------------------------------------------------------
## restoreSeqlevels 
## ---------------------------------------------------------------------

## Restore seqlevels in a TranscriptDb to original values.
\dontrun{
restoreSeqlevels(txdb)
seqlevels(txdb)
}
}

\author{Valerie Obenchain \email{vobencha@fhcrc.org}}

\keyword{methods}
\keyword{utilities}
