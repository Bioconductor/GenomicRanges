\name{encodeOverlaps-methods}
\alias{encodeOverlaps-methods}

\alias{encodeOverlaps}
\alias{encodeOverlaps,GRangesList,GRangesList-method}

\alias{isCompatibleWithSplicing}
\alias{isCompatibleWithSplicing,character-method}
\alias{isCompatibleWithSplicing,factor-method}
\alias{isCompatibleWithSplicing,OverlapEncodings-method}

\title{encodeOverlaps methods and related utilities}

\description{
  encodeOverlaps methods and related utilities.

  THIS IS AN EXPERIMENTAL FEATURE. STILL A WORK-IN-PROGRESS!
}

\usage{
encodeOverlaps(query, subject)
isCompatibleWithSplicing(x)
}

\arguments{
  \item{query, subject}{
    \link{GRangesList}, \link[IRanges]{RangesList} or \link[IRanges]{Ranges}
    objects.
  }
  \item{x}{
    An \link[IRanges]{OverlapEncodings} object.
  }
}

\seealso{
  \link{galn4},
  \code{\link[GenomicFeatures]{sortExonsByRank}},
  \code{\link[GenomicFeatures]{makeTranscriptDbFromUCSC}},
  \code{\link[IRanges]{findOverlaps}},
  \link{GappedAlignments-class},
  \link{GRangesList-class}
}

\examples{
## ---------------------------------------------------------------------
## A. ENCODE THE OVERLAPS RETURNED BY findOverlaps()
## ---------------------------------------------------------------------

## Load some aligned reads (RNA-seq reads aligned against the dm3
## genome, see '?galn4' for more information):
data(galn4)
reads4 <- as(galn4, "GRangesList")
table(elementLengths(reads4))  # 20316 alignments with 1 gap or more

## Retrieve transcripts and their exons from UCSC and extract the
## exons grouped by transcript in a GRangesList object.
## IMPORTANT: The reference genome of the transcripts must be *exactly*
## the same as the reference genome used to align the reads.
library(GenomicFeatures)
txdb <- makeTranscriptDbFromUCSC(genome="dm3", tablename="refGene")
exbytx0 <- exonsBy(txdb, by="tx")

## IMPORTANT: Before 'encodeOverlapsFromRangesMatching()' or
## 'encodeOverlaps()' and related utilities can be used, exons in
## transcripts located on the minus strand *must* be reordered by
## decreasing rank.
exbytx <- sortExonsByRank(exbytx0, decreasing.rank.on.minus.strand=TRUE)
exbytx0[[20]]
exbytx[[20]]

## Compute the overlaps:
overlaps <- findOverlaps(reads4, exbytx)

## Encode the overlaps:
ovenc <- encodeOverlapsFromRangesMatching(overlaps, reads4, exbytx)
ovenc
unique_ovenc <- levels(encoding(ovenc))
unique_ovenc  # 104 unique encodings
table(encoding(ovenc))

## Encodings are sort of cryptic but utilities are provided to extract
## specific meaning from them. Use of these utilities is covered below.

## ---------------------------------------------------------------------
## B. isCompatibleWithSplicing()
## ---------------------------------------------------------------------
## Use 'isCompatibleWithSplicing()' on the OverlapEncodings object to
## detect reads that overlap the transcript in a compatible way (i.e.
## in a way that is compatible with the splicing of the transcript).
## 'isCompatibleWithSplicing()' can also be used on the character vector
## containing the unique encodings to see those encodings that are
## considered compatible with the transcript.

## 7 unique encodings are compatible with the transcript:
unique_ovenc[isCompatibleWithSplicing(unique_ovenc)]

is_compat <- isCompatibleWithSplicing(ovenc)
table(is_compat)  # 241237 "compatible overlaps"

## For each alignment in 'galn4', compute the number of "compatible
## hits" i.e. the nb of transcripts with a compatible overlap:
galn4_compat_nhits <- tabulate(queryHits(overlaps)[is_compat],
                               nbins=length(galn4))
names(galn4_compat_nhits) <- names(galn4)
table(galn4_compat_nhits)
sum(galn4_compat_nhits != 0)  # 73973 alignments with at least one
                              # compatible hit

## 'galn4_compat_nhits != 0' can be used to subset either 'galn4'
## or 'reads4':
galn4[galn4_compat_nhits != 0]
reads4[galn4_compat_nhits != 0]

## For each transcript in 'exbytx' (or 'exbytx0'), compute the number
## of "compatible hits" i.e. the nb of alignments in 'galn4' (or
## 'reads4') with a compatible overlap:
exbytx_compat_nhits <- tabulate(subjectHits(overlaps)[is_compat],
                                nbins=length(exbytx))
names(exbytx_compat_nhits) <- names(exbytx)
table(exbytx_compat_nhits)

## ---------------------------------------------------------------------
## C. isCompatibleWithSkippedExons()
## ---------------------------------------------------------------------
## Coming soon...
}

\keyword{methods}
\keyword{utilities}
