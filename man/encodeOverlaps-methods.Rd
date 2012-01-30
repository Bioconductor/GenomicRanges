\name{encodeOverlaps-methods}
\alias{encodeOverlaps-methods}

\alias{encodeOverlaps}
\alias{isCompatibleWithSplicing}
\alias{isCompatibleWithSplicing,character-method}
\alias{isCompatibleWithSplicing,factor-method}
\alias{isCompatibleWithSplicing,OverlapEncodings-method}

\title{encodeOverlaps methods and related utilities}

\description{
  encodeOverlaps methods and related utilities.

  THIS IS AN EXPERIMENTAL FEATURE AND STILL A WORK-IN-PROGRESS!
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
## Load some aligned reads (RNA-seq reads aligned against the dm3
## genome, see ?galn4 for more information):
data(galn4)
reads4 <- as(galn4, "GRangesList")
table(elementLengths(reads4))  # 20316 alignments with 1 gap or more

## Retrieve transcripts and their exons from UCSC and extract the
## exons grouped by transcript in a GRangesList object:
library(GenomicFeatures)
txdb <- makeTranscriptDbFromUCSC(genome="dm3", tablename="refGene")
exbytx0 <- exonsBy(txdb, by="tx")

## IMPORTANT: Before encodeOverlaps()/encodeOverlapsFromRangesMatching()
## and related utilities can be used, exons *must* be reordered by
## decreasing rank in transcripts located on the minus strand:
exbytx <- sortExonsByRank(exbytx0, decreasing.rank.on.minus.strand=TRUE)
exbytx0[[20]]
exbytx[[20]]

## Compute the overlaps:
overlaps <- findOverlaps(reads4, exbytx)

## Compute the overlap encodings:
ovenc <- encodeOverlapsFromRangesMatching(overlaps, reads4, exbytx)
table(encoding(ovenc))
table(isCompatibleWithSplicing(ovenc))
}

\keyword{methods}
\keyword{utilities}
