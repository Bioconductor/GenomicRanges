\name{encodeOverlaps-methods}
\alias{encodeOverlaps-methods}

\alias{encodeOverlaps}
\alias{encodeOverlaps,GRangesList,GRangesList-method}

\alias{isCompatibleWithSplicing}
\alias{isCompatibleWithSplicing,character-method}
\alias{isCompatibleWithSplicing,factor-method}
\alias{isCompatibleWithSplicing,OverlapEncodings-method}

\alias{isCompatibleWithSkippedExons}
\alias{isCompatibleWithSkippedExons,character-method}
\alias{isCompatibleWithSkippedExons,factor-method}
\alias{isCompatibleWithSkippedExons,OverlapEncodings-method}

\title{encodeOverlaps methods and related utilities}

\description{
  encodeOverlaps methods and related utilities.

  THIS IS AN EXPERIMENTAL FEATURE. STILL A WORK-IN-PROGRESS!
}

\usage{
encodeOverlaps(query, subject)
isCompatibleWithSplicing(x)
isCompatibleWithSkippedExons(x, max.skipped.exons=NA)
}

\arguments{
  \item{query, subject}{
    \link{GRangesList}, \link[IRanges]{RangesList} or \link[IRanges]{Ranges}
    objects.
  }
  \item{x}{
    An \link[IRanges]{OverlapEncodings} object.
  }
  \item{max.skipped.exons}{
    Not supported yet. If \code{NA} (the default), the number of skipped
    exons must be 1 or more (there is no max).
  }
}

\author{
  H. Pages
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
## We are interested in a particular type of overlap where the read
## overlaps the transcript in a "compatible" way, that is, in a way that
## is compatible with the splicing of the transcript. We call this
## an overlap of type COMPATIBLE_WITH_SPLICING.
## Use 'isCompatibleWithSplicing()' on the OverlapEncodings object to
## detect overlaps of type COMPATIBLE_WITH_SPLICING.
## 'isCompatibleWithSplicing()' can also be used on the character vector
## containing the unique encodings to find those of type
## COMPATIBLE_WITH_SPLICING i.e. those corresponding to overlaps of this
## type.

## 7 unique encodings are of type COMPATIBLE_WITH_SPLICING:
table(encoding(ovenc))[isCompatibleWithSplicing(unique_ovenc)]

is_compat <- isCompatibleWithSplicing(ovenc)
table(is_compat)  # 241237 overlaps of type COMPATIBLE_WITH_SPLICING

## For each alignment in 'galn4', compute the number of hits (i.e.
## overlaps) of type COMPATIBLE_WITH_SPLICING:
galn4_compat_nhits <- tabulate(queryHits(overlaps)[is_compat],
                               nbins=length(galn4))
names(galn4_compat_nhits) <- names(galn4)
table(galn4_compat_nhits)
sum(galn4_compat_nhits != 0)  # 73973 alignments with at least one hit
                              # of type COMPATIBLE_WITH_SPLICING

## 'galn4_compat_nhits != 0' can be used to subset either 'galn4'
## or 'reads4':
galn4[galn4_compat_nhits != 0]
reads4[galn4_compat_nhits != 0]

## For each transcript in 'exbytx' (or 'exbytx0'), compute the number
## of hits of type COMPATIBLE_WITH_SPLICING:
exbytx_compat_nhits <- tabulate(subjectHits(overlaps)[is_compat],
                                nbins=length(exbytx))
names(exbytx_compat_nhits) <- names(exbytx)
table(exbytx_compat_nhits)

## ---------------------------------------------------------------------
## C. isCompatibleWithSkippedExons()
## ---------------------------------------------------------------------
## In many aspects, overlaps of type COMPATIBLE_WITH_SPLICING correspond
## to an ideal situation but in practise many reads don't fall into that
## category.
## Here we are interested in a "less perfect" type of overlap where the
## read overlaps the transcript in a way that *would* be compatible if 1
## or more exons were removed from the transcript. We call this an
## overlap of type COMPATIBLE_WITH_SKIPPED_EXONS.
## Use 'isCompatibleWithSkippedExons()' on the OverlapEncodings object
## to detect overlaps of type COMPATIBLE_WITH_SKIPPED_EXONS.
## 'isCompatibleWithSkippedExons()' can also be used on the character
## vector containing the unique encodings to find those of type
## COMPATIBLE_WITH_SKIPPED_EXONS i.e. those corresponding to overlaps of
## this type.

## 8 unique encodings are of type COMPATIBLE_WITH_SKIPPED_EXONS:
table(encoding(ovenc))[isCompatibleWithSkippedExons(unique_ovenc)]

## Encoding "2:jm:am:af:" (398 occurences in 'ovenc') corresponds to the
## following overlap:
##
##   read (1 gap):          ooooo------------ooo
##     transcript:   ...  ooooooo    oooo    oooooooo ...
## Note that the transcript can have more than the 3 exons represented
## here (hence the ... on each side of the 3 exons).

is_compat2 <- isCompatibleWithSkippedExons(ovenc)
table(is_compat2)  # 484 overlaps of type COMPATIBLE_WITH_SKIPPED_EXONS

## For each alignment in 'galn4', compute the number of hits (i.e.
## overlaps) of type COMPATIBLE_WITH_SKIPPED_EXONS:
galn4_compat2_nhits <- tabulate(queryHits(overlaps)[is_compat2],
                                nbins=length(galn4))
names(galn4_compat2_nhits) <- names(galn4)
table(galn4_compat2_nhits)
sum(galn4_compat2_nhits != 0)  # 216 alignments with at least one hit
                               # of type COMPATIBLE_WITH_SKIPPED_EXONS

## For each transcript in 'exbytx' (or 'exbytx0'), compute the number
## of hits of type COMPATIBLE_WITH_SKIPPED_EXONS:
exbytx_compat2_nhits <- tabulate(subjectHits(overlaps)[is_compat2],
                                 nbins=length(exbytx))
names(exbytx_compat2_nhits) <- names(exbytx)
table(exbytx_compat2_nhits)
}

\keyword{methods}
\keyword{utilities}
