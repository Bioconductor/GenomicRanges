\name{encodeOverlaps-methods}
\alias{encodeOverlaps-methods}

\alias{encodeOverlaps,GRangesList,GRangesList-method}

\alias{isCompatibleWithSplicing}
\alias{isCompatibleWithSplicing,character-method}
\alias{isCompatibleWithSplicing,factor-method}
\alias{isCompatibleWithSplicing,OverlapEncodings-method}

\alias{isCompatibleWithSkippedExons}
\alias{isCompatibleWithSkippedExons,character-method}
\alias{isCompatibleWithSkippedExons,factor-method}
\alias{isCompatibleWithSkippedExons,OverlapEncodings-method}

\alias{extractSkippedExonRanks}
\alias{extractSkippedExonRanks,character-method}
\alias{extractSkippedExonRanks,factor-method}
\alias{extractSkippedExonRanks,OverlapEncodings-method}

\title{encodeOverlaps methods and related utilities}

\description{
  encodeOverlaps methods and related utilities.

  THIS IS AN EXPERIMENTAL FEATURE. STILL A WORK-IN-PROGRESS!
}

\usage{
\S4method{encodeOverlaps}{GRangesList,GRangesList}(query, subject)
isCompatibleWithSplicing(x)
isCompatibleWithSkippedExons(x, max.skipped.exons=NA)
extractSkippedExonRanks(x)
}

\arguments{
  \item{query, subject}{
    \link{GRangesList}, \link[IRanges]{RangesList} or \link[IRanges]{Ranges}
    objects.
  }
  \item{x}{
    An \link[IRanges]{OverlapEncodings} object, a factor or a character
    vector.
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
  \code{\link[pasillaBamSubset]{untreated1_chr4}},
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

## Load some aligned reads (RNA-seq single-end reads aligned against the
## dm3 genome, see '?untreated1_chr4' in the pasillaBamSubset data
## package for more information):
library(pasillaBamSubset)
gal4 <- readGappedAlignments(untreated1_chr4(), use.names=TRUE)
reads4 <- as(gal4, "GRangesList")
table(elementLengths(reads4))  # 20316 alignments with 1 gap or more

## Retrieve transcripts and their exons from UCSC and extract the
## exons grouped by transcript in a GRangesList object.
## IMPORTANT: The reference genome of the transcripts must be *exactly*
## the same as the reference genome used to align the reads.
library(GenomicFeatures)
txdb <- makeTranscriptDbFromUCSC(genome="dm3", tablename="refGene")
exbytx0 <- exonsBy(txdb, by="tx")

## IMPORTANT: Before 'encodeOverlaps()' or related utilities can be used,
## exons in transcripts located on the minus strand *must* be reordered by
## decreasing rank.
exbytx <- sortExonsByRank(exbytx0, decreasing.rank.on.minus.strand=TRUE)
exbytx0[[20]]
exbytx[[20]]

## Compute the overlaps:
overlaps <- findOverlaps(reads4, exbytx)

## Encode the overlaps:
ovenc <- encodeOverlaps(reads4, exbytx, overlaps)
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
## overlaps the transcript in a "compatible" way, that is, in a way
## compatible with the splicing of the transcript. We call this an
## overlap of type COMPATIBLE_WITH_SPLICING.
## Use 'isCompatibleWithSplicing()' on the OverlapEncodings object to
## detect overlaps of type COMPATIBLE_WITH_SPLICING.
## 'isCompatibleWithSplicing()' can also be used on the character vector
## containing the unique encodings to find those of type
## COMPATIBLE_WITH_SPLICING i.e. those corresponding to overlaps of this
## type.

## 7 unique encodings are of type COMPATIBLE_WITH_SPLICING:
table(encoding(ovenc))[isCompatibleWithSplicing(unique_ovenc)]

## Encodings "1:i:" (204386 occurences in 'ovenc'), "2:jm:af:" (35236
## occurences) and "3:jmm:agm:aaf:" (267 occurences) correspond to the
## following overlaps:
##
## "1:i:"
##     read (no gap):            oooooo
##        transcript:   ...  oooooooooooo ...
##
## "2:jm:af:"
##      read (1 gap):              oooo------oo
##        transcript:   ...  oooooooooo      ooooooo ...
##
## "3:jmm:agm:aaf:"
##     read (2 gaps):            ooo-----oo---------ooo
##        transcript:   ...  ooooooo     oo         oooooooooo ...
##
## Note that the exons represented in the 2nd and 3rd drawings are
## consecutive in the transcript. The transcript can have more exons,
## which is denoted by the ... on each side of the drawings.

is_compat <- isCompatibleWithSplicing(ovenc)
table(is_compat)  # 241237 overlaps of type COMPATIBLE_WITH_SPLICING

## For each alignment in 'gal4', compute the number of hits (i.e.
## overlaps) of type COMPATIBLE_WITH_SPLICING:
gal4_compat_nhits <- tabulate(queryHits(overlaps)[is_compat],
                              nbins=length(gal4))
table(gal4_compat_nhits)
sum(gal4_compat_nhits != 0)  # 73973 alignments with at least one hit
                             # of type COMPATIBLE_WITH_SPLICING

## 'gal4_compat_nhits != 0' can be used to subset either 'gal4'
## or 'reads4':
gal4[gal4_compat_nhits != 0]
reads4[gal4_compat_nhits != 0]

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

is_compat2 <- isCompatibleWithSkippedExons(ovenc)
table(is_compat2)  # 484 overlaps of type COMPATIBLE_WITH_SKIPPED_EXONS

## For each alignment in 'gal4', compute the number of hits (i.e.
## overlaps) of type COMPATIBLE_WITH_SKIPPED_EXONS:
gal4_compat2_nhits <- tabulate(queryHits(overlaps)[is_compat2],
                               nbins=length(gal4))
table(gal4_compat2_nhits)
sum(gal4_compat2_nhits != 0)  # 216 alignments with at least one hit
                              # of type COMPATIBLE_WITH_SKIPPED_EXONS

## For each transcript in 'exbytx' (or 'exbytx0'), compute the number
## of hits of type COMPATIBLE_WITH_SKIPPED_EXONS:
exbytx_compat2_nhits <- tabulate(subjectHits(overlaps)[is_compat2],
                                 nbins=length(exbytx))
names(exbytx_compat2_nhits) <- names(exbytx)
table(exbytx_compat2_nhits)

## ---------------------------------------------------------------------
## D. COMBINING RESULTS OF isCompatibleWithSplicing() AND
##    isCompatibleWithSkippedExons() TO DETECT NOVEL SPLICE JUNCTIONS
## ---------------------------------------------------------------------
## An alignment in 'gal4' with hits of type
## COMPATIBLE_WITH_SKIPPED_EXONS but no hits of type
## COMPATIBLE_WITH_SPLICING suggests the presence of a transcript that
## is not in our annotations.

## First we extract the index of those alignments:
aln_shows_nov_splice_jct <- gal4_compat2_nhits != 0L &
                            gal4_compat_nhits == 0L
which(aln_shows_nov_splice_jct)

## We make this an index into 'overlaps' (Hits object):
is_nov_splice_jct <- queryHits(overlaps) %in%
                         which(aln_shows_nov_splice_jct)

## We intersect with 'is_compat2' to keep only the overlaps of interest:
is_nov_splice_jct <- is_nov_splice_jct & is_compat2

## For each overlap of interest, we extract the ranks of the skipped
## exons (we use a list for this as there might be more than 1 skipped
## exon per overlap):
skpexrk <- extractSkippedExonRanks(ovenc)[is_nov_splice_jct]
table(elementLengths(skpexrk))

## Finally, we split 'skpexrk' by transcript TxDb internal id.
names(skpexrk) <- queryHits(overlaps)[is_nov_splice_jct]
f <- names(exbytx)[subjectHits(overlaps)[is_nov_splice_jct]]
tx2skpexrk <- split(skpexrk, f)

## 'tx2skpexrk' is a named list of named lists of integer vectors.
## The first level of names (outer names) are transcript TxDb internal
## ids and the second level of names (inner names) are alignment
## indices into 'gal4':
names(tx2skpexrk)  # transcript TxDb internal ids

## Transcript "10" has 4 hits that skip exons 9 and 10:
tx2skpexrk[["10"]]
## Transcript "58" has 7 hits. One hit skips exons 6 to 10, five hits
## skip exons 9 and 10, and one hit skips exon 10 only:
tx2skpexrk[["58"]]


## FINAL NOTES:
##
## (1) The exon ranks returned by extractSkippedExonRanks() are exon
##     positions within 'exbytx', not 'exbytx0'. For exons on the plus
##     strand, it makes no difference i.e. the "real" exon ranks are
##     returned. But for exons on the minus strand, we get the
##     "anti-rank" i.e. N + 1 - rk where N is the total nb of exons in
##     the transcript and rk the "real" exon rank. So for example in
##     the case of transcript "58" (located on the minus strand, see
##     'exbytx[["58"]]'), exon with "anti-rank" 10 is actually at rank
##     2 (the transcript has 11 exons).
##
## (2) Let's be cautious about the interpretation of 'tx2skpexrk'.
##     Because of how we've conducted this analysis, the aligments
##     reported in 'tx2skpexrk' are guaranteed to NOT have any overlaps
##     of type COMPATIBLE_WITH_SPLICING with other known transcripts.
##     All we can say, for example in the case of transcript "10", is
##     that the 4 reported hits show evidence of one or more unknown
##     transcripts with a splice junction that corresponds to the gap
##     between exons 8 and 11. But without further analysis, we can't
##     make any assumption about the exons structure of those unknown
##     transcripts. In particular, we cannot assume the existence of an
##     unknown transcript made of the same exons as transcript "10"
##     minus exons 9 and 10!
}

\keyword{methods}
\keyword{utilities}
