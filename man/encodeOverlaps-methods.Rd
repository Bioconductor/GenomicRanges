\name{encodeOverlaps-methods}
\alias{encodeOverlaps-methods}

\alias{flipQuery}

\alias{encodeOverlaps,GRangesList,GRangesList,missing-method}

\alias{selectEncodingWithCompatibleStrand}

\alias{isCompatibleWithSplicing}
\alias{isCompatibleWithSplicing,character-method}
\alias{isCompatibleWithSplicing,factor-method}
\alias{isCompatibleWithSplicing,OverlapEncodings-method}

\alias{isCompatibleWithSkippedExons}
\alias{isCompatibleWithSkippedExons,character-method}
\alias{isCompatibleWithSkippedExons,factor-method}
\alias{isCompatibleWithSkippedExons,OverlapEncodings-method}

\alias{extractSteppedExonRanks}
\alias{extractSteppedExonRanks,character-method}
\alias{extractSteppedExonRanks,factor-method}
\alias{extractSteppedExonRanks,OverlapEncodings-method}

\alias{extractSpannedExonRanks}
\alias{extractSpannedExonRanks,character-method}
\alias{extractSpannedExonRanks,factor-method}
\alias{extractSpannedExonRanks,OverlapEncodings-method}

\alias{extractSkippedExonRanks}
\alias{extractSkippedExonRanks,character-method}
\alias{extractSkippedExonRanks,factor-method}
\alias{extractSkippedExonRanks,OverlapEncodings-method}

\alias{extractQueryStartInTranscript}

\alias{findCompatibleOverlaps}
\alias{findCompatibleOverlaps,GappedAlignments,GRangesList-method}
\alias{findCompatibleOverlaps,GappedAlignmentPairs,GRangesList-method}
\alias{countCompatibleOverlaps}


\title{encodeOverlaps method for GRangesList objects, and related utilities}

\description{
  \code{encodeOverlaps} method for \link{GRangesList}, and related utilities.
}

\usage{
\S4method{encodeOverlaps}{GRangesList,GRangesList,missing}(query, subject, hits=NULL,
               flip.query.if.wrong.strand=FALSE)

## Low-level utils:

flipQuery(x, i)

selectEncodingWithCompatibleStrand(ovencA, ovencB,
                                   query.strand, subject.strand, hits=NULL)

isCompatibleWithSplicing(x)
isCompatibleWithSkippedExons(x, max.skipped.exons=NA)

extractSteppedExonRanks(x, for.query.right.end=FALSE)
extractSpannedExonRanks(x, for.query.right.end=FALSE)
extractSkippedExonRanks(x, for.query.right.end=FALSE)

extractQueryStartInTranscript(query, subject, hits=NULL, ovenc=NULL,
                              flip.query.if.wrong.strand=FALSE,
                              for.query.right.end=FALSE)

## High-level convenience wrappers:

findCompatibleOverlaps(query, subject)
countCompatibleOverlaps(query, subject)
}

\arguments{
  \item{x}{
    For \code{flipQuery}: a \link{GRangesList} object.

    For \code{isCompatibleWithSplicing}, \code{isCompatibleWithSkippedExons},
    \code{extractSteppedExonRanks}, \code{extractSpannedExonRanks}, and
    \code{extractSkippedExonRanks}:
    an \link[IRanges]{OverlapEncodings} object, a factor or a character
    vector.
  }
  \item{i}{
    Subscript specifying the elements in \code{x} to flip. If missing, all
    the elements are flipped.
  }
  \item{ovencA, ovencB, ovenc}{
    \link[IRanges]{OverlapEncodings} objects.
  }
  \item{query, subject}{
    \link{GRangesList} objects except for \code{findCompatibleOverlaps}
    and \code{countCompatibleOverlaps} where \code{query} must be a
    \link{GappedAlignments} or \link{GappedAlignmentPairs} object.
  }
  \item{hits}{
    A \link[IRanges]{Hits} object.
    See \code{?`\link[IRanges]{encodeOverlaps,ANY,ANY,Hits-method}`} for
    a description of how a supplied \link[IRanges]{Hits} object is handled.
  }
  \item{flip.query.if.wrong.strand}{
    See the "Overlap encodings" vignette in the GenomicRanges package.
  }
  \item{query.strand, subject.strand}{
    Vector-like objects containing the strand of the query and subject,
    respectively.
  }
  \item{max.skipped.exons}{
    Not supported yet. If \code{NA} (the default), the number of skipped
    exons must be 1 or more (there is no max).
  }
  \item{for.query.right.end}{
    If \code{TRUE}, then the information reported in the output is for
    the right ends of the paired-end reads.
    Using \code{for.query.right.end=TRUE} with single-end reads is an error.
  }
}

\details{
  In the context of an RNA-seq experiment, encoding the overlaps
  between 2 \link{GRangesList} objects, one containing the reads (the query),
  and one containing the transcripts (the subject), can be used for
  detecting hits between reads and transcripts that are \emph{compatible}
  with the splicing of the transcript.

  The topic of working with overlap encodings is covered in details
  in the "Overlap encodings" vignette in the GenomicRanges package.
}

\author{
  H. Pages
}

\seealso{
  \itemize{
    \item The "Overlap encodings" vignette in the GenomicRanges package.
    \item The \code{\link[IRanges]{findOverlaps}} generic function defined
          in the IRanges package.
    \item The \link[IRanges]{OverlapEncodings} class defined in the IRanges
          package.
    \item The \link{GRangesList}, \link{GappedAlignments}, and
          \link{GappedAlignmentPairs} classes.
  }
}

\examples{
## Here we only show a simple example illustrating the use of
## countCompatibleOverlaps() on a very small data set. Please
## refer to the "Overlap encodings" vignette in the GenomicRanges
## package for a more comprehensive presentation of "overlap
## encodings" and related tools/concepts (e.g. "compatible"
## overlaps, "almost compatible" overlaps etc...), and for more
## examples.

## sm_treated1.bam contains a small subset of treated1.bam, a BAM
## file containing single-end reads from the "Pasilla" experiment
## (RNA-seq, Fly, see the pasilla data package for the details)
## and aligned to reference genome BDGP Release 5 (aka dm3 genome on
## the UCSC Genome Browser):
sm_treated1 <- system.file("extdata", "sm_treated1.bam",
                           package="GenomicRanges", mustWork=TRUE)

## Load the alignments:
library(Rsamtools)
flag0 <- scanBamFlag(isDuplicate=FALSE, isNotPassingQualityControls=FALSE)
param0 <- ScanBamParam(flag=flag0)
gal <- readGappedAlignments(sm_treated1, use.names=TRUE, param=param0)

## Load the transcripts (IMPORTANT: Like always, the reference genome
## of the transcripts must be *exactly* the same as the reference
## genome used to align the reads):
library(GenomicFeatures)
dm3_ensGene_txdb <- makeTranscriptDbFromUCSC(genome="dm3", tablename="ensGene")
dm3tx <- exonsBy(dm3_ensGene_txdb, by="tx", use.names=TRUE)

## Number of "compatible" transcripts per alignment in 'gal':
gal_ncompat <- countCompatibleOverlaps(gal, dm3tx)
elementMetadata(gal)$ncompat <- gal_ncompat
table(gal_ncompat)
mean(gal_ncompat >= 1)
## --> 33% of the alignments in 'gal' are "compatible" with at least
## 1 transcript in 'dm3tx'.

## Keep only alignments compatible with at least 1 transcript in
## 'dm3tx':
compatgal <- gal[gal_ncompat >= 1]
head(compatgal)
}

\keyword{methods}
\keyword{utilities}
