\name{encodeOverlaps-methods}
\alias{encodeOverlaps-methods}

\alias{encodeOverlaps,GRangesList,GRangesList,missing-method}
\alias{encodeOverlaps,GappedAlignments,GRangesList,missing-method}
\alias{encodeOverlaps,GappedAlignmentPairs,GRangesList,missing-method}

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
  \code{encodeOverlaps} methods and related utilities.
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
    \item \code{\link[GenomicFeatures]{sortExonsByRank}}.
    \item \code{\link[GenomicFeatures]{makeTranscriptDbFromUCSC}}.
    \item \code{\link[IRanges]{findOverlaps}}.
    \item \link{GappedAlignments-class}.
    \item \link{GappedAlignmentPairs-class}.
    \item \link{GRangesList-class}.
  }
}

\examples{
## See the "Overlap encodings" vignette in the GenomicRanges package for
## some examples.
}

\keyword{methods}
\keyword{utilities}
