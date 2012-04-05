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

\alias{extractSpannedExonRanks}
\alias{extractSpannedExonRanks,character-method}
\alias{extractSpannedExonRanks,factor-method}
\alias{extractSpannedExonRanks,OverlapEncodings-method}

\alias{extractSkippedExonRanks}
\alias{extractSkippedExonRanks,character-method}
\alias{extractSkippedExonRanks,factor-method}
\alias{extractSkippedExonRanks,OverlapEncodings-method}

\alias{extractQueryStartInTranscript}

\title{encodeOverlaps methods and related utilities}

\description{
  \code{encodeOverlaps} methods and related utilities.
}

\usage{
flipQuery(x, i)

\S4method{encodeOverlaps}{GRangesList,GRangesList,missing}(query, subject, hits=NULL,
               flip.query.if.wrong.strand=FALSE)

selectEncodingWithCompatibleStrand(x, y,
                                   query.strand, subject.strand, hits=NULL)

isCompatibleWithSplicing(x)
isCompatibleWithSkippedExons(x, max.skipped.exons=NA)

extractSpannedExonRanks(x)
extractSkippedExonRanks(x)

extractQueryStartInTranscript(query, subject, hits=NULL, ovenc=NULL,
                              flip.query.if.wrong.strand=FALSE)
}

\arguments{
  \item{x}{
    For \code{flipQuery}: A \link{GRangesList} object.

    For \code{selectEncodingWithCompatibleStrand}:
    An \link[IRanges]{OverlapEncodings} object.

    For \code{isCompatibleWithSplicing}, \code{isCompatibleWithSkippedExons},
    and \code{extractSkippedExonRanks}:
    An \link[IRanges]{OverlapEncodings} object, a factor or a character
    vector.
  }
  \item{i}{
    Subscript specifying the elements in \code{x} to flip. If missing, all
    the elements are flipped.
  }
  \item{query, subject}{
    \link{GRangesList} objects.
  }
  \item{hits}{
    A \link[IRanges]{Hits} object.
    See \code{?`\link[IRanges]{encodeOverlaps,ANY,ANY,Hits-method}`} for
    a description of how a supplied \link[IRanges]{Hits} object is handled.
  }
  \item{flip.query.if.wrong.strand}{
    See the "Overlap encodings" vignette in the GenomicRanges package.
  }
  \item{y}{
    An \link[IRanges]{OverlapEncodings} object.
  }
  \item{query.strand, subject.strand}{
    Vector-like objects containing the strand of the query and subject,
    respectively.
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
    \item \code{\link[IRanges]{findOverlaps}}.
    \item \link[IRanges]{OverlapEncodings-class}.
    \item \link{GRangesList-class}.
  }
}

\examples{
## See the "Overlap encodings" vignette in the GenomicRanges package for
## some examples.
}

\keyword{methods}
\keyword{utilities}
