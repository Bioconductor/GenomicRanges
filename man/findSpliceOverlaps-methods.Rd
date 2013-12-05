\name{findSpliceOverlaps-methods}

\alias{findSpliceOverlaps-methods}

\alias{findSpliceOverlaps}
\alias{findSpliceOverlaps,GRangesList,GRangesList-method}


\title{Classify ranges (reads) as compatible with existing genomic annotations 
       or as having novel splice events}

\description{
  The \code{findSpliceOverlaps} function identifies ranges (reads) that are 
  compatible with a specific transcript isoform. The non-compatible ranges are
  analyzed for the presence of novel splice events.

  NOTE: \code{findSpliceOverlaps} methods for \link{GAlignments},
  \link[GenomicAlignments]{GAlignmentPairs}, and \link[Rsamtools]{BamFile}
  objects are defined and documented in the \pkg{GenomicAlignments} package.
  Methods for \link[GenomicFeatures]{TranscriptDb} objects are defined and
  documented in the \pkg{GenomicFeatures} package.
}

\usage{
findSpliceOverlaps(query, subject, ignore.strand=FALSE, ...)

\S4method{findSpliceOverlaps}{GRangesList,GRangesList}(query, subject, ignore.strand=FALSE, ..., cds=NULL)

## Low-level utils:

## High-level convenience wrappers (coming soon):
#summarizeSpliceOverlaps(query, subject, ignore.strand=FALSE, ...)
}

\arguments{
  \item{query}{
    A \link{GRangesList} object containing the reads.
  }
  \item{subject}{
    A \link{GRangesList} object containing the annotations. This list is
    expected to be exons by transcripts. 
  }
  \item{ignore.strand}{
    When set to \code{TRUE}, strand information is ignored in the
    overlap calculations.
  }
  \item{...}{
    Additional arguments used in other methods.
  }
  \item{cds}{
    Optional \link{GRangesList} of coding regions for each transcript
    in the \code{subject}. If provided, the "coding" output column
    will be a \code{logical} vector indicating if the read falls 
    in a coding region. When not provided, the "coding" output is 
    \code{NA}.
  }
}

\details{
  When a read maps compatibly and uniquely to a transcript isoform we 
  can quantify the expression and look for shifts in the balance of 
  isoform expression. If a read does not map in compatible way, novel 
  splice events such as splice junctions, novel exons or retentions 
  can be quantified and compared across samples.

  \code{findSpliceOverlaps} detects which reads (query) match to 
  transcripts (subject) in a compatible fashion. Compatibility is based 
  on both the transcript bounds and splicing pattern. Assessing the 
  splicing pattern involves comparision of the read splices (i.e., the 
  "N" gaps in the cigar) with the transcript introns. For paired-end 
  reads, the inter-read gap is not considered a splice. The analysis 
  of non-compatible reads for novel splice events is under construction.
}

\value{
  The output is a \link[IRanges]{Hits} object with
  the metadata columns defined below. Each column is a \code{logical} 
  indicating if the read (query) met the criteria. 

  \itemize{
    \item{compatible: }{Every splice (N) in a read alignment matches
      an intron in an annotated transcript. The read does not
      extend into an intron or outside the transcript bounds.
    }
    \item{unique: }{The read is compatible with only one annotated 
      transcript.
    }
    \item{strandSpecific: }{The query (read) was stranded. 
    }
  }
}

\author{
  Michael Lawrence and Valerie Obenchain <vobencha@fhcrc.org> 
}
 
\seealso{
  \itemize{
    \item The \link{GRangesList} class.
    \item \link[GenomicAlignments]{findSpliceOverlaps-methods} in the
          \pkg{GenomicAlignments} package.
    \item \link[GenomicFeatures]{findSpliceOverlaps-methods} in the
          \pkg{GenomicFeatures} package.
  }
}

\examples{
## See ?`findSpliceOverlaps-methods` in the GenomicAlignments package
## for some examples.
}

\keyword{methods}
\keyword{utilities}
