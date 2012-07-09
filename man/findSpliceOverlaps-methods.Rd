\name{findSpliceOverlaps}
\alias{findSpliceOverlaps}
\alias{findSpliceOverlaps-methods}

\alias{findSpliceOverlaps,GappedAlignments,GRangesList-method}
\alias{findSpliceOverlaps,GappedAlignmentPairs,GRangesList-method}
\alias{findSpliceOverlaps,GRangesList,GRangesList-method}

\title{findSpliceOverlaps ** UNDER CONSTRUCTION **}

\description{
  Identify reads compatible with transcript isoforms and novel splice events
}

\usage{
findSpliceOverlaps(query, subject, ignore.strand=FALSE, ...)

\S4method{findSpliceOverlaps}{GappedAlignments,GRangesList}(query, subject, ignore.strand=FALSE, ..., cds=NULL)
\S4method{findSpliceOverlaps}{GappedAlignmentPairs,GRangesList}(query, subject, ignore.strand=FALSE, ..., cds=NULL)
\S4method{findSpliceOverlaps}{GRangesList,GRangesList}(query, subject, ignore.strand=FALSE, ..., cds=NULL)

## Low-level utils:

## High-level convenience wrappers (coming soon):
#summarizeSpliceOverlaps(query, subject, ignore.strand=FALSE, ...)
}

\section{Additional methods}{
  ## Methods in GenomicFeatures : 
  \describe{
    \item{}{
      findSpliceOverlaps,GappedAlignments,TranscriptDb(query, subject, ignore.strand=FALSE, ...)
    }
    \item{}{
      findSpliceOverlaps,GappedAlignmentPairs,TranscriptDb(query, subject, ignore.strand=FALSE, ...)
    }
    \item{}{
      findSpliceOverlaps,GRangesList,TranscriptDb(query, subject, ignore.strand=FALSE, ...)
    }
  }
  ## Methods in Rsamtools : 
  \describe{
    \item{}{
      findSpliceOverlaps,character,ANY(query, subject, ignore.strand=FALSE, ..., param=ScanBamParam(), singleEnd=TRUE, cds=NULL)
    }
    \item{}{
      findSpliceOverlaps,BamFile,ANY(query, subject, ignore.strand=FALSE, ..., param=ScanBamParam(), singleEnd=TRUE, cds=NULL)
    }
  }
}

\arguments{
  \item{query}{
    \code{character} name of a Bam file, a \link[Rsamtools]{BamFile}, 
    \linkS4class{GappedAlignments}, \linkS4class{GappedAlignmentPairs} 
    or a \linkS4class{GRangesList} object containing the reads.

    Single or paired-end reads are specified with the \code{singleEnd} 
    argument (default TRUE). Paired-end reads can be supplied in a Bam 
    file or \linkS4class{GappedAlignmentPairs} object. Single-end are
    expected to be in a Bam file, \linkS4class{GappedAlignments} or 
    \linkS4class{GRanges} object.
  }
  \item{subject}{
    A \link[GenomicFeatures]{TranscriptDb}, or \link{GRangesList} containing
    the annotations. 
  }
  \item{ignore.strand}{
    When set to \code{TRUE}, strand information is ignored in the
    overlap calculations.
  }
  \item{cds}{
    Optional \link{GRangesList} of coding regions for each transcript
    in the \code{subject}. If provided, the "coding" output column
    will contain \code{logical} values indicating if the read falls 
    in a coding region.
  }
  \item{\dots}{
  }
}

\details{
  When assigning a fragment of mRNA to a transcript both the bounds
  and splicing pattern must be checked. If a read maps compatibly and
  uniquely to a transcript isoform we can quantify the expression of 
  that isoform and look for shifts in the balance of isoform expression
  When a read does not map in a compatible way to a known isoform, novel
  splice events such as splice junctions, novel exons or retentions can
  be quantified and compared across samples.

  \code{findSpliceOverlaps} detects which reads (query) match to which 
  transcripts (subject) in a compatible fashion. The non-compatible 
  reads are further analyzed to identify novel splice events and 
  transcript bounds. The output is a \link[IRanges]{Hits} object with
  metadata columns defined below. Each column is a \code{logical} 
  indicating if the read (query) met the criteria. 

  \itemize{
    \item{compatible }{Every splice (N) in a read alignment matches
      an intron in an annotated transcript. The read does not
      extend into an intron or outside the transcript bounds.
    }
    \item{unique }{The read is compatible with only one annotated 
      transcript.
    }
    \item{novelTSS }{Novel transcription site start. No matching 
      transcript has a TSS further 5' than the bound of the read.
    }
    \item{novelTSE }{Novel transcription site end. No matching 
      transcript has a TSE further 3' than the bound of the read.
    }
    \item{novelSite }{A donor or acceptor site is not present in any
      known transcript. Put another way, the start or end of a read
      splice site does not match any start or end position of a known 
      intron.
    }
    \item{novelJunction }{A combination of donor and acceptor sites is 
      not present in any known transcript. Put another way, the
      start/end combination of a read splice does not match any
      start/end combination of a known intron.
    }
    \item{novelExon }{An internal segment of a read (i.e., a segment
      separated by splices) is completely within a region that is 
      intronic in all known transcripts.
    }
    \item{novelRetention }{A read spans a region that is intronic in all 
     known transcripts.
    }
  }
}

\author{
  Michael Lawrence and Valerie Obenchain <vobencha@fhcrc.org> 
}
 
\seealso{
  \itemize{
    \item The \link{GRangesList}, \link{GappedAlignments}, and
          \link{GappedAlignmentPairs} classes.
  }
}

\examples{
  ## ----------------------------------------------------------------------- 
  ## 'subject' as TranscriptDb
  ## ----------------------------------------------------------------------- 

  ## sm_treated1.bam contains a small subset of treated1.bam, a BAM
  ## file containing single-end reads from the "Pasilla" experiment
  ## (RNA-seq, Fly, see the pasilla data package for the details).
  fl <- system.file("extdata", "sm_treated1.bam", package="GenomicRanges")
  library(TxDb.Dmelanogaster.UCSC.dm3.ensGene)
  txdb <- TxDb.Dmelanogaster.UCSC.dm3.ensGene
  hits <- findSpliceOverlaps(fl, txdb) 

  ## How many reads were compatible with a transcript?
  table(values(hits)$compatible)
  ## How many of these compatible reads map to only one transcript?
  table(values(hits)$unique)
  ## A closer look at the non-compatible reads
  ## These aren't very interesting - no splices - 
  ## need a better example
  nc <- hits[values(hits)$compatible == FALSE]
  nc

  ## ----------------------------------------------------------------------- 
  ## 'subject' as GRangesList
  ## ----------------------------------------------------------------------- 

  ## If 'cds' is not specified  when 'subject' is a GRangesList
  ## all values in the "coding" output will be NA

  exbytx <- exonsBy(txdb, by="tx", use.names=TRUE)
  hits <- findSpliceOverlaps(fl, exbytx) 

  ## Alternatively, supply the 'cds'
  cdsbytx <- cdsBy(txdb)
  hits <- findSpliceOverlaps(fl, exbytx, cds=cdsbytx) 
}

\keyword{methods}
\keyword{utilities}
