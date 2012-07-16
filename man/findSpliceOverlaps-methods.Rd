\name{findSpliceOverlaps}
\alias{findSpliceOverlaps}
\alias{findSpliceOverlaps-methods}

\alias{findSpliceOverlaps,GappedAlignments,GRangesList-method}
\alias{findSpliceOverlaps,GappedAlignmentPairs,GRangesList-method}
\alias{findSpliceOverlaps,GRangesList,GRangesList-method}

\title{Classify ranges (reads) as compatible with existing genomic annotations 
       or as having novel splice events}

\description{
  The \code{findSpliceOverlaps} function identifies ranges (reads) that are 
  compatible with a specific transcript isoform. The non-compatible ranges are
  analyzed for the presence of novel splice events.
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
    A \link{GRangesList} containing the annotations. This list is
    expected to be exons by transcripts. 
  }
  \item{ignore.strand}{
    When set to \code{TRUE}, strand information is ignored in the
    overlap calculations.
  }
  \item{cds}{
    Optional \link{GRangesList} of coding regions for each transcript
    in the \code{subject}. If provided, the "coding" output column
    will be a \code{logical} vector indicating if the read falls 
    in a coding region. When not provided, the "coding" output is 
    \code{NA}.
  }
  \item{\dots}{
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
  reads, the inter-read gap is not considered a splice.

  Reads identified as non-compatible are further analyzed for novel 
  splice events. The output is a \link[IRanges]{Hits} object with
  elementMetadata columns defined below. Each column is a \code{logical} 
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
      Applies to reads with or without gaps.
    }
    \item{novelTSE }{Novel transcription site end. No matching 
      transcript has a TSE further 3' than the bound of the read.
      Applies to reads with or without gaps.
    }
    \item{novelSite }{A donor or acceptor site is not present in any
      known transcript. Put another way, the start or end of a read
      splice site does not match any start or end position of a known 
      intron. Applies to reads with gaps only.
    }
    \item{novelJunction }{A combination of donor and acceptor sites is 
      not present in any known transcript. Put another way, the
      start/end combination of a read splice does not match any
      start/end combination of a known intron. Applies to reads
      with gaps only.
    }
    \item{novelExon }{An internal segment of a read (i.e., a segment
      separated by splices) is completely within a region that is 
      intronic in all known transcripts. Applies to reads with
      gaps only.
    }
    \item{novelRetention }{A read fragment overlaps part or all of 
     a region that is intronic in all known transcripts. Applies 
     to reads with or without gaps.
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
  ## Isoform expression :
  ## ----------------------------------------------------------------------- 
  ## findSpliceOverlaps() can assist in quantifying isoform expression
  ## by identifying reads that map compatibly and uniquely to a 
  ## transcript isoform.
  \dontrun{
  library(Rsamtools)
  library(TxDb.Dmelanogaster.UCSC.dm3.ensGene)
  txdb <- TxDb.Dmelanogaster.UCSC.dm3.ensGene
  library(pasillaBamSubset)
  ## untreated1_chr4() contains single-end reads
  se <- untreated1_chr4()
  exbytx <- exonsBy(txdb, "tx")
  cdsbytx <- cdsBy(txdb, "tx")
  param <- ScanBamParam(which=GRanges("chr4", IRanges(1e5,3e5)))
  sehits <- findSpliceOverlaps(se, exbytx, cds=cdsbytx, param=param)

  ## Tally the reads by category to get an idea of read distribution.
  lst <- lapply(values(sehits), table)
  nms <- names(lst) %in% c("compatible", "unique")
  tbl <- do.call(rbind, lst[nms]) 
  tbl

  ## Reads compatible with one or more transcript isoforms.
  rnms <- rownames(tbl)
  tbl[rnms == "compatible","TRUE"]/sum(tbl[rnms == "compatible",])

  ## Reads compatible with a single isoform. 
  tbl[rnms == "unique","TRUE"]/sum(tbl[rnms == "unique",])

  ## All reads fall in a coding region as defined by 
  ## the txdb annotation.
  lst[["coding"]] 

  ## Check : Total number of reads should be the same across categories.
  lapply(lst, sum)
 
  ## ----------------------------------------------------------------------- 
  ## Novel splice events :
  ## ----------------------------------------------------------------------- 
  ## When a read does not map in a compatible way to a known isoform, novel
  ## splice events such as splice junctions, novel exons or retentions may
  ## be present. 

  ## Subset the data on the non-compatible reads and summarize
  ## the categories.
  nc <- sehits[values(sehits)$compatible == FALSE]
  lapply(values(nc), table)

  ## Isolate reads with a novelTSS or novelTSE
  novelTSS <- sehits[values(sehits)$novelTSS]
  novelTSE <- sehits[values(sehits)$novelTSE]
  ## Locate the annotation region of the novelTSS reads.
  exbytx <- exonsBy(txdb, "tx")
  exbytx[subjectHits(novelTSS)]

  ## None of the reads with a novelSite or novelJunction map to 
  ## only one isoform.
  sehits[values(sehits)$novelSite & values(sehits)$unique]
  sehits[values(sehits)$novelJunction & values(sehits)$unique]

  ## ----------------------------------------------------------------------- 
  ## Paired-end reads :
  ## ----------------------------------------------------------------------- 
  ## The 'singleEnd' argument is set to FALSE for a Bam file with
  ## paired-end reads.
  \dontrun{
  pe <- untreated3_chr4()
  hits2 <- findSpliceOverlaps(pe, txdb, singleEnd=FALSE, param=param)
  } 
  ## In addition to Bam files, paired-end reads can be supplied in a 
  ## GappedAlignmentPairs object. In this example the difference between 
  ## novelJunction and novelSite is demonstrated. All novelSites result in 
  ## novelJunctions but not vice versa.
  genes <- GRangesList(
      GRanges("chr1", IRanges(c(5, 20), c(10, 25)), "+"),
      GRanges("chr1", IRanges(c(5, 22), c(15, 25)), "+"))
  galp <- GappedAlignmentPairs(
      GappedAlignments("chr1", 5L, "11M4N6M", strand("+")),
      GappedAlignments("chr1", 50L, "6M", strand("-")),
      isProperPair=TRUE)
  findSpliceOverlaps(galp, genes)
}
}

\keyword{methods}
\keyword{utilities}
