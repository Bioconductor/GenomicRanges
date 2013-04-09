\name{findSpliceOverlaps}
\alias{findSpliceOverlaps}
\alias{findSpliceOverlaps-methods}

\alias{findSpliceOverlaps,GAlignments,GRangesList-method}
\alias{findSpliceOverlaps,GAlignmentPairs,GRangesList-method}
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

\S4method{findSpliceOverlaps}{GAlignments,GRangesList}(query, subject, ignore.strand=FALSE, ..., cds=NULL)
\S4method{findSpliceOverlaps}{GAlignmentPairs,GRangesList}(query, subject, ignore.strand=FALSE, ..., cds=NULL)
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
    \linkS4class{GAlignments}, \linkS4class{GAlignmentPairs} 
    or a \linkS4class{GRangesList} object containing the reads.

    Single or paired-end reads are specified with the \code{singleEnd} 
    argument (default FALSE). Paired-end reads can be supplied in a Bam 
    file or \linkS4class{GAlignmentPairs} object. Single-end are
    expected to be in a Bam file, \linkS4class{GAlignments} or 
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
  reads, the inter-read gap is not considered a splice. The analysis 
  of non-compatible reads for novel splice events is under construction.

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
    \item The \link{GRangesList}, \link{GAlignments}, and
          \link{GAlignmentPairs} classes.
  }
}

\examples{
## ----------------------------------------------------------------------- 
## Isoform expression :
## ----------------------------------------------------------------------- 
## findSpliceOverlaps() can assist in quantifying isoform expression
## by identifying reads that map compatibly and uniquely to a 
## transcript isoform.
library(Rsamtools)
library(TxDb.Dmelanogaster.UCSC.dm3.ensGene)
library(pasillaBamSubset)
se <- untreated1_chr4() ## single-end reads
txdb <- TxDb.Dmelanogaster.UCSC.dm3.ensGene
exbytx <- exonsBy(txdb, "tx")
cdsbytx <- cdsBy(txdb, "tx")
param <- ScanBamParam(which=GRanges("chr4", IRanges(1e5,3e5)))
sehits <- findSpliceOverlaps(se, exbytx, cds=cdsbytx, param=param)

## Tally the reads by category to get an idea of read distribution.
lst <- lapply(mcols(sehits), table)
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
## Paired-end reads :
## ----------------------------------------------------------------------- 
## 'singleEnd' is set to FALSE for a Bam file with paired-end reads.
pe <- untreated3_chr4()
hits2 <- findSpliceOverlaps(pe, exbytx, singleEnd=FALSE, param=param)
 
## In addition to Bam files, paired-end reads can be supplied in a 
## GAlignmentPairs object. 
genes <- GRangesList(
    GRanges("chr1", IRanges(c(5, 20), c(10, 25)), "+"),
    GRanges("chr1", IRanges(c(5, 22), c(15, 25)), "+"))
galp <- GAlignmentPairs(
    GAlignments("chr1", 5L, "11M4N6M", strand("+")),
    GAlignments("chr1", 50L, "6M", strand("-")),
    isProperPair=TRUE)
findSpliceOverlaps(galp, genes)
}

\keyword{methods}
\keyword{utilities}
