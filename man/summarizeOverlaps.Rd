\name{summarizeOverlaps}

\alias{summarizeOverlaps}
\alias{summarizeOverlaps,GRanges,GappedAlignments-method}
\alias{summarizeOverlaps,GRangesList,GappedAlignments-method}
\alias{summarizeOverlaps,GRanges,GappedAlignmentPairs-method}
\alias{summarizeOverlaps,GRangesList,GappedAlignmentPairs-method}
\alias{Union}
\alias{IntersectionStrict}
\alias{IntersectionNotEmpty}


\title{Perform overlap queries between reads and genomic features} 

\description{
  \code{summarizeOverlaps} extends \code{findOverlaps} by providing 
  options to resolve reads that overlap multiple features. Each read 
  is counted a maximum of once.
}

\usage{
  \S4method{summarizeOverlaps}{GRanges,GappedAlignments}(
    features, reads, mode, ignore.strand=FALSE, ...) 
  \S4method{summarizeOverlaps}{GRangesList,GappedAlignments}(
    features, reads, mode, ignore.strand=FALSE, ...) 
  \S4method{summarizeOverlaps}{GRanges,GappedAlignmentPairs}(
    features, reads, mode, ignore.strand=FALSE, ...) 
  \S4method{summarizeOverlaps}{GRangesList,GappedAlignmentPairs}(
    features, reads, mode, ignore.strand=FALSE, ...) 
}

\arguments{
  \item{features}{
    A \linkS4class{GRanges} or a \linkS4class{GRangesList} containing genomic 
    regions of interest. This will commonly be a \code{GRanges} of exons or
    transcripts or a \linkS4class{GRangesList} of exons by gene or transcripts
    by gene.

  }
  \item{reads}{
    A \linkS4class{GappedAlignments} (single-end), 
    \linkS4class{GappedAlignmentPairs} (paired-end), or
    \link[Rsamtools]{BamFileList} or \link[Rsamtools]{BamViews} objects
    containing the reads to be mapped to the genomic regions of interest.

    The \link[Rsamtools]{BamFileList} and \link[Rsamtools]{BamViews} methods 
    can handle single or paired-end reads by appropriately specifying the
    \code{singleEnd} argument. Currently \code{singleEnd} is set once for 
    all files in the list (i.e., the list cannot contain both single and 
    paired-end read files.). 
    See ?\code{summerizeOverlaps,GRanges,BamFileList-method}.
  }
  \item{mode}{
    Character name of a function that defines the counting method to be used. 
    Modes include `Union', `IntersectionStrict', or `IntersectionNotEmpty' and 
    are designed after the counting modes in the HTSeq package by Simon Anders 
    (see references). All methods are gap-aware and count a read a maximum of 
    once. A user can provide their own count function as the \code{mode} 
    argument when using the \link[Rsamtools]{BamFileList} or 
    \link[Rsamtools]{BamViews} methods.

    \itemize{
      \item Union : (Default) Reads that overlap any portion of exactly one 
            feature are counted. Reads that overlap multiple features are 
            discarded. For mode `Union' gapped reads are handled the same as 
            simple reads. If any portion of the gapped read hits >1 feature 
            the read is discarded. 

            The number of reads counted depends on the quality of the features
            (i.e., do the features overlap, have gaps, etc.). Of the three
            modes `Union' tends to be the most conservative and 
            often resuts in the lowest number of read counts. This is 
            because the method does not attempt to resolve a read that hits 
            multiple subjects but simply discards them. Both 
            `IntersectionStrict' and `IntersectionNotEmpty' have rules that
            attempt to assign a `multi-hit read' to a single feature.

      \item IntersectionStrict : The read must fall completely within a 
            single feature to be counted. A read can overlap multiple 
            features but must fall within only one. In the case of gapped reads, 
            all portions of the read fragment must fall within the same 
            feature for the read to be counted. The fragments can overlap 
            multiple features but collectively they must fall within 
            only one. 

      \item IntersectionNotEmpty : For this counting mode, the features are
            partitioned into unique disjoint regions. This is accomplished
            by disjoining the feature ranges then removing ranges shared by
            more than one feature. The result is a group of non-overlapping 
            regions each of which belong to a single feature. 
            Simple and gapped reads are counted if,

            \itemize{
            \item  the read or exactly 1 of the read fragments overlaps a unique
                   disjoint region
            \item the read or >1 read fragments overlap >1 unique disjoint
                  region from the same feature
            }
   }
  }
  \item{ignore.strand}{
    A logical value indicating if strand should be considered when matching.
  }
  \item{...}{Additional arguments for other methods.  If using multiple
    cores, you can pass arguments in here to be used by mclapply to
    indicate the number of cores to use etc.
  }
}

\details{
  \describe{
    \item{features :}{
      A `feature' can be any portion of a genomic region such as a gene, 
      transcript, exon etc. When the \code{features} argument is a 
      \linkS4class{GRanges} the rows define the features. The result
      will be the same length as the \linkS4class{GRanges}. When 
      \code{features} is a \linkS4class{GRangesList} the highest list-level 
      defines the features and the result will be the same length as the 
      \linkS4class{GRangesList}. 

      Each count `mode' attempts to assign a read that overlapps multiple 
      features to a single feature. If there are ranges that should be
      considered together (e.g., exons by transcript or cds regions by gene)
      the \linkS4class{GRangesList} would be appropriate. If there is no 
      grouping in the data then a \linkS4class{GRanges} would be appropriate. 
    }
    \item{paired-end reads :}{
      Paired-end reads should be provided in a 
      \linkS4class{GappedAlignmentPairs} container. Paired-end reads
      are counted the same as single-end reads with gaps. 

      The \link[Rsamtools]{BamFileList} and \link[Rsamtools]{BamViews} 
      methods have an additional argument, \code{singleEnd}, to indicate if 
      the bam files contain single or paired-end reads. 
      See ?\code{summarizeOverlaps,GRanges,BamFileList-method} for
      details.
    }
  }
}

\value{
  A \linkS4class{SummarizedExperiment} object. The \code{assays} slot holds 
  the counts, \code{rowData} holds the \code{features}, \code{colData}
  will either be \code{NULL} or hold any metadata that was present in the
  \code{reads}.
}

\references{
  HTSeq :
  \url{http://www-huber.embl.de/users/anders/HTSeq/doc/overview.html}

  htseq-count :
  \url{http://www-huber.embl.de/users/anders/HTSeq/doc/count.html}
}

\author{Valerie Obenchain <vobencha@fhcrc.org>}

\seealso{
  \itemize{
    \item \code{DESeq}, \code{DEXSeq} and \code{edgeR} packages 
    \item \link[Rsamtools]{BamFileList} and \link[Rsamtools]{BamViews} classes
    \item \link{GappedAlignments} and \link{GappedAlignmentPairs} classes
    \item \link{readGappedAlignments} and \link{readGappedAlignmentPairs} 
  }
}

\examples{
reads <- GappedAlignments(
    names = c("a","b","c","d","e","f","g"),
    seqnames = Rle(c(rep(c("chr1", "chr2"), 3), "chr1")),
    pos = as.integer(c(1400, 2700, 3400, 7100, 4000, 3100, 5200)),
    cigar = c("500M", "100M", "300M", "500M", "300M", 
        "50M200N50M", "50M150N50M"),
    strand = strand(rep("+", 7)))

## ---------------------------------------------------------------------
## 'features' as GRanges
##

features <- GRanges(
    seqnames = c(rep("chr1", 7), rep("chr2", 4)), strand = "+", 
    ranges = IRanges(c(1000, 3000, 3600, 4000, 4000, 5000, 5400, 2000, 
        3000, 7000, 7500), width = c(500, 500, 300, 500, 900, 500, 500, 
        900, 500, 600, 300)), 
    group_id = c("A", "B", "C", "C", "D", "D", "E", "F", "G", "H", "H"))


## Each row is a feature the read can 'hit'.
rowsAsFeatures <- data.frame(
    union = assays(summarizeOverlaps(features, reads))$counts, 
    intStrict = assays(summarizeOverlaps(features, reads, 
        mode="IntersectionStrict"))$counts,
    intNotEmpty = assays(summarizeOverlaps(features, reads,
        mode="IntersectionNotEmpty"))$counts,
    countOverlaps = countOverlaps(features, reads))

## Results from countOverlaps() are included to highlight how
## summarizeOverlaps() counts a read only once. For these 7
## reads, 'Union' is the most conservative counting mode, followed 
## by 'Intersectionstrict' and then 'IntersectionNotEmpty'. 

colSums(rowsAsFeatures)

## ---------------------------------------------------------------------
## 'features' as GRangesList
##

## Each highest list-level is a feature the read can 'hit'.
lst <- split(features, mcols(features)[["group_id"]])
listAsFeatures <- data.frame(
    union = assays(summarizeOverlaps(lst, reads))$counts, 
    intStrict = assays(summarizeOverlaps(lst, reads, 
        mode="IntersectionStrict"))$counts,
    intNotEmpty = assays(summarizeOverlaps(lst, reads,
        mode="IntersectionNotEmpty"))$counts,
    countOverlaps = countOverlaps(lst, reads))

## ---------------------------------------------------------------------
## 'reads' as BamFileList 
##

## Count BAM files and prepare output for DESeq or edgeR analysis
library(Rsamtools)
library(DESeq)
library(edgeR)

fls <- list.files(system.file("extdata",package="GenomicRanges"),
    recursive=TRUE, pattern="*bam$", full=TRUE)
names(fls) <- basename(fls)
bfl <- BamFileList(fls, index=character())

features <- GRanges(
    seqnames = c(rep("chr2L", 4), rep("chr2R", 5), rep("chr3L", 2)),
    ranges = IRanges(c(1000, 3000, 4000, 7000, 2000, 3000, 3600, 4000, 
        7500, 5000, 5400), width=c(rep(500, 3), 600, 900, 500, 300, 900, 
        300, 500, 500)), "-",
    group_id=c(rep("A", 4), rep("B", 5), rep("C", 2)))

solap <- summarizeOverlaps(features, bfl)

deseq <- newCountDataSet(assays(solap)$counts, rownames(colData(solap)))

edger <- DGEList(assays(solap)$counts, group=rownames(colData(solap)))

## ---------------------------------------------------------------------
## paired-end reads 
##

library("TxDb.Dmelanogaster.UCSC.dm3.ensGene")
exbygene <- exonsBy(TxDb.Dmelanogaster.UCSC.dm3.ensGene, "gene")
fl <- system.file("extdata", "untreated3_chr4.bam", 
    package="pasillaBamSubset") 

## Paired-end reads are stored in a GappedAlignmentPairs object
reads <- readGappedAlignmentPairs(fl)

res <- summarizeOverlaps(exbygene, reads, "Union")
stopifnot(length(assays(res)$counts) == length(exbygene))
}

\keyword{methods}
\keyword{utilities}

