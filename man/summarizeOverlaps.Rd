\name{summarizeOverlaps}

\alias{summarizeOverlaps}
\alias{summarizeOverlaps,GRanges,GAlignments-method}
\alias{summarizeOverlaps,GRangesList,GAlignments-method}
\alias{summarizeOverlaps,GRanges,GAlignmentsList-method}
\alias{summarizeOverlaps,GRangesList,GAlignmentsList-method}
\alias{summarizeOverlaps,GRanges,GAlignmentPairs-method}
\alias{summarizeOverlaps,GRangesList,GAlignmentPairs-method}
\alias{Union}
\alias{IntersectionStrict}
\alias{IntersectionNotEmpty}


\title{Perform overlap queries between reads and genomic features} 

\description{
  \code{summarizeOverlaps} extends \code{findOverlaps} by providing 
  options to resolve reads that overlap multiple features. 
}

\usage{
  \S4method{summarizeOverlaps}{GRanges,GAlignments}(
    features, reads, mode, ignore.strand=FALSE, ..., inter.feature=TRUE) 
  \S4method{summarizeOverlaps}{GRangesList,GAlignments}(
    features, reads, mode, ignore.strand=FALSE, ..., inter.feature=TRUE) 
  \S4method{summarizeOverlaps}{GRanges,GAlignmentPairs}(
    features, reads, mode, ignore.strand=FALSE, ..., inter.feature=TRUE) 
  \S4method{summarizeOverlaps}{GRangesList,GAlignmentPairs}(
    features, reads, mode, ignore.strand=FALSE, ..., inter.feature=TRUE) 
}

\arguments{
  \item{features}{
    A \linkS4class{GRanges} or a \linkS4class{GRangesList} containing genomic 
    regions of interest. This will commonly be a \code{GRanges} of exons or
    transcripts or a \linkS4class{GRangesList} of exons by gene or transcripts
    by gene.

  }
  \item{reads}{
    A \linkS4class{GAlignments} (single-end), 
    \linkS4class{GAlignmentPairs} (paired-end), or
    \link[Rsamtools]{BamFileList} or \link[Rsamtools]{BamViews} objects
    containing the reads to be mapped to the genomic regions of interest.

    The \link[Rsamtools]{BamFileList} and \link[Rsamtools]{BamViews} methods 
    can handle single or paired-end reads by appropriately specifying the
    \code{singleEnd} argument. Currently \code{singleEnd} is set once for 
    all files in the list (i.e., the list cannot contain both single and 
    paired-end read files.). For examples of iterating over Bam files 
    See ?\code{summarizeOverlaps,GRanges,BamFileList-method}.
  }
  \item{mode}{
    Character name of a function that defines the counting method to be used. 
    Modes include \sQuote{Union}, \sQuote{IntersectionStrict}, or 
    \sQuote{IntersectionNotEmpty} and are designed after the counting modes 
    in the HTSeq package by Simon Anders (see references). All methods are 
    gap-aware and count a read a maximum of once. Alternatively, a user can 
    provide their own count function as the \code{mode}. See details below.

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

            \enumerate{
              \item the read or exactly 1 of the read fragments overlaps a 
                    unique disjoint region
              \item the read or >1 read fragments overlap >1 unique disjoint
                    region from the same feature
            }

      \item User supplied function : The function must meet the following
            requirements:

            \enumerate{
              \item have 3 arguments that correspond to `features', `reads' and 
                    `ignore.strand', in that order
              \item return a vector of counts the same length as `features'
            }
    } 
  }
  \item{ignore.strand}{
    A logical indicating if strand should be considered when matching.
  }
  \item{inter.feature}{
    A logical indicating if the counting \code{mode} should be aware of
    overlapping features. When TRUE (default), reads mapping to multiple 
    features are dropped (i.e., not counted). When FALSE, these reads are 
    retained and a count is assigned to each feature they map to.

    There are 6 possible combinations of the \code{mode} and 
    \code{inter.feature} arguments. When \code{inter.feature=FALSE} the
    behavior of modes \sQuote{Union} and \sQuote{IntersectionNotEmpty} are the
    same resulting in 5 distinct ways to count.
  }
  \item{...}{
    Additional arguments such as \code{fragments} or \code{param}.  If 
    using multiple cores, arguments can be passed through to the mclapply 
    used when counting Bam files.
    \describe{
      \item{fragments}{
      A logical value indicating if singletons, reads with unmapped pairs
      and other fragments should be included in the counting. When
      \code{fragments=FALSE} only reads paired with the algorithm described
      at ?findMateAlignment are counted. When \code{fragments=TRUE} (default)
      all singletons, reads with unmapped pairs and other fragments are
      counted in addition to the reads paired with the ?findMateAlignment
      algorithm. This argument applies to paired-end reads only so
      \code{singleEnd} must be FALSE.
      }
      \item{param}{An optional \code{\link[Rsamtools]{ScanBamParam}} instance to
         further influence scanning, counting, or filtering.
      }
    }
  }
}

\details{
  \describe{
    \item{}{\code{summarizeOverlaps} offers counting modes to resolve reads 
      that overlap multiple features. The \code{mode} argument defines a
      set of rules to resolve the read to a single feature such that each read 
      is counted a maximum of once. New to GenomicRanges >= 1.13.9 is the
      \code{inter.feature} argument which allows reads to be counted for
      each feature they overlap. When \code{inter.feature=TRUE} the counting
      modes are aware of feature overlap and reads overlapping multiple 
      features are dropped and not counted. (This was the standard behavior 
      prior to GenomicRanges 1.13.9.) When \code{inter.feature=FALSE} multiple
      feature overlap is ignored and reads are counted once for each feature
      they map to. This essentially reduces modes \sQuote{Union} and
      \sQuote{IntersectionStrict} to \code{countOverlaps} with 
      \code{type="any"}, and \code{type="within"}, respectively.
      \sQuote{IntersectionNotEmpty} is not reduced to a derivative of
      \code{countOverlaps} because the shared regions are removed before 
      counting.
    }
    \item{features :}{
      A \sQuote{feature} can be any portion of a genomic region such as a gene, 
      transcript, exon etc. When the \code{features} argument is a 
      \linkS4class{GRanges} the rows define the features. The result
      will be the same length as the \linkS4class{GRanges}. When 
      \code{features} is a \linkS4class{GRangesList} the highest list-level 
      defines the features and the result will be the same length as the 
      \linkS4class{GRangesList}. 

      When \code{inter.feature=TRUE}, each count \code{mode} attempts to 
      assign a read that overlaps multiple features to a single feature. If 
      there are ranges that should be considered together (e.g., exons by 
      transcript or cds regions by gene) the \linkS4class{GRangesList} would 
      be appropriate. If there is no grouping in the data then a 
      \linkS4class{GRanges} would be appropriate. 
    }
    \item{paired-end reads :}{
      Paired-end reads should be provided in a 
      \linkS4class{GAlignmentPairs} container. Paired-end reads
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
    \item \link{GAlignments} and \link{GAlignmentPairs} classes
    \item \link{readGAlignments} and \link{readGAlignmentPairs} 
  }
}

\examples{
reads <- GAlignments(
    names = c("a","b","c","d","e","f","g"),
    seqnames = Rle(c(rep(c("chr1", "chr2"), 3), "chr1")),
    pos = as.integer(c(1400, 2700, 3400, 7100, 4000, 3100, 5200)),
    cigar = c("500M", "100M", "300M", "500M", "300M", 
              "50M200N50M", "50M150N50M"),
    strand = strand(rep("+", 7)))

gr <- GRanges(
    seqnames = c(rep("chr1", 7), rep("chr2", 4)), strand = "+", 
    ranges = IRanges(c(1000, 3000, 3600, 4000, 4000, 5000, 5400, 
                       2000, 3000, 7000, 7500), 
                     width = c(500, 500, 300, 500, 900, 500, 500, 
                               900, 500, 600, 300),
                     names=c("A", "B", "C1", "C2", "D1", "D2", "E", "F",
                             "G", "H1", "H2"))) 
groups <- factor(c(1,2,3,3,4,4,5,6,7,8,8))
grl <- splitAsList(gr, groups)
names(grl) <- LETTERS[seq_along(grl)]

## ---------------------------------------------------------------------
## Counting modes. 
## ---------------------------------------------------------------------

## First we count with a GRanges as the 'features'. Note that
## 'Union' is the most conservative counting mode followed by 
## 'IntersectionStrict' then 'IntersectionNotEmpty'.
counts1 <- 
    data.frame(union=assays(summarizeOverlaps(gr, reads))$counts, 
               intStrict=assays(summarizeOverlaps(gr, reads, 
                                mode="IntersectionStrict"))$counts,
               intNotEmpty=assays(summarizeOverlaps(gr, reads,
                                  mode="IntersectionNotEmpty"))$counts)

colSums(counts1)

## Split the 'features' into a GRangesList and count again.
counts2 <- 
    data.frame(union=assays(summarizeOverlaps(grl, reads))$counts, 
               intStrict=assays(summarizeOverlaps(grl, reads, 
                                mode="IntersectionStrict"))$counts,
               intNotEmpty=assays(summarizeOverlaps(grl, reads,
                                  mode="IntersectionNotEmpty"))$counts)
colSums(counts2)

## The GRangesList ('grl' object) has 8 features whereas the GRanges 
## ('gr' object) has 11. The affect on counting can be seen by looking
## at feature 'H' with mode 'Union'. In the GRanges this feature is 
## represented by ranges 'H1' and 'H2',
gr[c("H1", "H2")]

## and by list element 'H' in the GRangesList, 
grl["H"]
 
## Read "d" hits both 'H1' and 'H2'. This is considered a multi-hit when
## using a GRanges (each range is a separate feature) so the read was 
## dropped and not counted.
counts1[c("H1", "H2"), ]

## When using a GRangesList, each list element is considered a feature.
## The read hits multiple ranges within list element 'H' but only one 
## list element. This is not considered a multi-hit so the read is counted.
counts2["H", ]

## ---------------------------------------------------------------------
## Counting multi-hit reads.
## ---------------------------------------------------------------------

## The goal of the counting modes is to provide a set of rules that
## resolve reads hitting multiple features so each read is counted
## a maximum of once. However, sometimes it may be desirable to count 
## a read for each feature it overlaps. This can be accomplished by 
## setting 'inter.feature' to FALSE.

## When 'inter.feature=FALSE', modes 'Union' and 'IntersectionStrict'
## essentially reduce to countOverlaps() with type="any" and 
## type="within", respectively.

## When 'inter.feature=TRUE' only features "A", "F" and "G" have counts.
se1 <- summarizeOverlaps(gr, reads, mode="Union", inter.feature=TRUE)
assays(se1)$counts

## When 'inter.feature=FALSE' all 11 features have a count. There are 
## 7 total reads so one or more reads were counted more than once.
se2 <- summarizeOverlaps(gr, reads, mode="Union", inter.feature=FALSE)
assays(se2)$counts

## ---------------------------------------------------------------------
## Counting Bam files.
## ---------------------------------------------------------------------

library(Rsamtools)
library(pasillaBamSubset)
library("TxDb.Dmelanogaster.UCSC.dm3.ensGene")
exbygene <- exonsBy(TxDb.Dmelanogaster.UCSC.dm3.ensGene, "gene")

## (i) Single-end reads:

## 'yieldSize' can be set to iterate over large files in chunks.
## When counting, 'singleEnd' should be TRUE (default).
bf_s <- BamFile(untreated1_chr4(), yieldSize=50000)
se_s <- summarizeOverlaps(exbygene, bf_s, singleEnd=TRUE)
table(assays(se_s)$counts > 0)

## (ii) Paired-end reads:

## A paired-end file may contain singletons, reads with unmapped
## pairs or reads with more than two fragments. When 'fragments=FALSE'
## only the reads paired by the ?findMateAlignment algorithm will
## be included in the counting. In this case, 'yieldSize' cannot be used.
bf_p <- BamFile(untreated3_chr4())
se_p <- summarizeOverlaps(exbygene, bf_p, singleEnd=FALSE, fragments=FALSE)
table(assays(se_p)$counts > 0)

## When 'fragments=TRUE' (default) all singletons, reads with
## unmapped pairs and other fragments will be included in the
## counting as well as the reads paired with ?findMateAlignment.
## When 'fragments=TRUE', 'yieldSize' can be used and the file
## must be sorted by qname.
fl <- sortBam(untreated3_chr4(), tempfile(), byQname=TRUE)
bf_psort <- BamFile(fl, index=character(), yieldSize=90000, obeyQname=TRUE)
se_psort <- summarizeOverlaps(exbygene, bf_psort, singleEnd=FALSE)
table(assays(se_psort)$counts > 0)

## Both the qname-sorted and unsorted examples used the same Bam file,
## annotation, and counting method. As expected, using 'fragments=TRUE' 
## results in a larger number of total counts because singletons,
## unmapped pairs etc. are included in the counting.

## Total reads in the file:
countBam(untreated3_chr4())

## Reads counted with 'fragments=TRUE':
sum(assays(se_psort)$counts)

## Reads counted with 'fragments=TRUE':
sum(assays(se_p)$counts)

## ---------------------------------------------------------------------
## Count tables for DESeq or edgeR.
## ---------------------------------------------------------------------

fls <- list.files(system.file("extdata",package="GenomicRanges"),
                  recursive=TRUE, pattern="*bam$", full=TRUE)
names(fls) <- basename(fls)
bf <- BamFileList(fls, index=character(), yieldSize=1000)
genes <- GRanges(
    seqnames = c(rep("chr2L", 4), rep("chr2R", 5), rep("chr3L", 2)),
    ranges = IRanges(c(1000, 3000, 4000, 7000, 2000, 3000, 3600, 
                       4000, 7500, 5000, 5400), 
                     width=c(rep(500, 3), 600, 900, 500, 300, 900, 
                             300, 500, 500))) 

se <- summarizeOverlaps(genes, bf)

library(DESeq)
deseq <- newCountDataSet(assays(se)$counts, rownames(colData(se)))
library(edgeR)
edger <- DGEList(assays(se)$counts, group=rownames(colData(se)))

## ---------------------------------------------------------------------
## User supplied 'mode'. 
## ---------------------------------------------------------------------

## A user defined count function must have the same arguments as 
## the current counting modes.
\dontrun{
counter <- function(x, y,  ignore.strand, inter.feature) {
   ## count ...
}

se <- summarizeOverlaps(gr, reads, mode=counter) 
}
}

\keyword{methods}
\keyword{utilities}

