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

  ## mode funtions
  Union(features, reads, ignore.strand=FALSE, inter.feature=TRUE)
  IntersectionStrict(features, reads, ignore.strand=FALSE, inter.feature=TRUE)
  IntersectionNotEmpty(features, reads, ignore.strand=FALSE, inter.feature=TRUE)
}

\arguments{

  \item{reads}{
    A \code{\link[Rsamtools]{BamFileList}} that represents the data to be
    counted by \code{summarizeOverlaps}.
  }
  \item{features}{
    A \link{GRanges} or a \link{GRangesList} object of genomic regions of
    interest. When a \link{GRanges} is supplied, each row is considered a
    feature. When a \link{GRangesList} is supplied, each higher list-level is
    considered a feature. This distinction is important when defining
    overlaps.
  } 
  \item{mode}{
    A function that defines the method to be used when a read overlaps
    more than one feature. Pre-defined options are "Union",
    "IntersectionStrict", or "IntersectionNotEmpty" and are designed
    after the counting modes available in the HTSeq package by Simon
    Anders (see references).

    \itemize{
      \item "Union" : (Default) Reads that overlap any portion of exactly one
            feature are counted. Reads that overlap multiple features are
            discarded. This is the most conservative of the 3 modes.
      \item "IntersectionStrict" : A read must fall completely "within" the
            feature to be counted. If a read overlaps multiple features but
            falls "within" only one, the read is counted for that feature.
            If the read is "within" multiple features, the read is discarded.
      \item "IntersectionNotEmpty" : A read must fall in a unique disjoint
            region of a feature to be counted. When a read overlaps multiple
            features, the features are partitioned into disjoint intervals.
            Regions that are shared between the features are discarded leaving
            only the unique disjoint regions. If the read overlaps one of
            these remaining regions, it is assigned to the feature the
            unique disjoint region came from.
      \item user supplied function : A function can be supplied as the
            \code{mode} argument. It must (1) have arguments that correspond
            to \code{features}, \code{reads}, \code{ignore.strand} and
            \code{inter.feature} arguments (as in the defined mode functions)
            and (2) return a vector of counts the same length as
            \code{features}. 
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
    behavior of modes \sQuote{Union} and \sQuote{IntersectionStrict} are 
    essentially \sQuote{countOverlaps} with \sQuote{type=any} and
    \code{type=within}, respectively. \sQuote{IntersectionNotEmpty} does
    not reduce to a simple countOverlaps because common (shared) regions 
    of the annotation are removed before counting.
  }
  \item{...}{
    Additional arguments such as \code{fragments}, \code{singleEnd} or 
    \code{param}. If using multiple cores, arguments can be passed through 
    to mclapply used when counting Bam files.
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
      \item{singleEnd}{
        A logical value indicating if reads are single or paired-end.
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
  the counts, \code{rowData} holds the annotation specified in 
  \code{features}. 

  \code{colData} is a DataFrame with columns of \sQuote{object} (class of 
  \code{reads}) and \sQuote{records} (length of \code{reads}). When \code{reads}
  is a BamFile or BamFileList the \code{colData} holds the output of a call
  to \code{countBam} with columns of \sQuote{records} (total records in file), 
  \sQuote{nucleotides} and \sQuote{mapped}. The number in \sQuote{mapped} is
  the number of records returned when \code{isUnmappedQuery=FALSE} in the 
  \sQuote{ScanBamParam}.
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

## When the reads are Bam files, the 'colData' contains summary 
## information from a call to countBam().
colData(se)

## Create count tables.
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

