\name{countFeatureHits}

\alias{countFeatureHits}
\alias{countFeatureHits,character,GRanges-method}
\alias{countFeatureHits,BamFileList,GRanges-method}
\alias{countFeatureHits,BamViews,GRanges-method}
\alias{countFeatureHits,character,GRangesList-method}
\alias{countFeatureHits,BamFileList,GRangesList-method}
\alias{countFeatureHits,BamViews,GRangesList-method}
\alias{countFeatureHits,GappedAlignments,GRanges-method}
\alias{countFeatureHits,GappedAlignments,GRangesList-method}


\title{Count reads that map to genomic features} 

\description{
  Count reads that map to genomic features with options to resolve reads that
  overlap multiple features 
}

\usage{
  \S4method{countFeatureHits}{GappedAlignments,GRanges}(
    reads, features, mode = "Union",
    ignore.strand = FALSE, ..., param = ScanBamParam()) 
}

\arguments{
  \item{reads}{
    A character() vector with the name(s) of the BAM file(s) or a 
    \link{GappedAlignments} object. 
  }
  \item{features}{
    A \link{GRanges} or a \link{GRangesList} object of genomic regions of 
    interest. When a \link{GRanges} is supplied, each row is considered a 
    feature. When a \link{GRangesList} is supplied, each higher list-level is
    considered a feature. This distinction is important when defining an overlap
    between a read and a feature. See examples for details. 
  }
  \item{mode}{
    A \code{character(1)} string indicating what resolution method should be
    used when a read overlaps more than one feature. Options are "Union", 
    "IntersectionStrict", or "IntersectionNotEmpty" and are designed after
    the counting modes available in the HTSeq package by Simon Anders 
    (see references).

    \itemize{
      \item "Union" : (Default) Reads that overlap any portion of exactly one 
            feature are counted. Reads that overlap multiple features are 
            discarded. 
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
   }
  }
  \item{param}{An optional \link[Rsamtools]{ScanBamParam} instance to
     further influence scanning, counting, or filtering of the BAM file.}

  \item{ignore.strand}{
    A logical value indicating if strand should be considered when matching.
  }
  \item{...}{Additional arguments for other methods
  }
}

\details{
  In the context of \code{countFeatureHits} a "feature" can be any portion of a 
  genomic region such as a gene, transcript, exon etc. When the \code{features} 
  argument is a \link{GRanges} the rows define the features to be overlapped.
  When \code{features} is a \link{GRangesList} the highest list-levels define
  the features. 

  \code{countFeatureHits} offers three modes to handle reads that overlap
  multiple features: "Union", "IntersectionStrict", and "IntersectionNotEmpty".
  These modes are patterned after the counting methods in the HTSeq package
  (see references). Each mode follows rules that dictate how the read hit should 
  be assigned. Reads are counted a maximum of once.

  Currently reads must be input as either a BAM file or a
  \link{GappedAlignments} object. The information in the CIGAR field is used
  to determine if gapped reads are present. Gapped reads are counted with
  the same \code{mode} as the simple reads.

  NOTE : \code{countFeatureHits} does not currently handle paired-end reads. 
}

\value{
    A \link{SummarizedExperiment} object. The \code{assays} slot holds the
    counts, \code{rowData} holds the \code{features}, \code{colData}
    will either be \code{NULL} or hold any metadata that was present in the
    \code{reads}.
}


\references{
  \url{http://www-huber.embl.de/users/anders/HTSeq/doc/overview.html} 
  home page for HTSeq

 \url{http://www-huber.embl.de/users/anders/HTSeq/doc/count.html}
  counting with htseq-count
}

\author{Valerie Obenchain <vobencha@fhcrc.org>}

\seealso{
  \code{DESeq} and \code{edgeR} packages
}

\examples{

group_id <- c("A", "B", "C", "C", "D", "D", "E", "F", "G", "G", "H", "H")
features <- GRanges(
    seqnames = Rle(c("chr1", "chr2", "chr1", "chr1", "chr2", "chr2", 
        "chr1", "chr1", "chr2", "chr2", "chr1", "chr1")),
    strand = strand(rep("+", length(group_id))),
    ranges = IRanges(
        start=c(1000, 2000, 3000, 3600, 7000, 7500, 4000, 4000, 3000, 3350, 5000, 5400),
        width=c(500, 900, 500, 300, 600, 300, 500, 900, 150, 200, 500, 500)),
   DataFrame(group_id)
)
 
reads <- GappedAlignments(
    names = c("a","b","c","d","e","f","g"),
    rname = Rle(c(rep(c("chr1", "chr2"), 3), "chr1")),
    pos = as.integer(c(1400, 2700, 3400, 7100, 4000, 3100, 5200)),
    cigar = c("500M", "100M", "300M", "500M", "300M", "50M200N50M", "50M150N50M"),
    strand = strand(rep.int("+", 7L)))

## Results from countOverlaps are included in these
## examples to highlight how the decision logic in 
## countFeatureHits allows a read to be counted a maximum of once.

## When the 'features' argument is a GRanges, each row 
## is treated as feature. 
rowsAsFeatures <- 
    data.frame(union = assays(countFeatureHits(reads, features))$counts, 
               intStrict = assays(countFeatureHits(reads, features, 
                   mode="IntersectionStrict"))$counts,
               intNotEmpty = assays(countFeatureHits(reads, features,
                   mode="IntersectionNotEmpty"))$counts,
               countOverlaps = countOverlaps(features, reads))

## When the 'features' argument is a GRangesList, each
## highest list-level is treated as a feature.
lst <- split(features, values(features)[["group_id"]])
listAsFeatures <- 
    data.frame(union = assays(countFeatureHits(reads, lst))$counts, 
               intStrict = assays(countFeatureHits(reads, lst, 
                   mode="IntersectionStrict"))$counts,
               intNotEmpty = assays(countFeatureHits(reads, lst,
                   mode="IntersectionNotEmpty"))$counts,
               countOverlaps = countOverlaps(lst, reads))

## FIXME : find sample data
## Use params to subset reads, package output for input to DESeq or edgeR
#seqs <- scanBamHeader(bamFile)[[1]]$targets
#param <- ScanBamParam(which=GRanges(names(seqs)[1], IRanges(1, seqs[1])))
#counts <- countFeatureHits(bamFile, features, param=param)
#
#library(DESeq)
#deseq <- newCountDataSet(countData=assays(counts)$counts, 
#                         conditions= ,
#                         phenoData= ,
#                         featureData=  )
#
#library(edgeR)
#edger <- DGEList(counts=assays(counts)$counts, 
#                 group= )
#
}

\keyword{methods}
\keyword{utilities}

