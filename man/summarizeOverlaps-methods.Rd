\name{summarizeOverlaps}

\alias{summarizeOverlaps}
\alias{summarizeOverlaps,GRanges,GappedAlignments-method}
\alias{summarizeOverlaps,GRangesList,GappedAlignments-method}
\alias{Union}
\alias{IntersectionStrict}
\alias{IntersectionNotEmpty}


\title{Count reads that map to genomic features} 

\description{
  Count reads that map to genomic features with options to resolve reads that
  overlap multiple features 
}

\usage{
  \S4method{summarizeOverlaps}{GRanges,GappedAlignments}(
    features, reads, mode, ignore.strand = FALSE, ..., param = ScanBamParam()) 
}

\arguments{
  \item{features}{
    A \linkS4class{GRanges} or a \linkS4class{GRangesList} object of genomic 
    regions of interest. When a \linkS4class{GRanges} is supplied, each row is 
    considered a different feature. When a \linkS4class{GRangesList} is 
    supplied, each highest list-level is considered a feature and the multiple 
    elements are considered portions of the same feature. See examples for 
    details. 
  }
  \item{reads}{
    A \linkS4class{GappedAlignments}, \link[Rsamtools]{BamFileList} or a
    \link[Rsamtools]{BamViews} object. 
  }
  \item{mode}{
    A character names of function that defines the counting method to be used. 
    Pre-defined methods are "Union", "IntersectionStrict", or 
    "IntersectionNotEmpty" and are designed after the counting modes available 
    in the HTSeq package by Simon Anders (see references).

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
  \item{...}{Additional arguments for other methods.  If using multiple
    cores, you can pass arguments in here to be used by mclapply to
    indicate the number of cores to use etc.
  }
}

\details{
  In the context of \code{summarizeOverlaps} a "feature" can be any portion of a 
  genomic region such as a gene, transcript, exon etc. When the \code{features} 
  argument is a \linkS4class{GRanges} the rows define the features to be 
  overlapped. When \code{features} is a \linkS4class{GRangesList} the highest 
  list-levels define the features. 

  \code{summarizeOverlaps} offers three \code{mode} functions to handle reads
  that overlap multiple features: "Union", "IntersectionStrict", and
  "IntersectionNotEmpty".  These functions are patterned after the counting
  methods in the HTSeq package (see references). Each mode has a set of rules
  that dictate how a read is assigned. Reads are counted a maximum of once. 
  Alternatively, users can provide their own counting function as the 
  \code{mode} argument and take advantage of the infrastructure in 
  \code{summarizeOverlaps} to count across multiple files and parse the results 
  into a \linkS4class{SummarizedExperiment} object.

  Currently reads must be input as either a BAM file or a 
  \linkS4class{GappedAlignments} object. The information in the CIGAR field 
  is used to determine if gapped reads are present. Gapped reads are counted 
  with the same \code{mode} as the simple reads.

  NOTE : \code{summarizeOverlaps} does not currently handle paired-end reads. 
}

\value{
  A \linkS4class{SummarizedExperiment} object. The \code{assays} slot holds the
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
  \link[Rsamtools]{BamFileList} 
  \link[Rsamtools]{BamViews} 
}

\examples{

  group_id <- c("A", "B", "C", "C", "D", "D", "E", "F", "G", "H", "H")
  features <- GRanges(
      seqnames = Rle(c("chr1", "chr2", "chr1", "chr1", "chr2", "chr2", 
          "chr1", "chr1", "chr2", "chr1", "chr1")),
      strand = strand(rep("+", length(group_id))),
      ranges = IRanges(
          start=c(1000, 2000, 3000, 3600, 7000, 7500, 4000, 4000, 3000, 
              5000, 5400),
          width=c(500, 900, 500, 300, 600, 300, 500, 900, 500, 500, 500)),
     DataFrame(group_id)
  )
 
  reads <- GappedAlignments(
      names = c("a","b","c","d","e","f","g"),
      rname = Rle(c(rep(c("chr1", "chr2"), 3), "chr1")),
      pos = as.integer(c(1400, 2700, 3400, 7100, 4000, 3100, 5200)),
      cigar = c("500M", "100M", "300M", "500M", "300M", 
          "50M200N50M", "50M150N50M"),
      strand = strand(rep.int("+", 7L)))

  ## Results from countOverlaps are included to highlight how the 
  ## modes in summarizeOverlaps count a read a maximum of once.

  ## When the 'features' argument is a GRanges, each row 
  ## is treated as a different feature. 
  rowsAsFeatures <- 
      data.frame(union = assays(summarizeOverlaps(features, reads))$counts, 
                 intStrict = assays(summarizeOverlaps(features, reads, 
                     mode="IntersectionStrict"))$counts,
                 intNotEmpty = assays(summarizeOverlaps(features, reads,
                     mode="IntersectionNotEmpty"))$counts,
                 countOverlaps = countOverlaps(features, reads))

  ## When the 'features' argument is a GRangesList, each
  ## highest list-level is a different feature.
  lst <- split(features, values(features)[["group_id"]])
  listAsFeatures <- 
      data.frame(union = assays(summarizeOverlaps(lst, reads))$counts, 
                 intStrict = assays(summarizeOverlaps(lst, reads, 
                     mode="IntersectionStrict"))$counts,
                 intNotEmpty = assays(summarizeOverlaps(lst, reads,
                     mode="IntersectionNotEmpty"))$counts,
                 countOverlaps = countOverlaps(lst, reads))

  ## Read across BAM files and package output for DESeq or edgeR analysis
  library(Rsamtools)
  library(DESeq)
  library(edgeR)

  fls = list.files(system.file("extdata",package="GenomicRanges"),
      recursive=TRUE, pattern="*bam$", full=TRUE)
  bfl <- BamFileList(fls)
  features <- GRanges(
      seqnames = Rle(c("chr2L", "chr2R", "chr2L", "chr2R", "chr2L", "chr2R",
          "chr2L", "chr2R", "chr2R", "chr3L", "chr3L")),
      strand = strand(rep("+", 11)),
      ranges = IRanges(start=c(1000, 2000, 3000, 3600, 7000, 7500, 4000, 4000,
          3000, 5000, 5400), width=c(500, 900, 500, 300, 600, 300, 500, 900, 
          500, 500, 500))
  )

  solap <- summarizeOverlaps(features, bfl)

  deseq <- newCountDataSet(countData=assays(solap)$counts, 
                           conditions=rownames(colData(solap)))

  edger <- DGEList(counts=assays(solap)$counts, group=rownames(colData(solap)))
}

\keyword{methods}
\keyword{utilities}

