\name{utils}
\docType{methods}
\alias{tabulateSummarizedExperiment}
\alias{tabulateSummarizedExperiment,GenomicRanges,character-method}
\alias{tabulateSummarizedExperiment,GenomicRanges,BamViews-method}
\alias{tabulateSummarizedExperiment,GenomicRanges,BamFileList-method}

\title{
  tabulateSummarizedExperiment helper methods
}

\description{
  Create a \link{SummarizedExperiment} object from a list of .bam files
  or from an object that knows where .bam files are located.
}

\usage{
  \S4method{tabulateSummarizedExperiment}{GenomicRanges,character}(x, ...) 
  \S4method{tabulateSummarizedExperiment}{GenomicRanges,BamViews}(x, ...) 
  \S4method{tabulateSummarizedExperiment}{GenomicRanges,BamFileList}(x, ...) 
}

\arguments{
  \item{x}{
    A character vector to files, a \code{BamFileList} or a
    \code{BamViews} object. If the character vector is named, the names
    will be used as sample labels.
  }

}

\value{
  A \code{SummarizedExperiment} object with aligned reads tabulated
  according to the provided annotation object.
}

\examples{
\dontrun{
## Requirement #1 is that we want the user to provide a GRanges object
library(TxDb.Hsapiens.UCSC.hg18.knownGene)
library(GenomicFeatures)
library(ShortRead)
txdb <- Hsapiens_UCSC_hg18_knownGene_TxDb
## and they can do whatever filtering is needed to produce the annot
keep <- c(paste("chr", c(1:22), sep=""), "chrX", "chrY")
isActiveSeq(txdb)[seqlevels(txdb)] <- FALSE
isActiveSeq(txdb)[keep] <- TRUE
exons <- exons(txdb, columns=c("exon_id","gene_id","tx_id","tx_name"))
names(exons) <- values(exons)$exon_id


## AND requirement #2 is a dir is needed that points to some data
dir <- "/mnt/cpl/data/Solexa/SOC/101101/Samples"

## We will wrap this with a method where dir is character later, but for now...
## indexed BAM files
fls <- sub(".bai$", "", list.files(dir, recursive=TRUE,
                                   pattern="*bai$", full=TRUE))
##TEMP restrict so it doesn't take too long to test:
fls <- fls[1:2]
fileNames = sub("tophat_", "", basename(dirname(fls)))
names(fls) = fileNames

## fls: named character vector to specify the files (and sample names = names)
## annot: GRanges Object to specify the annotations.
## strandAgnostic: Ignore the strand of the data sequence?
## annotType: type of annotation for the rows. (exon_id, tx_id, cds_id, null)

## Requirement #3 is a data.frame of the colData where the rownames are the
## sample names.
expData <- data.frame(condition=sub("\\_.+","",fileNames, perl=TRUE))
rownames(expData) <- fileNames


## the core function will take a bam file list 
sx <- tabulateBamFiles(fls, annot=exons, columnData = expData,
                       annotType="exon_id", strandAgnostic=TRUE)

}
}


\author{Marc Carlson}

\keyword{methods}
\keyword{utilities}

