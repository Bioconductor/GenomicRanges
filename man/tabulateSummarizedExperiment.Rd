\name{tabulateBamFiles}
\docType{methods}
\alias{tabulateBamFiles}
\alias{tabulateBamFiles,character-method}
\alias{tabulateBamFiles,BamViews-method}
\alias{tabulateBamFiles,BamFileList-method}

\title{
  tabulateBamFiles helper methods
}

\description{
  Create a \link{SummarizedExperiment} object from a list of .bam files
  or from an object that knows where .bam files are located.
}

\usage{
  \S4method{tabulateBamFiles}{character}(x,annot,columnData,annotType,strandAgnostic) 
  \S4method{tabulateBamFiles}{BamViews}(x,annot,columnData,annotType,strandAgnostic) 
  \S4method{tabulateBamFiles}{BamFileList}(x,annot,columnData,annotType,strandAgnostic) 
}

\arguments{
  \item{x}{
    A character vector to files, a \code{BamFileList} or a
    \code{BamViews} object. If the character vector is named, the names
    will be used as sample labels.
  }
  \item{annot}{
    An annotation object in the form of a GRanges object.  Counts will
    be based on the individual elements of this element such that the
    number of rows in the counts table will be equal to the length of the
    GRanges object.
  }
  \item{columnData}{
  A \code{data.frame} that contains columns describing the details of
  the experiment. The rownames of this \code{data.frame} object must match the
  names of x (or the values of x for the case where the names of x are null).
  }
  \item{annotType}{
  This specifies the id type to use for uniquely labelling the rows of
  the tabulated data.  Acceptable values can be "tx_id", "exon_id"
  "gene_id" or  NULL.  The idea is that this is the id in your annot
  object that you are using to count with, so if you a \code{GRanges}
  object that has an element for each exon, then you will want to use
  "exon_id" etc.
  }
  \item{strandAgnostic}{TRUE of FALSE value used to decide whether to
  ignore strand information that is present in the .bam files when counting.
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

## example with GRangesList
txs <- transcriptsBy(txdb, by="gene")
sx2 <- tabulateBamFiles(fls, annot=txs, columnData = expData,
                       annotType="exon_id", strandAgnostic=TRUE)


}
}


\author{Marc Carlson}

\keyword{methods}
\keyword{utilities}

