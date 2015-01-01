\name{coordinate-mapping-methods}
\alias{coordinate-mapping-methods}
\alias{coordinate-mapping}

\alias{mapToGenome}
\alias{mapToGenome,GenomicRanges,GenomicRanges-method}
\alias{mapToGenome,GenomicRanges,GRangesList-method}
\alias{pmapToGenome}
\alias{pmapToGenome,Ranges,GenomicRanges-method}
\alias{pmapToGenome,GenomicRanges,GenomicRanges-method}
\alias{pmapToGenome,GenomicRanges,GRangesList-method}

\alias{mapToTranscript}
\alias{mapToTranscript,GenomicRanges,GenomicRanges-method}
\alias{mapToTranscript,GenomicRanges,GRangesList-method}
\alias{pmapToTranscript}
\alias{pmapToTranscript,GenomicRanges,GenomicRanges-method}
\alias{pmapToTranscript,GenomicRanges,GRangesList-method}


\title{Map range coordinates between genome and transcript space}

\description{
  Coordinate mapping methods for \linkS4class{GenomicRanges} objects.

  See \code{?`\link[GenomicAlignments]{coordinate-mapping-methods}`} in the 
  \pkg{GenomicAlignments} package for methods on 
  \link[GenomicAlignment]{GAlignment} objects.
}

\usage{

## mapToGenome, pmapToGenome
\S4method{mapToGenome}{GenomicRanges,GenomicRanges}(x, alignment, 
          ignore.strand = TRUE, ...) 
\S4method{mapToGenome}{GenomicRanges,GRangesList}(x, alignment, 
          ignore.strand = TRUE, ...) 

\S4method{pmapToGenome}{Ranges,GenomicRanges}(x, alignment, ...)
\S4method{pmapToGenome}{GenomicRanges,GenomicRanges}(x, alignment, 
          ignore.strand = TRUE, ...) 
\S4method{pmapToGenome}{GenomicRanges,GRangesList}(x, alignment, 
          ignore.strand = TRUE, ...) 

## mapToTranscript, pmapToTranscript
\S4method{mapToTranscript}{GenomicRanges,GenomicRanges}(x, alignment, 
          ignore.strand = TRUE, ...) 
\S4method{mapToTranscript}{GenomicRanges,GRangesList}(x, alignment, 
          ignore.strand = TRUE, ...) 
\S4method{pmapToTranscript}{GenomicRanges,GenomicRanges}(x, alignment, 
          ignore.strand = TRUE, ...) 
\S4method{pmapToTranscript}{GenomicRanges,GRangesList}(x, alignment, 
          ignore.strand = TRUE, ...) 
}

\arguments{
  \item{x}{
    \code{\linkS4class{GenomicRanges}} object of positions to be mapped.

    \itemize{
      \item{} For \code{mapToTranscript} \code{x} is in genome-based 
        coordinates 
      \item{} For \code{mapToGenome} \code{x} is in transcript-based 
        coordinates 
    }
  }
  \item{alignment}{
    A \code{\linkS4class{GenomicRanges}} object that represents the alignment 
    (aka mapping) between the sequences in \code{x} and the sequences in the
    result.

    The current implementation assumes that positions in \code{x} 
    occur within the \code{alignment} region.
  }
  \item{ignore.strand}{
    When TRUE, strand is ignored in overlap operations. When \code{alignment} 
    is a \code{GRangesList} and \code{ignore.strand=FALSE} all inner list 
    elements of a common outer list element have the same strand. 
  }
  \item{\dots}{
    Arguments passed to other methods.
  }
}

\details{

  In this documentation the term `parallel` is used to describe the mapping
  of i-th element with i-th element and is not related to code execution
  with multiple cores or threads.

  \itemize{
    \item{Mapping strategy:}{
      \itemize{
        \item \code{mapToTranscript} attempts to map all elements of \code{x} to
          all elements of \code{alignment}. The rational is that \code{x} is in
          genomic (global) coordinate and could potentially map to any of the
          transcript-centric ranges in \code{alignment}.

        \item \code{mapToGenome} uses name matching between \code{x} and
          \code{alignment} to determine which pairs should be mapped. In this 
          case \code{x} is in transcript-based coordinates and is only mapped
          to ranges in \code{alignment} with the same name (i.e., from the same 
          gene or transcript).

        \item Both `parallel` methods, \code{pmapToTranscript} and
          \code{pmapToGenome} map only the i-th element of \code{x} with the 
          i-th element of \code{alignment}.
      }
    }
    \item{Unmapped ranges:}{
      In the element-wise (aka "parallel") methods, ranges in \code{x} that 
      do not map to \code{alignment} are returned as zero-width ranges
      starting at 1 with a seqname of "unmapped". In the non-parallel
      methods unmapped ranges are not returned.
    }
  }
}

\value{
  An object the same class as \code{x}.

  Parallel methods return an object the same shape as \code{x}. Ranges that
  are non-hits (out of bounds) or a strand-mismatch are returned as zero-width
  ranges starting at 1.

  Non-parallel methods return an object that varies in length similar to a
  Hits object. The result only contains mapped records, strand mismatch
  and non-hits are not returned.
}

\seealso{
  \itemize{
    \item \code{?`coordinate-mapping-methods`} in the \pkg{GenomicAlignments} 
      package for methods on \link[GenomicAlignments]{GAlignments} objects.
  }
}

\examples{

## ---------------------------------------------------------------------
## A. Map local sequence locations to the genome
## ---------------------------------------------------------------------

## NAGNAG alternative splicing plays an essential role in biological processes 
## and represents a highly adaptable system for posttranslational regulation 
## of gene function. The majority of NAGNAG studies largely focus on messenger 
## RNA. A study by Sun, Lin, and Yan 
## (http://www.hindawi.com/journals/bmri/2014/736798/) demonstrated that
## NAGNAG splicing is also operative in large intergenic noncoding RNA
## (lincRNA). 

## One finding of interest was that linc-POLR3G-10 exhibited two NAGNAG 
## acceptors located in two distinct transcripts: TCONS_00010012 and 
## TCONS_00010010. 

## Extract the exon coordinates of TCONS_00010012 and TCONS_00010010: 
lincrna <- c("TCONS_00010012", "TCONS_00010010")
library(TxDb.Hsapiens.UCSC.hg19.lincRNAsTranscripts)
txdb <- TxDb.Hsapiens.UCSC.hg19.lincRNAsTranscripts
exons <- exonsBy(txdb, by="tx", use.names=TRUE)[lincrna]
exons

## The two NAGNAG acceptors were identified in the upstream region of 
## the fourth and fifth exons located in TCONS_00010012.
## Extract the sequences for transcript TCONS_00010012:
library(BSgenome.Hsapiens.UCSC.hg19)
genome <- BSgenome.Hsapiens.UCSC.hg19
exons_seq <- getSeq(genome, exons[[1]])

## TCONS_00010012 has 4 exons:
exons_seq

## The most common triplet among the lincRNA sequences was CAG. Identify
## the location of this pattern in all exons.
cag_loc <- vmatchPattern("CAG", exons_seq)

## Convert the first occurance of CAG in each exon back to genome coordinates.
## pmapToGenome() maps the i-th element of 'x' to the i-th element of 
## 'alignment'.
first_loc <- do.call(c, sapply(cag_loc, "[", 1, simplify=TRUE))
pmapToGenome(first_loc, exons[[1]])


## -----------------------------------------------------------------------
## B. Map 3'UTR variants to genome coordinates
## -----------------------------------------------------------------------

## A study by Skeeles et. al (PLoS ONE 8(3): e58609. doi:
## 10.1371/journal.pone.0058609) investigated the impact of 3'UTR variants 
## on the expression of cancer susceptibility genes.

## 8 candidate miRNA genes on chromosome 12 were used to test for 
## differential luciferase expression in mice. In Table 2 of the manuscript
## variant locations are given as nucleotide position within the gene.
geneNames <- c("Bcap29", "Dgkb", "Etv1", "Hbp1", "Hbp1", "Ifrd1", "Ifrd1", 
               "Pik3cg", "Pik3cg", "Tspan13", "Twistnb")
starts <- c(1409, 3170, 3132, 2437, 2626, 3239, 3261, 4947, 4979, 958, 1489) 
snps <- GRanges("chr12", IRanges(starts, width=1, names=geneNames))

## To map these transcript-space coordinates to the genome we need genome-based
## gene ranges.
library(org.Mm.eg.db)
geneid <- select(org.Mm.eg.db, unique(geneNames), "ENTREZID", "SYMBOL")
geneid

## Extract the gene regions:
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
genes <- genes(txdb)[geneid$ENTREZID]

## mapToGenome() determines which pairs to match by comparing names in 'x' 
## and 'alignment'. Ranges in 'snps' will be mapped to all ranges in 'genes'
## with the same name. Currently the names of 'genes' are internal gene ids.
## Rename 'genes' with the appropriate gene symbol.
names(genes) <- geneid$SYMBOL 

## The xHits and alignmentHits columns in the output describe which ranges in
## 'snps' was mapped to which range in 'alignment'.
mapToGenome(snps, genes)


## -----------------------------------------------------------------------
## C. Map dbSNP variants to CDS or cDNA coordinates
## -----------------------------------------------------------------------

## The GIPR gene encodes a G-protein coupled receptor for gastric inhibitory 
## polypeptide (GIP). Originally GIP was identified to inhibited gastric acid 
## secretion and gastrin release but was later demonstrated to stimulate 
## insulin release in the presence of elevated glucose.

## In this example 5 SNPs located in the GIPR gene are mapped to cDNA 
## coordinates. A list of SNPs in GIPR can be downloaded from dbSNP or NCBI.
rsids <- c("rs4803846", "rs139322374", "rs7250736", "rs7250754", "rs9749185")

## Extract genomic coordinates with a SNPlocs package.
library(SNPlocs.Hsapiens.dbSNP141.GRCh38)
snps <- snpid2grange(SNPlocs.Hsapiens.dbSNP141.GRCh38, rsids)

## Gene regions of GIPR can be extracted from a TxDb package of compatible
## build. The TxDb package uses Entrez gene identifiers and GIPR is a gene 
## symbol. Conversion between gene symbols and Entrez gene IDs is done by 
## calling select() on an organism db package.
library(org.Hs.eg.db)
geneid <- select(org.Hs.eg.db, "GIPR", "ENTREZID", "SYMBOL")

## The transcriptsBy() extractor returns a single range for each transcript
## which includes the UTR and exon regions (i.e., cDNA).
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
txbygene <- transcriptsBy(txdb, "gene")
cDNA <- txbygene[geneid$ENTREZID]
cDNA

## Before mapping, the chromosome names (seqlevels) in the two objects must 
## be harmonized. The style for 'snps' is dbSNP and 'cDNA' is UCSC.
seqlevelsStyle(snps)
seqlevelsStyle(cDNA)

## Modify the style and genome in 'snps' to match 'cDNA'.
seqlevelsStyle(snps) <- seqlevelsStyle(cDNA)
genome(snps) <- genome(cDNA)

## When mapping to transcript-centric coordinates we unlist 'cDNA' because
## each range is a separte transcript. In the CDS mappings 'alignment' is 
## left as a GRangesList because all ranges belong to the CDS region.

## Map all 5 SNPS to all 4 transcripts:
mapToTranscript(snps, unlist(cDNA))

## Map the first SNP to transcript uc002pct.1 and the second to 
## transcript uc002pcu.1.
pmapToTranscript(snps[1:2], unlist(cDNA)[1:2])

## The cdsBy() extractor returns coding regions by gene or by transcript.
## Extract the coding regions for transcript uc002pct.1.
cds <- cdsBy(txdb, "tx", use.names=TRUE)["uc002pct.1"]
cds

## Map all 5 SNPs to the CDS region: 
mapToTranscript(snps, cds)

## Only the second SNP could be mapped. Unlisting the 'cds' object maps the
## SNPs to the individual CDS ranges (vs the concatenated range).
mapToTranscript(snps[2], unlist(cds))

## The location is the same because the SNP hit the first CDS range. If the
## transcript had been on the negative strand the difference in mapping
## would be more obvious.
strand(snps) <- strand(cds) <- "-"
mapToTranscript(snps[2], cds, ignore.strand=FALSE)
mapToTranscript(snps[2], unlist(cds), ignore.strand=FALSE)
}

\keyword{methods}
\keyword{utilities}
