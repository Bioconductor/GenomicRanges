## Simple functions for reading and summarizing CIGARs from BAM files

## Creates a table of CIGAR values from a BAM file
##
## Input:
##   file: The file name of the BAM file to be parsed.
##
## Output:
##   A DataFrame with two columns:
##     cigar (CompressedRleList) and count (integer)
readCigarTable <- function(file)
{
    fileData <- scanBam(file, param = ScanBamParam(what=c("cigar")))
    if (length(fileData) != 1)
        stop("scanBam generated more than 1 list element")
    cigar <-
      factor(cigars(fileData[[1]][["cigar"]])[!is.na(cigars(fileData[[1]][["cigar"]]))])
    cigarToCigarTable(cigar)
}


## Examples
if (FALSE) {

suppressMessages(library(GenomicRanges))
suppressMessages(library(Rsamtools))
dataDir <- "/home/biocdev/data_store/1000genomes/extdata"
file1 <- file.path(dataDir, "NA19239.SLX.maq.SRP000033.2009_09.subset.bam")
file2 <- file.path(dataDir, "NA19240.chrom6.SLX.maq.SRP000032.2009_07.subset.bam")
file3 <- file.path(dataDir, "NA19240.chrom1.454.ssaha2.SRP000032.2009_10.bam")

cigarTable1 <- readCigarTable(file1) ## Takes approx 1 minute
cigarTable2 <- readCigarTable(file2) ## Takes approx 1 minute
cigarTable3 <- readCigarTable(file3) ## Takes approx 6.5 minutes

summarizeCigarTable(cigarTable1)
summarizeCigarTable(cigarTable2)
summarizeCigarTable(cigarTable3)

}
