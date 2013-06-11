\name{cigar-utils}

\alias{cigar-utils}

\alias{validCigar}

\alias{CIGAR_OPS}

\alias{explodeCigarOps}
\alias{explodeCigarOpLengths}
\alias{cigarToRleList}
\alias{splitCigar}

\alias{cigarToQWidth}
\alias{cigarToWidth}

\alias{cigarQNarrow}
\alias{cigarNarrow}

\alias{cigarRangesOnReference}
\alias{cigarRangesOnQuery}
\alias{cigarToIRanges}
\alias{cigarToIRangesListByAlignment}
\alias{cigarToIRangesListByRName}

\alias{cigarOpTable}

\alias{queryLoc2refLoc}
\alias{queryLocs2refLocs}

% Old stuff:
\alias{cigarToCigarTable}
\alias{summarizeCigarTable}


\title{
  CIGAR utility functions
}

\description{
  Utility functions for low-level CIGAR manipulation.
}

\usage{
## Supported CIGAR operations:
CIGAR_OPS

## Transform CIGARs into other useful representations:
explodeCigarOps(cigar)
explodeCigarOpLengths(cigar)
cigarToRleList(cigar)
splitCigar(cigar)

## From CIGARs to sequence lengths:
cigarToQWidth(cigar, before.hard.clipping=FALSE)
cigarToWidth(cigar)

## Narrow CIGARs:
cigarQNarrow(cigar, start=NA, end=NA, width=NA)
cigarNarrow(cigar, start=NA, end=NA, width=NA)

## From CIGARs to ranges:
cigarRangesOnReference(cigar, pos, flag=NULL, ops=CIGAR_OPS,
                       drop.empty.ranges=FALSE, reduce.ranges=FALSE,
                       with.ops=FALSE)

cigarRangesOnQuery(cigar, flag=NULL, ops=CIGAR_OPS,
                   drop.empty.ranges=FALSE, reduce.ranges=FALSE,
                   with.ops=FALSE)

cigarToIRanges(cigar,
               drop.D.ranges=FALSE, drop.empty.ranges=FALSE,
               reduce.ranges=TRUE)

cigarToIRangesListByAlignment(cigar, pos, flag=NULL,
                              drop.D.ranges=FALSE, drop.empty.ranges=FALSE,
                              reduce.ranges=TRUE)

cigarToIRangesListByRName(cigar, rname, pos, flag=NULL,
                          drop.D.ranges=FALSE, drop.empty.ranges=FALSE,
                          reduce.ranges=TRUE)

## Summarize CIGARs:
cigarOpTable(cigar)

## Translate coordinates between query-based and reference-based:
queryLoc2refLoc(qloc, cigar, pos=1)
queryLocs2refLocs(qlocs, cigar, pos, flag=NULL)
}

\arguments{
  \item{cigar}{
    A character vector or factor containing the extended CIGAR strings.
    It can be of arbitrary length except for \code{cigarToIRanges} and
    \code{queryLoc2refLoc} which only accept a single CIGAR (as a character
    vector or factor of length 1).
  }
  \item{before.hard.clipping}{
    Should the returned widths be the lengths of the reads before
    or after "hard clipping"? Hard clipping of a read is
    encoded with an H in the CIGAR.
    If NO (\code{before.hard.clipping=FALSE}, the default),
    then the returned widths are the lengths of the query
    sequences stored in the SAM/BAM file.
    If YES (\code{before.hard.clipping=TRUE}), then the returned
    widths are the lengths of the original reads.
  }
  \item{start,end,width}{
    Vectors of integers. NAs and negative values are accepted and
    "solved" according to the rules of the SEW (Start/End/Width)
    interface (see \code{?\link[IRanges]{solveUserSEW}} for the details).
  }
  \item{pos}{
    An integer vector containing the 1-based leftmost
    position/coordinate for each (eventually clipped) read
    sequence.
  }
  \item{flag}{
    \code{NULL} or an integer vector containing the SAM flag for
    each read.

    According to the SAM Spec v1.4, flag bit 0x4 is the only reliable place
    to tell whether a segment (or read) is mapped (bit is 0) or not (bit
    is 1). If \code{flag} is supplied, then \code{cigarRangesOnReference},
    \code{cigarRangesOnQuery}, and \code{cigarToIRangesListByAlignment},
    translate unmapped reads to empty list elements in the returned
    \link[IRanges]{IRangesList} object.
    If \code{flag} is supplied, \code{cigarToIRangesListByRName} just
    ignore these reads.
  }
  \item{ops}{
    Character vector containing the CIGAR operations to translate to ranges.
    By default, all CIGAR operations are translated.
  }
  \item{with.ops}{
    \code{TRUE} or \code{FALSE} indicating whether the returned ranges should
    be named with their corresponding CIGAR operation.
  }
  \item{drop.empty.ranges}{
    Should empty ranges be dropped?
  }
  \item{reduce.ranges}{
    Should adjacent ranges coming from the same cigar be merged or not?
    Using \code{TRUE} can significantly reduce the size of the returned
    object.
  }
  \item{drop.D.ranges}{
    Should the ranges corresponding to a deletion from the
    reference (encoded with a D in the CIGAR) be dropped?
    By default we keep them to be consistent with the pileup tool
    from SAMtools.
    Note that, when \code{drop.D.ranges} is \code{TRUE}, then Ds
    and Ns in the CIGAR are equivalent.
  }
  \item{rname}{
    A character vector/factor containing the name of the
    reference sequence associated with each read (i.e. the
    name of the sequence the read has been aligned to).
  }
  \item{qloc}{
    An integer vector containing "query-based locations" i.e.
    1-based locations relative to the query sequence
    stored in the SAM/BAM file.
  }
  \item{qlocs}{
    A list of the same length as \code{cigar} where each
    element is an integer vector containing "query-based
    locations" i.e. 1-based locations relative to the corresponding
    query sequence stored in the SAM/BAM file.
  }
}

\value{
  \code{CIGAR_OPS} is a predefined character vector containing the supported
  CIGAR operations: M, I, D, N, S, H, P, =, X. See p. 4 of the SAM Spec v1.4
  at \url{http://samtools.sourceforge.net/} for the list of CIGAR operations
  and their meanings.

  For \code{explodeCigarOps} and \code{explodeCigarOpLengths}:
  Both functions return a list of the same length as \code{cigar} where each
  list element is a character vector (for \code{explodeCigarOps}) or an integer
  vector (for \code{explodeCigarOpLengths}). The 2 lists have the same shape,
  that is, same \code{length()} and same \code{elementLengths()}. The i-th
  character vector in the list returned by \code{explodeCigarOps} contains one
  single-letter string per CIGAR operation in \code{cigar[i]}. The i-th integer
  vector in the list returned by \code{explodeCigarOpLengths} contains the
  corresponding CIGAR operation lengths. Zero-length operations are ignored.

  For \code{cigarToRleList}: A \link[IRanges]{CompressedRleList} object.

  For \code{splitCigar}: A list of the same length as \code{cigar}
  where each element is itself a list with 2 elements of the same
  lengths, the 1st one being a raw vector containing the CIGAR
  operations and the 2nd one being an integer vector containing
  the lengths of the CIGAR operations.

  For \code{cigarToQWidth}: An integer vector of the same length
  as \code{cigar} where each element is the width of the query (i.e.
  the length of the query sequence) as inferred from the corresponding
  element in \code{cigar} (NAs in \code{cigar} will produce NAs in the
  returned vector).

  For \code{cigarToWidth}: An integer vector of the same length
  as \code{cigar} where each element is the width of the alignment
  (i.e. its total length on the reference, gaps included)
  as inferred from the corresponding element in \code{cigar}
  (NAs in \code{cigar} will produce NAs in the returned vector).

  For \code{cigarQNarrow} and \code{cigarNarrow}: A character vector
  of the same length as \code{cigar} containing the narrowed cigars.
  In addition the vector has an "rshift" attribute which is an integer
  vector of the same length as \code{cigar}. It contains the values that
  would need to be added to the POS field of a SAM/BAM file as a
  consequence of this cigar narrowing.

  For \code{cigarRangesOnReference} and \code{cigarRangesOnQuery}:
  An \link[IRanges]{IRangesList} object (more precisely a
  \link[IRanges]{CompressedIRangesList} object) of the same length
  as \code{cigar}.

  For \code{cigarToIRanges}: An \link[IRanges]{IRanges} object
  describing where the bases in the read align with respect to
  an imaginary reference sequence assuming that the leftmost
  aligned base is at position 1 in the reference (i.e. at the
  first position).

  For \code{cigarToIRangesListByAlignment}:
  An \link[IRanges]{IRangesList} object (more precisely a
  \link[IRanges]{CompressedIRangesList} object) of the same length
  as \code{cigar}.

  For \code{cigarToIRangesListByRName}: A named \link[IRanges]{IRangesList}
  object with one element (\link[IRanges]{IRanges}) per unique reference
  sequence.

  For \code{cigarOpTable}: An integer matrix with number of rows equal
  to the length of \code{cigar} and seven columns, one for each extended
  CIGAR operation.

  For \code{queryLoc2refLoc}: An integer vector of the same length as
  \code{qloc} containing the "reference-based locations" (i.e. the
  1-based locations relative to the reference sequence) corresponding
  to the "query-based locations" passed in \code{qloc}.

  For \code{queryLocs2refLocs}: A list of the same length as
  \code{qlocs} where each element is an integer vector containing
  the "reference-based locations" corresponding to the "query-based
  locations" passed in the corresponding element in \code{qlocs}.
}

\references{
  \url{http://samtools.sourceforge.net/}
}

\author{
  H. Pages and P. Aboyoun
}

\seealso{
  \itemize{
    \item The \link{GAlignments} container for storing a set of genomic
          alignments.
    \item The \link[IRanges]{IRanges}, \link[IRanges]{IRangesList}, and
          \link[IRanges]{RleList} classes in the IRanges package.
    \item The \code{\link[IRanges]{coverage}} generic and methods for
          computing the coverage across a set of ranges or genomic ranges.
  }
}

\examples{
CIGAR_OPS

## ---------------------------------------------------------------------
## A. ON A SINGLE CIGAR
## ---------------------------------------------------------------------
cigar1 <- "3H15M55N4M2I6M2D5M6S"

## cigarToQWidth()/cigarToWidth():
cigarToQWidth(cigar1)
cigarToQWidth(cigar1, before.hard.clipping=TRUE)
cigarToWidth(cigar1)

## cigarQNarrow():
cigarQNarrow(cigar1, start=4, end=-3)
cigarQNarrow(cigar1, start=10)
cigarQNarrow(cigar1, start=19)
cigarQNarrow(cigar1, start=24)

## cigarNarrow():
cigarNarrow(cigar1)  # only drops the soft/hard clipping
cigarNarrow(cigar1, start=10)
cigarNarrow(cigar1, start=15)
cigarNarrow(cigar1, start=15, width=57)
cigarNarrow(cigar1, start=16)
#cigarNarrow(cigar1, start=16, width=55)  # ERROR! (empty cigar)
cigarNarrow(cigar1, start=71)
cigarNarrow(cigar1, start=72)
cigarNarrow(cigar1, start=75)

## cigarToIRanges():
cigarToIRanges(cigar1)
cigarToIRanges(cigar1, reduce.ranges=FALSE)
cigarToIRanges(cigar1, drop.D.ranges=TRUE)

## ---------------------------------------------------------------------
## B. ON A VECTOR CONTAINING MORE THAN ONE CIGAR
## ---------------------------------------------------------------------
cigar2 <- c("40M", cigar1, "2S10M2000N15M", "3H25M5H")

## Transform CIGARs into other useful representations:
explodeCigarOps(cigar2)
explodeCigarOpLengths(cigar2)
cigarToRleList(cigar2)
splitCigar(cigar2)

## From CIGARs to ranges:
pos <- c(1, 1001, 1,  351)
cigarRangesOnReference(cigar2, pos, with.ops=TRUE)
cigarRangesOnQuery(cigar2, with.ops=TRUE)

res1a <- cigarToIRangesListByAlignment(cigar2, pos)
res1b <- cigarRangesOnReference(cigar2, pos,
                                ops=setdiff(CIGAR_OPS, "N"),
                                reduce.ranges=TRUE)
stopifnot(identical(res1a, res1b))

res2a <- cigarToIRangesListByAlignment(cigar2, pos, drop.D.ranges=TRUE)
res2b <- cigarRangesOnReference(cigar2, pos,
                                ops=setdiff(CIGAR_OPS, c("D", "N")),
                                reduce.ranges=TRUE)
stopifnot(identical(res2a, res2b))

rname <- c("chr6", "chr6", "chr2", "chr6")
cigarToIRangesListByRName(cigar2, rname, pos)

## Summarize CIGARs:
cigarOpTable(cigar2)

## ---------------------------------------------------------------------
## C. PERFORMANCE
## ---------------------------------------------------------------------

if (interactive()) {
  ## We simulate 20 millions aligned reads, all 40-mers. 95% of them
  ## align with no indels. 5% align with a big deletion in the
  ## reference. In the context of an RNAseq experiment, those 5% would
  ## be suspected to be "junction reads".
  set.seed(123)
  nreads <- 20000000L
  njunctionreads <- nreads * 5L / 100L
  cigar3 <- character(nreads)
  cigar3[] <- "40M"
  junctioncigars <- paste(
      paste(10:30, "M", sep=""),
      paste(sample(80:8000, njunctionreads, replace=TRUE), "N", sep=""),
      paste(30:10, "M", sep=""), sep="")
  cigar3[sample(nreads, njunctionreads)] <- junctioncigars
  some_fake_rnames <- paste("chr", c(1:6, "X"), sep="")
  rname <- sample(some_fake_rnames, nreads, replace=TRUE)
  pos <- sample(80000000L, nreads, replace=TRUE)

  ## The following takes < 5 sec. to complete:
  system.time(rglist <- cigarToIRangesListByAlignment(cigar3, pos))

  ## The following takes < 10 sec. to complete:
  system.time(irl <- cigarToIRangesListByRName(cigar3, rname, pos))

  ## Internally, cigarToIRangesListByRName() turns 'rname' into a factor
  ## before starting the calculation. Hence it will run sligthly
  ## faster if 'rname' is already a factor.
  rname2 <- as.factor(rname)
  system.time(irl2 <- cigarToIRangesListByRName(cigar3, rname2, pos))

  ## The sizes of the resulting objects are about 240M and 160M,
  ## respectively:
  object.size(rglist)
  object.size(irl)
}

## ---------------------------------------------------------------------
## D. COMPUTE THE COVERAGE OF THE READS STORED IN A BAM FILE
## ---------------------------------------------------------------------
## The information stored in a BAM file can be used to compute the
## "coverage" of the mapped reads i.e. the number of reads that hit any
## given position in the reference genome.
## The following function takes the path to a BAM file and returns an
## object representing the coverage of the mapped reads that are stored
## in the file. The returned object is an RleList object named with the
## names of the reference sequences that actually receive some coverage.

extractCoverageFromBAM <- function(file)
{
  ## This ScanBamParam object allows us to load only the necessary
  ## information from the file.
  param <- ScanBamParam(flag=scanBamFlag(isUnmappedQuery=FALSE,
                                         isDuplicate=FALSE),
                        what=c("rname", "pos", "cigar"))
  bam <- scanBam(file, param=param)[[1]]
  ## Note that unmapped reads and reads that are PCR/optical duplicates
  ## have already been filtered out by using the ScanBamParam object above.
  irl <- cigarToIRangesListByRName(bam$cigar, bam$rname, bam$pos)
  irl <- irl[elementLengths(irl) != 0]  # drop empty elements
  coverage(irl)
}

library(Rsamtools)
f1 <- system.file("extdata", "ex1.bam", package="Rsamtools")
extractCoverageFromBAM(f1)
}

\keyword{manip}
