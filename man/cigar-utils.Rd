\name{cigar-utils}
\Rdversion{1.1}
\alias{cigar-utils}

\alias{validCigar}
\alias{cigarOpTable}
\alias{cigarToQWidth}
\alias{cigarToWidth}
\alias{cigarQNarrow}
\alias{cigarNarrow}

\alias{cigarToIRanges}
\alias{cigarToIRangesListByAlignment}
\alias{cigarToIRangesListByRName}

\alias{queryLoc2refLoc}
\alias{queryLocs2refLocs}

\alias{splitCigar}
\alias{cigarToRleList}
\alias{cigarToCigarTable}
\alias{summarizeCigarTable}

\title{
  CIGAR utility functions
}
\description{
  Utility functions for low-level CIGAR manipulation.
}
\usage{
cigarOpTable(cigar)

cigarToQWidth(cigar, before.hard.clipping=FALSE)
cigarToWidth(cigar)

cigarQNarrow(cigar, start=NA, end=NA, width=NA)
cigarNarrow(cigar, start=NA, end=NA, width=NA)

cigarToIRanges(cigar,
               drop.D.ranges=FALSE, drop.empty.ranges=FALSE,
               reduce.ranges=TRUE)

cigarToIRangesListByAlignment(cigar, pos, flag=NULL,
                              drop.D.ranges=FALSE, drop.empty.ranges=FALSE,
                              reduce.ranges=TRUE)

cigarToIRangesListByRName(cigar, rname, pos, flag=NULL,
                          drop.D.ranges=FALSE, drop.empty.ranges=FALSE,
                          reduce.ranges=TRUE)

queryLoc2refLoc(qloc, cigar, pos=1)
queryLocs2refLocs(qlocs, cigar, pos, flag=NULL)

splitCigar(cigar)
cigarToRleList(cigar)
cigarToCigarTable(cigar)
summarizeCigarTable(x)
}

\arguments{
  \item{cigar}{
    A character vector/factor containing the extended CIGAR
    string for each read.
    For \code{cigarToIRanges} and \code{queryLoc2refLoc},
    this must be a single string (i.e. a character vector/factor
    of length 1).
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
  \item{drop.D.ranges}{
    Should the ranges corresponding to a deletion from the
    reference (encoded with a D in the CIGAR) be dropped?
    By default we keep them to be consistent with the pileup tool
    from SAMtools.
    Note that, when \code{drop.D.ranges} is \code{TRUE}, then Ds
    and Ns in the CIGAR are equivalent.
  }
  \item{drop.empty.ranges}{
    Should empty ranges be dropped?
  }
  \item{reduce.ranges}{
    Should adjacent ranges coming from the same cigar be merged or not?
    Using \code{TRUE} (the default) can significantly reduce the
    size of the returned object.
  }
  \item{pos}{
    An integer vector containing the 1-based leftmost
    position/coordinate for each (eventually clipped) read
    sequence.
  }
  \item{flag}{
    \code{NULL} or an integer vector containing the SAM flag for
    each read.
    According to the SAM specs, flag bit 0x004 has the following
    meaning: when bit 0x004 is ON then "the query sequence itself
    is unmapped".
    When \code{flag} is provided,
    \code{cigarToIRangesListByAlignment} and
    \code{cigarToIRangesListByRName} ignore these reads.
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
  \item{x}{
    A \link[IRanges]{DataFrame} produced by \code{cigarToCigarTable}.
  }
}

\value{
  For \code{cigarOpTable}: An integer matrix with number of rows equal
  to the length of \code{cigar} and seven columns, one for each extended
  CIGAR operation.

  For \code{cigarToQWidth}: An integer vector of the same length
  as \code{cigar} where each element is the width of the query (i.e.
  the length of the query sequence) as inferred from the corresponding
  element in \code{cigar} (NAs in \code{cigar} will produce NAs in the
  returned vector).

  For \code{cigarQNarrow} and \code{cigarNarrow}: A character vector
  of the same length as \code{cigar} containing the narrowed cigars.
  In addition the vector has an "rshift" attribute which is an integer
  vector of the same length as \code{cigar}. It contains the values that
  would need to be added to the POS field of a SAM/BAM file as a
  consequence of this cigar narrowing.

  For \code{cigarToWidth}: An integer vector of the same length
  as \code{cigar} where each element is the width of the alignment
  (i.e. its total length on the reference, gaps included)
  as inferred from the corresponding element in \code{cigar}
  (NAs in \code{cigar} will produce NAs in the returned vector).

  For \code{cigarToIRanges}: An \link[IRanges]{IRanges} object
  describing where the bases in the read align with respect to
  an imaginary reference sequence assuming that the leftmost
  aligned base is at position 1 in the reference (i.e. at the
  first position).

  For \code{cigarToIRangesListByAlignment}: A
  \link[IRanges]{CompressedIRangesList} object of the same length
  as \code{cigar}.

  For \code{cigarToIRangesListByRName}: A named \link[IRanges]{IRangesList}
  object with one element (\link[IRanges]{IRanges}) per unique reference
  sequence.

  For \code{queryLoc2refLoc}: An integer vector of the same length as
  \code{qloc} containing the "reference-based locations" (i.e. the
  1-based locations relative to the reference sequence) corresponding
  to the "query-based locations" passed in \code{qloc}.

  For \code{queryLocs2refLocs}: A list of the same length as
  \code{qlocs} where each element is an integer vector containing
  the "reference-based locations" corresponding to the "query-based
  locations" passed in the corresponding element in \code{qlocs}.

  For \code{splitCigar}: A list of the same length as \code{cigar}
  where each element is itself a list with 2 elements of the same
  lengths, the 1st one being a raw vector containing the CIGAR
  operations and the 2nd one being an integer vector containing
  the lengths of the CIGAR operations.

  For \code{cigarToRleList}: A \link[IRanges]{CompressedRleList} object.

  For \code{cigarToCigarTable}: A frequency table of the CIGARs
  in the form of a \link[IRanges]{DataFrame} with two columns:
  \code{cigar} (a \link[IRanges]{CompressedRleList}) and \code{count}
  (an integer).

  For \code{summarizeCigarTable}: A list with two elements:
  \code{AlignedCharacters} (integer) and \code{Indels} (matrix)
}

\references{
  \url{http://samtools.sourceforge.net/}
}

\author{
  H. Pages and P. Aboyoun
}

\seealso{
  \link[IRanges]{IRanges-class},
  \link[IRanges]{IRangesList-class},
  \code{\link[IRanges]{coverage}},
  \link[IRanges]{RleList-class}
}

\examples{
## ---------------------------------------------------------------------
## A. SIMPLE EXAMPLES
## ---------------------------------------------------------------------

## With a cigar vector of length 1:
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

## With a cigar vector of length 4:
cigar2 <- c("40M", cigar1, "2S10M2000N15M", "3H25M5H")
pos <- c(1, 1001, 1,  351)
cigarToIRangesListByAlignment(cigar2, pos)
rname <- c("chr6", "chr6", "chr2", "chr6")
cigarToIRangesListByRName(cigar2, rname, pos)

cigarOpTable(cigar2)

splitCigar(cigar2)
cigarToRleList(cigar2)

cigarToCigarTable(cigar2)
cigarToCigarTable(cigar2)[,"cigar"]
cigarToCigarTable(cigar2)[,"count"]

summarizeCigarTable(cigarToCigarTable(cigar2))

## ---------------------------------------------------------------------
## B. PERFORMANCE
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
## C. COMPUTE THE COVERAGE OF THE READS STORED IN A BAM FILE
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
