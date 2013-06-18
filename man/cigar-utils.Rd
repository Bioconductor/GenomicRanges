\name{cigar-utils}

\alias{cigar-utils}

\alias{validCigar}

\alias{CIGAR_OPS}

\alias{explodeCigarOps}
\alias{explodeCigarOpLengths}
\alias{cigarToRleList}
\alias{splitCigar}

\alias{cigarRangesOnReferenceSpace}
\alias{cigarRangesOnQuerySpace}
\alias{cigarRangesOnPairwiseSpace}
\alias{cigarToIRangesListByAlignment}
\alias{cigarToIRangesListByRName}

\alias{cigarWidthOnReferenceSpace}
\alias{cigarWidthOnQuerySpace}
\alias{cigarWidthOnPairwiseSpace}

\alias{cigarNarrow}
\alias{cigarQNarrow}

\alias{cigarOpTable}

\alias{queryLoc2refLoc}
\alias{queryLocs2refLocs}

% Old stuff:
\alias{cigarToIRanges}
\alias{cigarToWidth}
\alias{cigarToQWidth}
\alias{cigarToCigarTable}
\alias{summarizeCigarTable}


\title{
  CIGAR utility functions
}

\description{
  Utility functions for low-level CIGAR manipulation.
}

\usage{
## -=-= Supported CIGAR operations =-=-
CIGAR_OPS

## -=-= Transform CIGARs into other useful representations =-=-
explodeCigarOps(cigar)
explodeCigarOpLengths(cigar)
cigarToRleList(cigar)
splitCigar(cigar)

## -=-= From CIGARs to ranges =-=-
cigarRangesOnReferenceSpace(cigar, flag=NULL,
        remove.N.regions=FALSE, pos=1L, f=NULL,
        ops=CIGAR_OPS, drop.empty.ranges=FALSE, reduce.ranges=FALSE,
        with.ops=FALSE)

cigarRangesOnQuerySpace(cigar, flag=NULL,
        before.hard.clipping=FALSE, after.soft.clipping=FALSE,
        ops=CIGAR_OPS, drop.empty.ranges=FALSE, reduce.ranges=FALSE,
        with.ops=FALSE)

cigarRangesOnPairwiseSpace(cigar, flag=NULL,
        remove.N.regions=FALSE,
        ops=CIGAR_OPS, drop.empty.ranges=FALSE, reduce.ranges=FALSE,
        with.ops=FALSE)

cigarToIRangesListByAlignment(cigar, pos=1L, flag=NULL,
                              drop.D.ranges=FALSE,
                              drop.empty.ranges=FALSE, reduce.ranges=TRUE)

cigarToIRangesListByRName(cigar, rname, pos=1L, flag=NULL,
                          drop.D.ranges=FALSE,
                          drop.empty.ranges=FALSE, reduce.ranges=TRUE)

## -=-= From CIGARs to sequence lengths =-=-
cigarWidthOnReferenceSpace(cigar, flag=NULL,
        remove.N.regions=FALSE)

cigarWidthOnQuerySpace(cigar, flag=NULL,
        before.hard.clipping=FALSE, after.soft.clipping=FALSE)

cigarWidthOnPairwiseSpace(cigar, flag=NULL,
        remove.N.regions=FALSE)

## -=-= Narrow CIGARs =-=-
cigarNarrow(cigar, start=NA, end=NA, width=NA)
cigarQNarrow(cigar, start=NA, end=NA, width=NA)

## -=-= Summarize CIGARs =-=-
cigarOpTable(cigar)

## -=-= Translate coordinates between query and reference spaces =-=-
queryLoc2refLoc(qloc, cigar, pos=1L)
queryLocs2refLocs(qlocs, cigar, pos=1L, flag=NULL)
}

\arguments{
  \item{cigar}{
    A character vector or factor containing the extended CIGAR strings.
    It can be of arbitrary length except for \code{queryLoc2refLoc} which
    only accepts a single CIGAR (as a character vector or factor of length 1).
  }
  \item{flag}{
    \code{NULL} or an integer vector containing the SAM flag for
    each read.

    According to the SAM Spec v1.4, flag bit 0x4 is the only reliable place
    to tell whether a segment (or read) is mapped (bit is 0) or not (bit
    is 1). If \code{flag} is supplied, then \code{cigarRangesOnReferenceSpace},
    \code{cigarRangesOnQuerySpace}, \code{cigarRangesOnPairwiseSpace},
    \code{cigarToIRangesListByAlignment}, and \code{cigarToIRangesListByRName}
    don't produce any range for unmapped reads i.e. they treat them as if
    their CIGAR was empty (independently of what their CIGAR is).
    If \code{flag} is supplied, then \code{cigarWidthOnReferenceSpace},
    \code{cigarWidthOnQuerySpace}, and \code{cigarWidthOnPairwiseSpace}
    return \code{NA}s for unmapped reads.
  }
  \item{remove.N.regions}{
    \code{TRUE} or \code{FALSE}.
    If \code{TRUE}, then \code{cigarRangesOnReferenceSpace} and
    \code{cigarWidthOnReferenceSpace} report ranges/widths with respect
    to the "reference" space from which the N regions have been removed,
    and \code{cigarRangesOnPairwiseSpace} and \code{cigarWidthOnPairwiseSpace}
    report them with respect to the "pairwise" space from which the N
    regions have been removed.
  }
  \item{pos}{
    An integer vector containing the 1-based leftmost
    position/coordinate for each (eventually clipped) read
    sequence. Must have length 1 (in which case it's recycled to the
    length of \code{cigar}), or the same length as \code{cigar}.
  }
  \item{f}{
    \code{NULL} or a factor of length \code{cigar}.
    If \code{NULL}, then the ranges are grouped by alignment i.e. the
    returned \link[IRanges]{IRangesList} object has 1 list element per
    element in \code{cigar}. Otherwise they are grouped by factor level
    i.e. the returned \link[IRanges]{IRangesList} object has 1 list element
    per level in \code{f} and is named with those levels.

    For example, if \code{f} is a factor containing the chromosome for each
    read, then the returned \link[IRanges]{IRangesList} object will have
    1 list element per chromosome and each list element will contain all
    the ranges on that chromosome.
  }
  \item{ops}{
    Character vector containing the CIGAR operations to translate to ranges.
    By default, all CIGAR operations are translated.
  }
  \item{drop.empty.ranges}{
    Should empty ranges be dropped?
  }
  \item{reduce.ranges}{
    Should adjacent ranges coming from the same cigar be merged or not?
    Using \code{TRUE} can significantly reduce the size of the returned
    object.
  }
  \item{with.ops}{
    \code{TRUE} or \code{FALSE} indicating whether the returned ranges should
    be named with their corresponding CIGAR operation.
  }
  \item{before.hard.clipping}{
    \code{TRUE} or \code{FALSE}.
    If \code{TRUE}, then \code{cigarRangesOnQuerySpace} and
    \code{cigarWidthOnQuerySpace} report ranges/widths with respect
    to the "query" space to which the H regions have been added.
    \code{before.hard.clipping} and \code{after.soft.clipping} cannot
    both be \code{TRUE}.
  }
  \item{after.soft.clipping}{
    \code{TRUE} or \code{FALSE}.
    If \code{TRUE}, then \code{cigarRangesOnQuerySpace} and
    \code{cigarWidthOnQuerySpace} report ranges/widths with respect
    to the "query" space from which the S regions have been removed.
    \code{before.hard.clipping} and \code{after.soft.clipping} cannot
    both be \code{TRUE}.
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
    A factor of length \code{cigar} containing the name of the
    reference sequence associated with each read (i.e. the
    name of the sequence the read has been aligned to).
    It's treated as the \code{f} argument above.
  }
  \item{start,end,width}{
    Vectors of integers. NAs and negative values are accepted and
    "solved" according to the rules of the SEW (Start/End/Width)
    interface (see \code{?\link[IRanges]{solveUserSEW}} for the details).
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

  For \code{cigarRangesOnReferenceSpace}, \code{cigarRangesOnQuerySpace},
  and \code{cigarRangesOnPairwiseSpace}: An \link[IRanges]{IRangesList}
  object (more precisely a \link[IRanges]{CompressedIRangesList} object)
  with 1 list element per element in \code{cigar}.
  However, if \code{f} is a factor, then the returned
  \link[IRanges]{IRangesList} object can be a \link[IRanges]{SimpleIRangesList}
  object (instead of \link[IRanges]{CompressedIRangesList}), and it has
  1 list element per level in \code{f} and is named with those levels.

  For \code{cigarToIRangesListByAlignment}:
  An \link[IRanges]{IRangesList} object (more precisely a
  \link[IRanges]{CompressedIRangesList} object) of the same length
  as \code{cigar}.

  For \code{cigarToIRangesListByRName}: A named \link[IRanges]{IRangesList}
  object with one list element per level in \code{rname}.

  For \code{cigarWidthOnReferenceSpace} and \code{cigarWidthOnPairwiseSpace}:
  An integer vector of the same length as \code{cigar} where each element is
  the width of the alignment with respect to the "reference" and "pairwise"
  space, respectively. More precisely, for \code{cigarWidthOnReferenceSpace},
  the returned widths are the lengths of the alignments on the reference,
  N gaps included (except if \code{remove.N.regions} is \code{TRUE}).
  NAs or \code{"*"} in \code{cigar} will produce NAs in the returned vector.

  For \code{cigarWidthOnQuerySpace}: An integer vector of the same
  length as \code{cigar} where each element is the length of the corresponding
  query sequence as inferred from the CIGAR string. Note that, by default
  (i.e. if \code{before.hard.clipping} and \code{after.soft.clipping} are
  \code{FALSE}), this is the length of the query sequence stored in the
  SAM/BAM file. If \code{before.hard.clipping} or \code{after.soft.clipping}
  is \code{TRUE}, the returned widths are the lengths of the query sequences
  before hard clipping or after soft clipping.
  NAs or \code{"*"} in \code{cigar} will produce NAs in the returned vector.

  For \code{cigarNarrow} and \code{cigarQNarrow}: A character vector
  of the same length as \code{cigar} containing the narrowed cigars.
  In addition the vector has an "rshift" attribute which is an integer
  vector of the same length as \code{cigar}. It contains the values that
  would need to be added to the POS field of a SAM/BAM file as a
  consequence of this cigar narrowing.

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
    \item The \link[Rsamtools]{sequenceLayer} function in the \pkg{Rsamtools}
          package for laying the query sequences to the "reference" or
          "pairwise" spaces.

    \item The \link{GAlignments} container for storing a set of genomic
          alignments.

    \item The \link[IRanges]{IRanges}, \link[IRanges]{IRangesList}, and
          \link[IRanges]{RleList} classes in the \pkg{IRanges} package.

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

## cigarRangesOnReferenceSpace():
cigarRangesOnReferenceSpace(cigar1, with.ops=TRUE)[[1]]

cigarRangesOnReferenceSpace(cigar1,
                            reduce.ranges=TRUE, with.ops=TRUE)[[1]]

ops <- setdiff(CIGAR_OPS, "N")

cigarRangesOnReferenceSpace(cigar1, ops=ops, with.ops=TRUE)[[1]]

cigarRangesOnReferenceSpace(cigar1, ops=ops,
                            reduce.ranges=TRUE, with.ops=TRUE)[[1]]

ops <- setdiff(CIGAR_OPS, c("D", "N"))

cigarRangesOnReferenceSpace(cigar1, ops=ops, with.ops=TRUE)[[1]]

## cigarWidthOnReferenceSpace()/cigarWidthOnQuerySpace():
cigarWidthOnReferenceSpace(cigar1)
cigarWidthOnQuerySpace(cigar1)
cigarWidthOnQuerySpace(cigar1, before.hard.clipping=TRUE)

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

## cigarQNarrow():
cigarQNarrow(cigar1, start=4, end=-3)
cigarQNarrow(cigar1, start=10)
cigarQNarrow(cigar1, start=19)
cigarQNarrow(cigar1, start=24)

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
cigarRangesOnReferenceSpace(cigar2, pos=pos, with.ops=TRUE)
cigarRangesOnQuerySpace(cigar2, with.ops=TRUE)
cigarRangesOnPairwiseSpace(cigar2, with.ops=TRUE)

res1a <- cigarToIRangesListByAlignment(cigar2, pos=pos)
res1b <- cigarRangesOnReferenceSpace(cigar2,
                                     pos=pos,
                                     ops=setdiff(CIGAR_OPS, "N"),
                                     reduce.ranges=TRUE)
stopifnot(identical(res1a, res1b))

res2a <- cigarToIRangesListByAlignment(cigar2, pos=pos,
                                       drop.D.ranges=TRUE)
res2b <- cigarRangesOnReferenceSpace(cigar2,
                                     pos=pos,
                                     ops=setdiff(CIGAR_OPS, c("D", "N")),
                                     reduce.ranges=TRUE)
stopifnot(identical(res2a, res2b))

rname <- factor(c("chr6", "chr6", "chr2", "chr6"),
                levels=c("chr2", "chr6"))
cigarToIRangesListByRName(cigar2, rname, pos=pos)

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
  rname <- factor(sample(some_fake_rnames, nreads, replace=TRUE),
                  levels=some_fake_rnames)
  pos <- sample(80000000L, nreads, replace=TRUE)

  ## The following takes < 3 sec. to complete:
  system.time(irl1 <- cigarToIRangesListByAlignment(cigar3, pos=pos))

  ## The following takes < 4 sec. to complete:
  system.time(irl2 <- cigarToIRangesListByRName(cigar3, rname, pos=pos))

  ## The sizes of the resulting objects are about 240M and 160M,
  ## respectively:
  object.size(irl1)
  object.size(irl2)
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
  irl <- cigarToIRangesListByRName(bam$cigar, bam$rname, pos=bam$pos)
  irl <- irl[elementLengths(irl) != 0]  # drop empty elements
  coverage(irl)
}

library(Rsamtools)
f1 <- system.file("extdata", "ex1.bam", package="Rsamtools")
extractCoverageFromBAM(f1)
}

\keyword{manip}
