### -------------------------------------------------------------------------
### Some low-level (non exported) utility functions to operate on transcripts
### represented as groups of exon ranges.
###
### These functions are implemented in C. This file only contains R wrappers
### for the .Call2 entry points. Those wrappers are not doing any argument
### checking and therefore are considered "unsafe". They are in turn called
### by "safe" and user-friendly higher level wrappers defined in
### GenomicFeatures. The reason why the "unsafe" wrappers are here and not in
### the GenomicFeatures package was to keep GenomicFeatures free of native
### code.
###
### For all the functions below:
###   o 'exonStarts', 'exonEnds' are assumed to be lists of integer vectors.
###     The two lists are assumed to have the "same shape" i.e.
###     elementNROWS() returns identical vectors on them;
###   o 'strand' is assumed to be a character vector with allowed values
###     "+" and "-" only;
###   o 'decreasing.rank.on.minus.strand' is assumed to be TRUE or FALSE.
###   o 'error.if.out.of.bounds' is assumed to be TRUE or FALSE.

unsafe.transcriptWidths <- function(exonStarts, exonEnds)
{
    .Call2("transcript_widths",
          exonStarts, exonEnds,
          PACKAGE="GenomicRanges")
}

### 'tlocs' is assumed to be a list of integer vectors of the same length (but
### not necessarily the "same shape") as 'exonStarts' and 'exonEnds'.
unsafe.transcriptLocs2refLocs <- function(tlocs,
                exonStarts, exonEnds, strand,
                decreasing.rank.on.minus.strand, error.if.out.of.bounds)
{
    .Call2("tlocs2rlocs",
          tlocs,
          exonStarts, exonEnds, strand,
          decreasing.rank.on.minus.strand,
          error.if.out.of.bounds,
          PACKAGE="GenomicRanges")
}

