### -------------------------------------------------------------------------
### Some low-level (non exported) utility functions to operate on transcripts
### represented as groups of exon ranges
###
### For all the functions below:
###   o 'exonStarts', 'exonEnds' are assumed to be lists of integer vectors.
###     The two lists are assumed to have the "same shape" i.e.
###     elementLengths() returns identical vectors on them;
###   o 'strand' is assumed to be a character vector with allowed values
###     "+" and "-" only;
###   o 'reorder.exons.on.minus.strand' is assumed to be TRUE or FALSE.


### The safe and user-friendly wrapper to this is in GenomicFeatures.
unsafe.transcriptWidths <- function(exonStarts, exonEnds)
{
    .Call("transcript_widths",
          exonStarts, exonEnds,
          PACKAGE="GenomicRanges")
}

### The safe and user-friendly wrapper to this is in Biostrings.
unsafe.extractTranscripts <- function(classname, x,
                exonStarts, exonEnds, strand,
                reorder.exons.on.minus.strand, lkup)
{
    .Call("extract_transcripts",
          classname, x,
          exonStarts, exonEnds, strand,
          reorder.exons.on.minus.strand, lkup,
          PACKAGE="GenomicRanges")
}

### The safe and user-friendly wrapper to this is in GenomicFeatures.
### 'tlocs' is assumed to be a list of integer vectors of the same length (but
### not necessarily the "same shape") as 'exonStarts' and 'exonEnds'.
unsafe.transcriptLocs2refLocs <- function(tlocs,
                exonStarts, exonEnds, strand,
                reorder.exons.on.minus.strand)
{
    .Call("tlocs2rlocs",
          tlocs,
          exonStarts, exonEnds, strand,
          reorder.exons.on.minus.strand,
          PACKAGE="GenomicRanges")
}

