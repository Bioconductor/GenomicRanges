### =========================================================================
### GNCList objects
### -------------------------------------------------------------------------


setClass("GNCList",
    contains="GenomicRanges",
    representation(
        nclists="NCLists",
        space="integer",
        granges="GRanges"
    )
)

setMethod("granges", "GNCList",
    function(x, use.mcols=FALSE)
    {
        if (!isTRUEorFALSE(use.mcols))
            stop("'use.mcols' must be TRUE or FALSE")
        ans <- x@granges
        if (use.mcols)
            mcols(ans) <- mcols(x)
        ans
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Accessors
###

setMethod("length", "GNCList", function(x) length(granges(x)))
setMethod("names", "GNCList", function(x) names(granges(x)))
setMethod("seqnames", "GNCList", function(x) seqnames(granges(x)))
setMethod("start", "GNCList", function(x, ...) start(granges(x)))
setMethod("end", "GNCList", function(x, ...) end(granges(x)))
setMethod("width", "GNCList", function(x) width(granges(x)))
setMethod("ranges", "GNCList", function(x, use.mcols=FALSE) ranges(granges(x)))
setMethod("strand", "GNCList", function(x) strand(granges(x)))
setMethod("seqinfo", "GNCList", function(x) seqinfo(granges(x)))

setAs("GNCList", "GRanges", function(from) granges(x))


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor
###

GNCList <- function(x)
{
    if (!is(x, "GenomicRanges"))
        stop("'x' must be a GenomicRanges object")
    if (!is(x, "GRanges"))
        x <- as(x, "GRanges")
    ans_space <- 3L * (as.integer(seqnames(x)) - 1L) + as.integer(strand(x))
    ans_nclists <- NCLists(split(ranges(x), ans_space))
    new2("GNCList", nclists=ans_nclists, space=ans_space, granges=x,
         check=FALSE)
}

