### =========================================================================
### subtract()
### -------------------------------------------------------------------------


### Similar to bedtools subtract.
### TODO: Move definition of generic to IRanges package and implement method
### for IntegerRanges objects.
setGeneric("subtract", signature=c("x", "y"),
    function(x, y, minoverlap=1L, ...) standardGeneric("subtract")
)

### Returns a GenomicRangesList object parallel to 'x'.
setMethod("subtract", c("GenomicRanges", "GenomicRanges"),
    function(x, y, minoverlap=1L, ignore.strand=FALSE)
    {
        y <- reduce(y, ignore.strand=ignore.strand)
        hits <- findOverlaps(x, y, minoverlap=minoverlap,
                                   ignore.strand=ignore.strand)
        ans <- psetdiff(x, extractList(y, as(hits, "IntegerList")))
        mcols(ans) <- mcols(x)
        setNames(ans, names(x))
    }
)

