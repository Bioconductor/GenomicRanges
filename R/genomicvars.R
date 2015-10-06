### =========================================================================
### Manipulate genomic variables i.e. data/variable defined along a genome
### -------------------------------------------------------------------------


setAs("RleList", "GRanges", function(from) {
  rd <- as(from, "RangedData")
  rd$strand <- "*"
  gr <- as(rd, "GRanges")
  seqlengths(gr) <- elementLengths(from)
  gr
})

setAs("RleViewsList", "GRanges", function(from) {
  as(as(from, "RangedData"), "GRanges")
})

### Represent a collection of named RleList objects as a GRanges with 1
### metadata column per RleList object.
### When used on a single RleList object 'x', 'bindAsGRanges(x)' is equivalent
### to 'as(x, "GRanges")' (except for the name of the metadata column).
bindAsGRanges <- function(...)
{
    args <- list(...)
    if (length(args) == 0L)
        stop("nothing to bind")
    ## TODO: Implement (in C) fast 'elementIs(objects, class)' in S4Vectors
    ## that does 'sapply(objects, is, class, USE.NAMES=FALSE)', and use it
    ## here.
    if (!all(sapply(args, is, "RleList", USE.NAMES=FALSE)))
        stop("the objects to bind must be RleList objects")
    ans_seqlevels <- names(args[[1L]])
    if (is.null(ans_seqlevels))
        stop("the RleList objects to combine must have names")
    if (any(ans_seqlevels %in% c(NA_character_, ""))
     || anyDuplicated(ans_seqlevels))
        stop(wmsg("the names on the RleList objects cannot contain NAs, ",
                  "empty strings, or duplicates"))
    if (!all(sapply(args[-1L], function(arg)
                               identical(names(arg), ans_seqlevels))))
        stop(wmsg("the RleList objects to combine must have the same length ",
                  "and the same names in the same order"))
    df <- cbind(...)

    ## Prepare 'ans_seqnames'.
    ans_seqnames <- Rle(factor(df[ , "group_name"], levels=ans_seqlevels))

    ## Prepare 'ans_ranges'.
    ans_width <- df[ , "runLength"]
    partitioning <- PartitioningByEnd(df[ , "group"], NG=length(ans_seqlevels),
                                      names=ans_seqlevels)
    width_list <- relist(ans_width, partitioning)
    ans_end <- unlist(lapply(width_list, cumsum), use.names=FALSE)
    ans_ranges <- IRanges(end=ans_end, width=ans_width)

    ## Prepare 'ans_seqlengths'.
    ans_seqlengths <- setNames(ans_end[end(partitioning)], names(partitioning))

    ## First 3 columns are group, group_name, and runLength. Get rid of them.
    ans_mcols <- as(df[-(1:3)], "DataFrame")
    newGRanges("GRanges", ans_seqnames, ans_ranges, mcols=ans_mcols,
               seqlengths=ans_seqlengths)
}

### Return a named RleList with 1 list-element per seqlevel in 'x'.
### Works on any metadata column that can be put in Rle form (i.e. any atomic
### vector or factor).
mcolAsRleList <- function(x, mcolname)
{
    if (!is(x, "GenomicRanges"))
        stop("'x' must be a GRanges object")
    V <- mcols(x)[ , mcolname]

    ## If 'V' is numeric, then we can use coverage().
    if (is.numeric(V))
        return(coverage(x, weight=V))

    ## Otherwise 'x' must be disjoint and we compute the RleList in a loop.
    ## This would also work on a numeric metadata column but would be slower
    ## than using coverage(), especially if 'x' has many seqlevels.
    if (!isDisjoint(x))
        stop(wmsg("cannot turn non-numeric metadata column into a ",
                  "named RleList when 'x' is not disjoint"))
    rg_per_chrom <- split(ranges(x), seqnames(x))
    V_per_chrom <- split(V, seqnames(x))
    rle_list <- mapply(
        function(seqlen, ir, v) {
            if (is.na(seqlen))
                seqlen <- max(end(ir))
            rle <- Rle(v[NA_integer_], seqlen)
            rle[ir] <- Rle(v, width(ir))
            rle
        },
        seqlengths(x),
        rg_per_chrom,
        V_per_chrom,
        SIMPLIFY=FALSE
    )
    as(rle_list, "SimpleRleList")
}

binnedAverage <- function(bins, numvar, mcolname)
{
    if (!is(bins, "GRanges"))
        stop("'x' must be a GRanges object")
    if (!is(numvar, "RleList"))
        stop("'numvar' must be an RleList object")
    if (!identical(seqlevels(bins), names(numvar)))
        stop("'seqlevels(bin)' and 'names(numvar)' must be identical")

    ## A version of viewMeans() that pads "out of limits" views with zeros.
    viewMeans2 <- function(v) {
        means <- viewMeans(v)
        w0 <- width(v)
        w1 <- width(trim(v))
        means <- means * w1 / w0
        means[w0 != 0L & w1 == 0L] <- 0
        means
    }

    bins_per_chrom <- split(ranges(bins), seqnames(bins))
    means_list <- lapply(names(numvar),
        function(seqname) {
            v <- Views(numvar[[seqname]], bins_per_chrom[[seqname]])
            viewMeans2(v)
        })
    new_mcol <- unsplit(means_list, as.factor(seqnames(bins)))
    mcols(bins)[[mcolname]] <- new_mcol
    bins
}

