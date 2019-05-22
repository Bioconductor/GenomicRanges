### =========================================================================
### Manipulate genomic variables i.e. data/variable defined along a genome
### -------------------------------------------------------------------------
###
### The concept of genomic variables could be formalized via a dedicated
### container. This container could be a simple extension of GRanges with no
### additional slots and the following constraints:
###   - The ranges are unstranded (i.e. strand is set to * for all ranges).
###   - The ranges are disjoint and ordered.
### Then the metadata columns are the genomic variables.
### For now, we just use a GRanges object. We make sure it's disjoint and we
### ignore its strand. We don't mind if it's not ordered and make sure that
### the code that operates on it works properly even if it's not ordered.
###

setAs("RleList", "GRanges",
    function(from)
    {
        what <- "RleList object to coerce to GRanges"
        from_names <- names(from)
        if (is.null(from_names))
            stop(what, " must have names")
        msg <- GenomeInfoDb:::.valid.Seqinfo.seqnames(
                                   from_names,
                                   what=paste("names of", what))
        if (!is.null(msg))
            stop(wmsg(msg))
        from_runlens <- runLength(from)
        nrun <- lengths(from_runlens, use.names=FALSE)
        ans_seqnames <- Rle(factor(from_names, levels=from_names), nrun)
        ans_width <- unlist(from_runlens, use.names=FALSE)
        ans_end <- unlist(cumsum(from_runlens), use.names=FALSE)
        ans_ranges <- IRanges(end=ans_end, width=ans_width)
        score <- unlist(runValue(from), use.names=FALSE)
        GRanges(ans_seqnames,
                ans_ranges,
                score=score,
                seqlengths=lengths(from))
    }
)

### Works only on a RleViews object with no out-of-limits views.
.from_RleViews_to_IRanges_with_score_and_view_mcols <- function(from)
{
    unlisted_from <- unlist(from, use.names=FALSE)  # will fail if 'from' has
                                                    # out-of-limits views
    from_ranges <- ranges(from)
    q <- PartitioningByWidth(runLength(unlisted_from))
    s <- PartitioningByWidth(from_ranges)
    hits <- findOverlaps(q, s, minoverlap=1L)
    ans <- from_ranges[subjectHits(hits)]
    Ltrim <- pmax(start(q)[queryHits(hits)] - start(s)[subjectHits(hits)], 0L)
    Rtrim <- pmax(end(s)[subjectHits(hits)] - end(q)[queryHits(hits)], 0L)
    ans <- windows(ans, start=1L+Ltrim, end=-1L-Rtrim)
    score <- runValue(unlisted_from)[queryHits(hits)]
    view <- subjectHits(hits)
    mcols(ans) <- DataFrame(score=score, view=view)
    ans
}

### Return a GRanges object with "score" and "views" metadata columns on it.
setAs("RleViewsList", "GRanges",
    function(from)
    {
        irl <- IRangesList(
            lapply(from, .from_RleViews_to_IRanges_with_score_and_view_mcols)
        )
        ans <- as(irl, "GRanges")
        mcols(ans) <- mcols(unlist(irl, use.names=FALSE), use.names=FALSE)
        ans
    }
)

### Represent a collection of named RleList objects as a GRanges with 1
### metadata column per RleList object.
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
    DFL <- cbind(...)  # named CompressedSplitDataFrameList
    unlisted_DFL <- unlist(DFL, use.names=FALSE)  # DataFrame
    DFL_partitioning <- PartitioningByEnd(DFL)

    ## Prepare 'ans_seqnames'.
    ans_seqnames <- Rle(factor(ans_seqlevels, levels=ans_seqlevels),
                        width(DFL_partitioning))

    ## Prepare 'ans_ranges'.
    ans_width <- unlisted_DFL[ , "runLength"]
    width_list <- relist(ans_width, DFL)
    ans_end <- unlist(lapply(width_list, cumsum), use.names=FALSE)
    ans_ranges <- IRanges(end=ans_end, width=ans_width)

    ## Prepare 'ans_seqlengths'.
    ans_seqlengths <- setNames(ans_end[end(DFL_partitioning)], names(DFL))

    ## First column is "runLength". Get rid of it.
    ans_mcols <- unlisted_DFL[-1L]
    ans <- new_GRanges("GRanges", ans_seqnames, ans_ranges,
                                  mcols=ans_mcols, seqinfo=ans_seqlengths)

    ## Keep only ranges for which at least one variable is not NA.
    keep_idx <- which(rowSums(!is.na(mcols(ans, use.names=FALSE))) != 0L)
    ans[keep_idx]
}

### Return a named RleList with 1 list element per seqlevel in 'x'.
### Works on any metadata column that can be put in Rle form (i.e. any atomic
### vector or factor).
mcolAsRleList <- function(x, varname)
{
    if (!is(x, "GenomicRanges"))
        stop("'x' must be a GRanges object")
    var <- mcols(x, use.names=FALSE)[ , varname]

    ## If 'var' is numeric, then we can use coverage().
    #if (is.numeric(var))
    #    return(coverage(x, weight=var))

    ## Otherwise 'x' must be disjoint and we compute the RleList in a loop.
    ## This would also work on a numeric metadata column but would be slower
    ## than using coverage(), especially if 'x' has many seqlevels.
    if (!isDisjoint(x, ignore.strand=TRUE))
        stop(wmsg("cannot turn non-numeric metadata column into a ",
                  "named RleList object when 'x' is not disjoint ",
                  "(ignoring the strand)"))
    rg_per_chrom <- split(ranges(x), seqnames(x))
    var_per_chrom <- split(var, seqnames(x))
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
        var_per_chrom,
        SIMPLIFY=FALSE
    )
    as(rle_list, "SimpleRleList")
}

binnedAverage <- function(bins, numvar, varname, na.rm=FALSE)
{
    if (!is(bins, "GRanges"))
        stop("'x' must be a GRanges object")
    if (!is(numvar, "RleList"))
        stop("'numvar' must be an RleList object")
    if (!identical(seqlevels(bins), names(numvar)))
        stop("'seqlevels(bin)' and 'names(numvar)' must be identical")

    ## A version of viewMeans() that pads "out of limits" views with zeros.
    viewMeans2 <- function(v, na.rm=FALSE) {
        if (!isTRUEorFALSE(na.rm))
            stop("'na.rm' must be TRUE or FALSE")
        means <- viewMeans(v, na.rm=na.rm)
        w0 <- width(v)
        v1 <- trim(v)
        w1 <- width(v1)
        if (na.rm) {
            na_count <- sum(is.na(v1))
            w0 <- w0 - na_count
            w1 <- w1 - na_count
        }
        means <- means * w1 / w0
        means[w0 != 0L & w1 == 0L] <- 0
        means
    }

    bins_per_chrom <- split(ranges(bins), seqnames(bins))
    means_list <- lapply(names(numvar),
        function(seqname) {
            v <- Views(numvar[[seqname]], bins_per_chrom[[seqname]])
            viewMeans2(v, na.rm=na.rm)
        })
    new_mcol <- unsplit(means_list, as.factor(seqnames(bins)))
    mcols(bins)[[varname]] <- new_mcol
    bins
}

