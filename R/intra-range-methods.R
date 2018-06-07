### =========================================================================
### Intra-range methods
### -------------------------------------------------------------------------
###
### The methods implemented in this file should behave consistently with
### those defined in IRanges intra-range-methods.R
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### update_ranges()
###

setMethod("update_ranges", "GenomicRanges",
    function(x, start=NULL, end=NULL, width=NULL, use.names=TRUE)
    {
        new_ranges <- update_ranges(ranges(x), start=start,
                                               end=end,
                                               width=width,
                                               use.names=use.names)
        update(x, ranges=new_ranges)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### resize()
###

setMethod("resize", "GenomicRanges",
    function(x, width, fix="start", use.names=TRUE, ignore.strand=FALSE)
    {
        if (!missing(fix) && length(x) > 0L && 
            (length(fix) > length(x) || length(x) %% length(fix) > 0L))
            stop("'x' is not a multiple of 'fix' length")
        if (!isTRUEorFALSE(ignore.strand))
            stop("'ignore.strand' must be TRUE or FALSE")
        if (ignore.strand) {
            fix <- Rle(fix, length(x))
        } else {
            revFix <- c(start="end", end="start", center="center")
            if (length(x) == 0L)
              fix <- character()
            else fix <- ifelse(strand(x) == "-", revFix[fix], fix)
        }
        ## For some unclear reason (likely a bug in callNextMethod) we need
        ## to explicitly pass the arguments to the call below, otherwise
        ## it seems that the original unmodified 'start' gets passed.
        #callNextMethod()
        callNextMethod(x, width, fix=fix, use.names=use.names)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### flank()
###

setMethod("flank", "GenomicRanges", 
    function(x, width, start=TRUE, both=FALSE, use.names=TRUE,
                ignore.strand=FALSE)
    {
        if (!isTRUEorFALSE(ignore.strand))
            stop("'ignore.strand' must be TRUE or FALSE")
        if (!is.logical(start) || S4Vectors:::anyMissing(start))
            stop("'start' must be logical without NA's")
        start <- S4Vectors:::recycleVector(unname(start), length(x))
        if (!ignore.strand)
            start <- as.vector(start != (strand(x) == "-"))
        ## For some unclear reason (likely a bug in callNextMethod) we need
        ## to explicitly pass the arguments to the call below, otherwise
        ## it seems that the original unmodified 'start' gets passed.
        #callNextMethod()
        callNextMethod(x, width, start=start, both=both, use.names=use.names)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### promoters()
###

setMethod("promoters", "GenomicRanges",
    function(x, upstream=2000, downstream=200, use.names=TRUE)
    {
        x_len <- length(x)
        upstream <- recycleIntegerArg(upstream, "upstream", x_len)
        downstream <- recycleIntegerArg(downstream, "downstream", x_len)
        if (x_len == 0) {
            if (!S4Vectors:::normargUseNames(use.names))
                names(x) <- NULL
            return(x)
        }
        if (min(upstream) < 0L || min(downstream) < 0L)
            stop("'upstream' and 'downstream' must be integers >= 0")
        x_ranges <- ranges(x)
        x_start <- start(x_ranges)
        x_end <- end(x_ranges)
        strand_is_minus <- strand(x) == "-"
        on_plus <- which(!strand_is_minus)
        on_plus_TSS <- x_start[on_plus]
        x_start[on_plus] <- on_plus_TSS - upstream[on_plus]
        x_end[on_plus] <- on_plus_TSS + downstream[on_plus] - 1L
        on_minus <- which(strand_is_minus)
        on_minus_TSS <- x_end[on_minus]
        x_end[on_minus] <- on_minus_TSS + upstream[on_minus]
        x_start[on_minus] <- on_minus_TSS - downstream[on_minus] + 1L
        new_ranges <- update_ranges(x_ranges, start=x_start,
                                              end=x_end,
                                              use.names=use.names)
        update(x, ranges=new_ranges)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### reflect()
###

### TODO: Add "reflect" method for GenomicRanges objects.


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### restrict()
###

.checkParms <- function(x, parm)
{
    if (!all(is.na(parm))) {
        if (!all(names(parm) %in% levels(seqnames(x))))
            stop("start should be a named numeric vector ",
                 "corresponding to seqnames")
    }
    temp <- structure(rep(NA_integer_, length(levels(seqnames(x)))), 
                      names=levels(seqnames(x)))
    temp[names(parm)] <- parm
    temp
}

.restrictRngs <- function(x, start, end, keep.all.ranges, use.names)
{
    tmp <- names(x)
    names(x) <- seq_len(length(x))
    rng <- ranges(x)
    res <- restrict(ranges(x), start, end, keep.all.ranges, use.names=TRUE)
    x <- x[as.numeric(names(res))]
    ranges(x) <- res
    if (!use.names)
        names(x) <- NULL
    else 
        names(x) <- tmp[as.numeric(names(res))]
    x
}

setMethod("restrict", "GenomicRanges",
    function(x, start=NA, end=NA, keep.all.ranges=FALSE, use.names=TRUE)
    {
        if (is.null(names(start)) && is.null(names(end)))
            return(.restrictRngs(x, start, end,keep.all.ranges, use.names))
        x_mcols <- mcols(x, use.names=FALSE)
        nms <- names(x_mcols)
        mcols(x) <- cbind(x_mcols, DataFrame(posIndx=seq_len(length(x))))
        splt <- split(x, seqnames(x))
        start <- .checkParms(x, start)
        end <- .checkParms(x, end) 
        res <- lapply(names(splt),
                      function(i) {
                          .restrictRngs(splt[[i]], start=start[i], end=end[i],
                                        keep.all.ranges, use.names)
                      })
        names(res) <- names(splt)
        ord <- unlist(GRangesList(res), use.names=FALSE)
        df <- data.frame(orig=mcols(ord)[["posIndx"]],
                         final=seq_len(length(ord)))
        indx <- order(df[["orig"]], df[["final"]])
        ord <- ord[indx, ]
        mcols(ord) <- subset(mcols(ord, use.names=FALSE), select=nms)
        ord
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### trim()
###

setMethod("trim", "GenomicRanges",
    function(x, use.names=TRUE)
    {
        ## We trim only out-of-bound ranges located on non-circular sequences
        ## whose length is not NA.
        ## See get_out_of_bound_index() in GenomicRanges-class.R
        idx <- get_out_of_bound_index(x)
        if (length(idx) == 0L)
            return(x)
        new_ranges <- ranges(x)
        seqnames_id <- as.integer(seqnames(x))[idx]
        new_end <- unname(seqlengths(x))[seqnames_id]
        new_ranges[idx] <- restrict(new_ranges[idx], start=1L, end=new_end,
                                                     keep.all.ranges=TRUE)
        update(x, ranges=new_ranges)
    }
)

