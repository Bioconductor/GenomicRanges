### =========================================================================
### Intra-range methods
### -------------------------------------------------------------------------
###
### The methods implemented in this file should behave consistently with
### those defined in IRanges intra-range-methods.R
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### shift()
###

setMethod("shift", "GenomicRanges",
    function(x, shift=0L, use.names=TRUE)
    {
        new_ranges <- shift(ranges(x), shift=shift, use.names=use.names)
        clone(x, ranges=new_ranges) 
    }
)

setMethod("shift", "GRangesList",
    function(x, shift=0L, use.names=TRUE)
    {
        ## Unlist 'x'.
        unlisted_x <- unlist(x, use.names=FALSE)
        ## Recycle and unlist 'shift'.
        if (!is(shift, "List"))
            shift <- as(shift, "List")
        shift <- S4Vectors:::VH_recycle(shift, x, "shift", "x")
        unlisted_shift <- unlist(shift, use.names=FALSE)
        ## Compute unlisted 'ans'.
        unlisted_ans <- shift(unlisted_x, shift=unlisted_shift,
                              use.names=use.names)
        ## Relist and return.
        relist(unlisted_ans, x)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### narrow()
###

setMethod("narrow", "GenomicRanges",
    function(x, start=NA, end=NA, width=NA, use.names=TRUE)
    {
        new_ranges <- narrow(ranges(x), start=start, end=end, width=width,
                                        use.names=use.names)
        clone(x, ranges=new_ranges) 
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
            fix <- Rle(rep.int(fix, length(x)))
        } else {
            revFix <- c(start="end", end="start", center="center")
            if (length(x) == 0L)
              fix <- character()
            else fix <- ifelse(strand(x) == "-", revFix[fix], fix)
        }
        new_ranges <- resize(ranges(x), width=width, fix=fix,
                                        use.names=use.names)
        clone(x, ranges=new_ranges)
    }
)

setMethod("resize", "GRangesList",
    function(x, width, fix="start", use.names=TRUE, ignore.strand=FALSE)
    {
        ## Unlist 'x'.
        unlisted_x <- unlist(x, use.names=FALSE)
        ## Recycle and unlist 'width'.
        if (!is(width, "List"))
            width <- as(width, "List")
        width <- S4Vectors:::VH_recycle(width, x, "width", "x")
        unlisted_width <- unlist(width, use.names=FALSE)
        ## Compute unlisted 'ans'.
        unlisted_ans <- resize(unlisted_x, unlisted_width, fix=fix, 
                               use.names=use.names, ignore.strand=ignore.strand)
        ## Relist and return.
        relist(unlisted_ans, x)
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
        if (ignore.strand)
            start <- rep.int(start, length(x))
        else 
            start <- as.vector(start == (strand(x) != "-"))
        new_ranges <- flank(ranges(x), width=width, start=start, both=both,
                                       use.names=use.names)
        clone(x, ranges=new_ranges) 
    }
)

setMethod("flank", "GRangesList",
    function(x, width, start=TRUE, both=FALSE, use.names=TRUE,
             ignore.strand=FALSE)
    {
        ## Unlist 'x'.
        unlisted_x <- unlist(x, use.names=FALSE)
        ## Recycle and unlist 'width'.
        if (!is(width, "List"))
            width <- as(width, "List")
        width <- S4Vectors:::VH_recycle(width, x, "width", "x")
        unlisted_width <- unlist(width, use.names=FALSE)
        ## Compute unlisted 'ans'.
        unlisted_ans <- flank(unlisted_x, unlisted_width,
                              start=start, both=both,
                              use.names=use.names, ignore.strand=ignore.strand)
        ## Relist and return.
        relist(unlisted_ans, x)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### promoters()
###

setMethod("promoters", "GenomicRanges",
    function(x, upstream=2000, downstream=200, ...)
    {
        if (!isSingleNumber(upstream))
            stop("'upstream' must be a single integer")
        if (!is.integer(upstream))
            upstream <- as.numeric(upstream)
        if (!isSingleNumber(downstream))
            stop("'downstream' must be a single integer")
        if (!is.integer(downstream))
            downstream <- as.numeric(downstream)
        if (upstream < 0 | downstream < 0)
            stop("'upstream' and 'downstream' must be integers >= 0")
        if (any(strand(x) == "*"))
            warning("'*' ranges were treated as '+'")
        on_plus <- which(strand(x) == "+" | strand(x) == "*")
        on_plus_TSS <- start(x)[on_plus]
        start(x)[on_plus] <- on_plus_TSS - upstream
        end(x)[on_plus] <- on_plus_TSS + downstream - 1L
        on_minus <- which(strand(x) == "-")
        on_minus_TSS <- end(x)[on_minus]
        end(x)[on_minus] <- on_minus_TSS + upstream
        start(x)[on_minus] <- on_minus_TSS - downstream + 1L
        x
    }
)

setMethod("promoters", "GRangesList",
    function(x, upstream=2000, downstream=200, ...)
    {
        x@unlistData <- promoters(x@unlistData, upstream=upstream,
                                  downstream=downstream)
        x
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
        nms <- names(mcols(x))
        mcols(x) <- cbind(mcols(x), DataFrame(posIndx=seq_len(length(x))))
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
        mcols(ord) <- subset(mcols(ord), select=nms)
        ord
    }
)

setMethod("restrict", "GRangesList",
    function(x, start = NA, end = NA, keep.all.ranges = FALSE, use.names = TRUE)
    {
        endoapply(x, restrict, start=start, end=end,keep.all.ranges=keep.all.ranges
               , use.names=use.names )

    })


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
        clone(x, ranges=new_ranges)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Zooming (symmetrically scales the width).
###

setMethod("Ops", c("GenomicRanges", "numeric"),
    function(e1, e2)
    {
        if (S4Vectors:::anyMissing(e2))
            stop("NA not allowed as zoom factor")
        e2 <- recycleNumericArg(e2, "e2", length(e1))
        if (.Generic == "*") {
            e2 <- ifelse(e2 < 0, abs(1/e2), e2)
            resize(e1, width(e1) / e2, fix="center")
        } else {
            if (.Generic == "-") {
                e2 <- -e2
                .Generic <- "+"
            }
            if (.Generic == "+") {
                if (any(-e2*2 > width(e1)))
                    stop("adjustment would result in ranges ",
                         "with negative widths")
                resize(e1, width(e1) + e2*2, fix="center")
            }
        }
    }
)

