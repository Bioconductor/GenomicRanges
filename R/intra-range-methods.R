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
        update(x, ranges=new_ranges)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### narrow()
###

### The default "narrow" method calls windows() so we only need to implement
### a "windows" method for GenomicRanges objects to make narrow() work on
### these objects.
setMethod("windows", "GenomicRanges",
    function(x, start=NA, end=NA, width=NA)
    {
        ranges(x) <- windows(ranges(x), start=start, end=end, width=width)
        x
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
        ranges(x) <- new_ranges
        x
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
        new_ranges <- flank(ranges(x), width=width, start=start, both=both,
                                       use.names=use.names)
        ranges(x) <- new_ranges
        x
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### promoters()
###

setMethod("promoters", "GenomicRanges",
    function(x, upstream=2000, downstream=200)
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
        ranges(x) <- new_ranges
        x
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

