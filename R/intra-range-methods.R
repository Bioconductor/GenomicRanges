### =========================================================================
### Intra-range methods
### -------------------------------------------------------------------------
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### shift()
###

setMethod("shift", "GenomicRanges",
    function(x, shift=0L, use.names=TRUE)
    {
        ranges <- shift(ranges(x), shift, use.names=use.names)

        withCallingHandlers({
            clone(x, ranges=ranges) 
        }, warning=function(warn) {
            msg <- conditionMessage(warn)
            exp <- gettext("'ranges' contains values outside of sequence bounds",
                           domain="R")
            if (msg == exp) {
                msg <- paste0(msg, ". See ?trim to subset ranges within",
                              " sequence bounds.")
                warning(simpleWarning(msg, conditionCall(warn)))
                invokeRestart("muffleWarning")
            } else {
                warn
            }
        })
    }
) 


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### narrow()
###

setMethod("narrow", "GenomicRanges",
    function(x, start=NA, end=NA, width=NA, use.names=TRUE)
    {
        rng <- narrow(ranges(x), start=start, end=end, width=width,
                      use.names=TRUE)
        ranges(x) <- rng
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
        if (ignore.strand)
            start <- rep.int(TRUE, length(x))
        else 
            start <- as.vector(start == (strand(x) != "-"))
        ranges <-
            flank(ranges(x), width=width, start=start, both=both,
                  use.names=use.names)
        if (!IRanges:::anyMissing(seqlengths(x))) {
            start(x) <- start(ranges)
            end(x) <- end(ranges)
        } else {
            x <- clone(x, ranges=ranges)
        }
        x
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### reflect()
###

### TODO: Add "reflect" method for GenomicRanges objects.


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### resize()
###

setMethod("resize", "GenomicRanges",
    function(x, width, fix="start", use.names=TRUE, ignore.strand=FALSE)
    {
        if (!missing(fix) &&
            (length(fix) > length(x) || length(x) %% length(fix) > 0L))
            stop("'x' is not a multiple of 'fix' length")
        if (!isTRUEorFALSE(ignore.strand))
            stop("'ignore.strand' must be TRUE or FALSE")
        if (ignore.strand) {
            fix <- Rle(rep.int(fix, length(x)))
        } else {
            revFix <- c(start="end", end="start", center="center")       
            fix <- ifelse(strand(x) == "-", revFix[fix], fix)
        }
        ranges <-
            resize(ranges(x), width=width, fix=fix, use.names=use.names)
        if (!IRanges:::anyMissing(seqlengths(x))) {
            start(x) <- start(ranges)
            end(x) <- end(ranges)
        } else {
            x <- clone(x, ranges=ranges)
        }
        x
    }
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### trim()
###

setMethod("trim", "GenomicRanges",
    function(x, use.names=TRUE)
    {
        end <- NA_integer_ 
        if (any(!is.na(seqlengths(x)))) {
            seqlen <- seqlengths(x)
            seqlen[isCircular(x) %in% TRUE] <- NA_integer_
            idx <- match(as.character(seqnames(x)), names(seqlen))
            end <- seqlen[idx]
        }
        x@ranges <- restrict(ranges(x), start=1L, end=end, 
                             keep.all.ranges=TRUE, use.names=use.names)
        x
    }
)


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
        df <- data.frame(orig=mcols(ord)$posIndx,
                         final=seq_len(length(ord)))
        indx <- with(df, order(orig, final))
        ord <- ord[indx, ]
        mcols(ord) <- subset(mcols(ord), select=-c(posIndx))
        ord
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Zooming (symmetrically scales the width).
###

setMethod("Ops", c("GenomicRanges", "numeric"),
    function(e1, e2)
    {
        if (IRanges:::anyMissing(e2))
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

