### =========================================================================
### tileGenome()
### -------------------------------------------------------------------------

### 'old_breakpoints' and 'new_breakpoints' must be non-negative non-decreasing
### numeric or integer vectors of length >= 1 with no NAs. In addition,
### 'old_breakpoints' must be named. None of this is checked (we trust the
### caller).
### Returns a NumericList (or IntegerList) object with one list element per
### unique new breakpoint.
.replace_breakpoints <- function(old_breakpoints, new_breakpoints)
{
    ## Compute 'offsets'.
    offsets <- c(0L, head(old_breakpoints, n=-1L))
    names(offsets) <- names(old_breakpoints)

    ## Set names on 'new_breakpoints'.
    new2old <- 1L + findInterval(new_breakpoints - 1L, old_breakpoints)
    names(new_breakpoints) <- names(old_breakpoints)[new2old]

    ## Compute 'all_breakpoints'.
    all_breakpoints <- c(old_breakpoints, new_breakpoints)
    unique_idx <- !duplicated(all_breakpoints)
    all_breakpoints <- all_breakpoints[unique_idx]

    ## Compute 'all2new' mapping.
    old2new <- 1L + findInterval(old_breakpoints - 1L, new_breakpoints)
    all2new <- c(old2new, seq_along(new_breakpoints))
    all2new <- all2new[unique_idx]

    ## Compute 'ans' and return it.
    ans <- unname(splitAsList(all_breakpoints, all2new))
    unlisted_ans <- unlist(ans, use.names=FALSE)
    unlisted_ans <- unlisted_ans - offsets[names(unlisted_ans)]
    relist(unlisted_ans, ans)
}

### 'seqlengths' must be a non-negative numeric or integer vector of length
### >= 1 with no NAs and with unique names.
### Only one of 'ntile' and 'tilesize' can be specified.
tileGenome <- function(seqlengths, ntile, tilesize)
{
    ## Check 'seqlengths'.
    if (!is.numeric(seqlengths) || length(seqlengths) == 0L)
        stop("'seqlengths' must be a non-empty numeric vector")
    seqlengths_names <- names(seqlengths)
    if (is.null(seqlengths_names))
        stop("'seqlengths' must be named")
    if (any(seqlengths_names %in% c(NA_character_, "")))
        stop("'seqlengths' has names that are NA or the empty string")
    if (any(duplicated(seqlengths_names)))
        stop("'seqlengths' has duplicated names")
    if (!is.integer(seqlengths))
        seqlengths <- setNames(as.integer(seqlengths), seqlengths_names)
    if (IRanges:::anyMissingOrOutside(seqlengths, lower=0L))
        stop("'seqlengths' contains NAs or negative values")

    old_breakpoints <- setNames(cumsum(as.numeric(seqlengths)),
                                seqlengths_names)
    genome_size <- old_breakpoints[length(old_breakpoints)]

    if (missing(ntile)) {
        if (missing(tilesize))
            stop("one of 'ntile' and 'tilesize' must be specified")
        ## Check 'tilesize'.
        if (!isSingleNumber(tilesize))
            stop("'tilesize' must be a single number")
        if (tilesize < 1 || tilesize > genome_size)
            stop("'tilesize' must be >= 1 and <= genome size")
        ntile <- ceiling(genome_size / tilesize)
        if (ntile > .Machine$integer.max)
            stop("this 'tilesize' is too small")
        ntile <- as.integer(ntile)
    } else {
        if (!missing(tilesize))
            stop("only one of 'ntile' and 'tilesize' can be specified")
        ## Check 'ntile'.
        if (!isSingleNumber(ntile))
            stop("'ntile' must be a single number")
        if (ntile < 1 || ntile > genome_size)
            stop("'ntile' must be >= 1 and <= genome size")
        if (!is.integer(ntile)) {
            if (ntile > .Machine$integer.max)
                stop("'ntile' must be <= .Machine$integer.max")
            ntile <- as.integer(ntile)
        }
    }
    tilesize <- genome_size / ntile
    new_breakpoints <- ceiling(tilesize * seq_len(ntile))
    breakpoint_list <- .replace_breakpoints(old_breakpoints, new_breakpoints)
    gr_end <- unlist(breakpoint_list, use.names=FALSE)
    gr_seqnames <- Rle(names(gr_end))
    gr_start <- rep.int(1L, length(gr_end))
    run_lens <- runLength(gr_seqnames)
    run_starts <- c(1L, cumsum(head(run_lens, n=-1L)) + 1L)
    idx <- IRanges:::fancy_mseq(run_lens - 1L, offset=run_starts)
    gr_start[idx] <- gr_end[idx - 1L] + 1L
    gr <- GRanges(gr_seqnames, IRanges(gr_start, gr_end),
                  seqlengths=seqlengths)
    split(gr, breakpoint_list)
}

