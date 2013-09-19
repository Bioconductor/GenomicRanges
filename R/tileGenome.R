### =========================================================================
### tileGenome()
### -------------------------------------------------------------------------


.make_breakpoints_for_cut_last_tile_in_chrom <- function(seqlengths, tilesize)
{
    tile_relative_breakpoints <- lapply(seqlengths,
        function(seqlength) {
            nfulltile <- seqlength %/% tilesize
            if (nfulltile == 0L)
                return(seqlength)
            relative_breakpoints <- ceiling(tilesize * seq_len(nfulltile))
            if (relative_breakpoints[[nfulltile]] >= seqlength)
                return(relative_breakpoints)
            c(relative_breakpoints, seqlength)
        })
    chrom_breakpoints <- cumsum(as.numeric(seqlengths))
    chrom_offsets <- c(0, head(chrom_breakpoints, n=-1L))
    ntile_per_chrom <- elementLengths(tile_relative_breakpoints)
    unlist(tile_relative_breakpoints, use.names=FALSE) +
        rep.int(chrom_offsets, ntile_per_chrom)
}

### 'old_breakpoints' and 'new_breakpoints' must be non-negative non-decreasing
### numeric or integer vectors of length >= 1 with no NAs. In addition,
### 'old_breakpoints' must be named. None of this is checked (we trust the
### caller).
### Returns a NumericList (or IntegerList) object with one list element per
### unique new breakpoint.
.superimpose_breakpoints <- function(old_breakpoints, new_breakpoints)
{
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

    ## Compute and return final result.
    unname(splitAsList(all_breakpoints, all2new))
}

.get_relative_ends <- function(absolute_ends, chrom_breakpoints)
{
    chrom_offsets <- c(0, head(chrom_breakpoints, n=-1L))
    names(chrom_offsets) <- names(chrom_breakpoints)
    absolute_ends <- unlist(absolute_ends, use.names=FALSE)
    absolute_ends - chrom_offsets[names(absolute_ends)]
}

.get_relative_starts <- function(relative_ends, seqnames)
{
    relative_starts <- rep.int(1L, length(relative_ends))
    run_lens <- runLength(seqnames)
    run_starts <- c(1L, cumsum(head(run_lens, n=-1L)) + 1L)
    idx <- IRanges:::fancy_mseq(run_lens - 1L, offset=run_starts)
    relative_starts[idx] <- relative_ends[idx - 1L] + 1L
    relative_starts
}

### 'seqlengths' must be a non-negative numeric or integer vector of length
### >= 1 with no NAs and with unique names.
### Only one of 'ntile' and 'tilesize' can be specified.
### If 'cut.last.tile.in.chrom' is FALSE (the default), returns a GRangesList
### object with 1 list element per tile, each of them containing a number of
### ranges equal to the number of chromosomes spanned by the tile (when tiles
### are small, most of them only span 1 chromosome).
### If 'cut.last.tile.in.chrom' is TRUE (only supported when 'tilesize' is
### specified), tiles never span more than 1 chromosome (the last tile in each
### chromosome is cut if necessary). In that case, a GRanges object is
### returned with 1 element per tile. 
tileGenome <- function(seqlengths, ntile, tilesize,
                       cut.last.tile.in.chrom=FALSE)
{
    if (!isTRUEorFALSE(cut.last.tile.in.chrom))
        stop("'cut.last.tile.in.chrom' must be TRUE or FALSE")

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

    chrom_breakpoints <- setNames(cumsum(as.numeric(seqlengths)),
                                  seqlengths_names)
    genome_size <- chrom_breakpoints[[length(chrom_breakpoints)]]

    if (!missing(ntile)) {
        if (!missing(tilesize))
            stop("only one of 'ntile' and 'tilesize' can be specified")
        if (cut.last.tile.in.chrom)
            stop("'cut.last.tile.in.chrom' must be FALSE ",
                 "when 'ntile' is supplied")
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
    } else {
        if (missing(tilesize))
            stop("one of 'ntile' and 'tilesize' must be specified")
        ## Check 'tilesize'.
        if (!isSingleNumber(tilesize))
            stop("'tilesize' must be a single number")
        if (tilesize < 1 || tilesize > genome_size)
            stop("'tilesize' must be >= 1 and <= genome size")
        if (cut.last.tile.in.chrom) {
            tile_breakpoints <-
                .make_breakpoints_for_cut_last_tile_in_chrom(seqlengths,
                                                             tilesize)
        } else {
            ntile <- ceiling(genome_size / tilesize)
            if (ntile > .Machine$integer.max)
                stop("this 'tilesize' is too small")
            ntile <- as.integer(ntile)
        }
    }
    if (!cut.last.tile.in.chrom) {
        tilesize <- genome_size / ntile
        tile_breakpoints <- ceiling(tilesize * seq_len(ntile))
    }
    absolute_ends <- .superimpose_breakpoints(chrom_breakpoints,
                                              tile_breakpoints)
    gr_end <- .get_relative_ends(absolute_ends, chrom_breakpoints)
    gr_seqnames <- Rle(factor(names(gr_end), levels=names(seqlengths)))
    gr_start <- .get_relative_starts(gr_end, gr_seqnames)
    gr <- GRanges(gr_seqnames, IRanges(gr_start, gr_end),
                  seqlengths=seqlengths)
    if (cut.last.tile.in.chrom)
        return(gr)
    relist(gr, absolute_ends)
}

