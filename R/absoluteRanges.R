### =========================================================================
### absoluteRanges() & related
### -------------------------------------------------------------------------


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### isSmallGenome()
###

### TODO: Maybe this could be re-used in tileGenome().
normarg_seqlengths <- function(seqlengths)
{
    if (!is.numeric(seqlengths)) 
        stop(wmsg("'seqlengths' must be a non-empty numeric vector"))
    if (length(seqlengths) == 0L)
        return(setNames(integer(0), character(0)))
    seqlengths_names <- names(seqlengths)
    if (is.null(seqlengths_names)) 
        stop(wmsg("'seqlengths' must be named"))
    if (any(seqlengths_names %in% c(NA_character_, ""))) 
        stop(wmsg("'seqlengths' has names that are NA or the empty string"))
    if (any(duplicated(seqlengths_names))) 
        stop(wmsg("'seqlengths' has duplicated names"))
    if (!is.integer(seqlengths)) 
        seqlengths <- setNames(as.integer(seqlengths), seqlengths_names)
    if (any(seqlengths < 0L, na.rm=TRUE))
        stop(wmsg("'seqlengths' contains negative values"))
    seqlengths
}

### 'seqlengths' can be an integer or numeric vector, or any object from which
### the sequence lengths can be extracted with seqlengths().
### Returns TRUE if the total length of the underlying sequences is <= 
### '.Machine$integer.max' (e.g. Fly genome), FALSE if not (e.g. Human genome),
### or NA if it cannot be computed (because some sequence lengths are NA).
isSmallGenome <- function(seqlengths)
{
    if (is.numeric(seqlengths)) {
        seqlengths <- normarg_seqlengths(seqlengths)
    } else {
        seqlengths <- seqlengths(seqlengths)
    }
    if (any(is.na(seqlengths)))
        return(NA)
    total_length <- sum(seqlengths)  # no more integer overflow in R >= 3.5
    total_length <= .Machine$integer.max
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### absoluteRanges()
###

### Transform the genomic ranges in 'x' into "absolute" ranges i.e. into
### ranges counted from the beginning of the virtual sequence obtained by
### concatenating all the sequences in the genome (in the order reported by
### 'seqlevels(x)'). Ignore the strand.
### ONLY WORK ON A SMALL GENOME! (see isSmallGenome() above)
absoluteRanges <- function(x)
{
    if (!is(x, "GenomicRanges"))
        stop(wmsg("'x' must be a GenomicRanges object"))
    x_seqlengths <- seqlengths(x)
    if (!isTRUE(isSmallGenome(x_seqlengths)))
        stop(wmsg("the total length of the underlying sequences is too big ",
                  "or couldn't be computed (because some lengths are NA)"))
    x_seqids <- as.integer(seqnames(x))
    idx <- which(start(x) < 1L | end(x) > x_seqlengths[x_seqids])
    if (length(idx) != 0L)
        stop(wmsg("Some ranges in 'x' are not within the bounds of ",
                  "the sequence that they belong to. Cannot convert ",
                  "them into absolute ranges."))
    offsets <- c(0L, cumsum(unname(x_seqlengths)[-length(x_seqlengths)]))
    x_ranges <- ranges(x)
    shift(x_ranges, shift=offsets[x_seqids])
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### relativeRanges()
###

### The reverse of absoluteRanges().
### ONLY WORK ON A SMALL GENOME! (see isSmallGenome() above)
relativeRanges <- function(x, seqlengths)
{
    if (!is(x, "IntegerRanges"))
        stop(wmsg("'x' must be an IntegerRanges object"))
    if (is.numeric(seqlengths)) {
        ans_seqlengths <- normarg_seqlengths(seqlengths)
        ans_seqinfo <- Seqinfo(seqnames=names(ans_seqlengths),
                               seqlengths=ans_seqlengths)
    } else {
        if (is(seqlengths, "Seqinfo")) {
            ans_seqinfo <- seqlengths
        } else {
            ans_seqinfo <- seqinfo(seqlengths)
        }
        ans_seqlengths <- seqlengths(ans_seqinfo)
    }
    if (!isTRUE(isSmallGenome(ans_seqlengths)))
        stop(wmsg("the total length of the sequences specified ",
                  "thru 'seqlengths' is too big or couldn't be ",
                  "computed (because some lengths are NA)"))
    offsets <- c(0L, cumsum(unname(ans_seqlengths)))

    ## Map each range in 'x' to a sequence in the genome.
    ticks <- offsets + 1L
    start2seqid <- findInterval(start(x), ticks)
    end2seqid <- findInterval(end(x),  ticks)
    if (!identical(start2seqid, end2seqid))
        stop(wmsg("Some ranges in 'x' cannot be mapped to a sequence in the ",
                  "genome because they cross sequence boundaries. ",
                  "Cannot convert them into relative ranges."))
    if (any(start2seqid < 1L) || any(start2seqid > length(ans_seqlengths)))
        stop(wmsg("Some ranges in 'x' cannot be mapped to a sequence in the ",
                  "genome because they are outside the boundaries of the ",
                  "genome. Cannot convert them into relative ranges."))

    ans_ranges <- shift(x, shift=-offsets[start2seqid])
    ans_seqnames <- names(ans_seqlengths)[start2seqid]
    GRanges(ans_seqnames, ans_ranges, seqinfo=ans_seqinfo)
}

