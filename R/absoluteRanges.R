### Transform the genomic ranges in 'x' into "absolute" ranges i.e. into
### ranges counted from the beginning of the virtual sequence obtained by
### concatenating all the sequences in the genome (in the order reported by
### 'seqlevels(x)'). Ignore the strand.
### IMPORTANT NOTE: The function might fail because of an integer overflow if
### the length of the genome is >= 2^31 (e.g. if it's the Human genome).
absoluteRanges <- function(x)
{
    if (!is(x, "GenomicRanges"))
        stop("'x' must be a GenomicRanges object")
    x_seqlengths <- seqlengths(x)
    if (any(is.na(x_seqlengths)))
        stop("'seqlengths(x)' cannot contain NAs")
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

### TODO: Use this in tileGenome().
normarg_seqlengths <- function(seqlengths)
{
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
    if (S4Vectors:::anyMissingOrOutside(seqlengths, lower=0L)) 
        stop("'seqlengths' contains NAs or negative values")
    seqlengths
}

### The reverse of absoluteRanges().
relativeRanges <- function(x, seqlengths)
{
    if (!is(x, "Ranges"))
        stop("'x' must be a Ranges object")
    if (is(seqlengths, "Seqinfo")) {
        ans_seqinfo <- seqlengths
        seqlengths <- seqlengths(ans_seqinfo)
    } else {
        seqlengths <- normarg_seqlengths(seqlengths)
        ans_seqinfo <- Seqinfo(seqnames=names(seqlengths),
                               seqlengths=seqlengths)
    }
    offsets <- c(0L, cumsum(unname(seqlengths)[-length(seqlengths)]))
    chrom_starts <- offsets + 1L
    start2chrom <- findInterval(start(x), chrom_starts)
    end2chrom <- findInterval(end(x), chrom_starts)
    if (!identical(start2chrom, end2chrom))
        stop(wmsg("Some ranges in 'x' are crossing sequence boundaries. ",
                  "Cannot convert them into relative ranges."))
    ans_ranges <- shift(x, shift=-offsets[start2chrom])
    ans_seqnames <- names(seqlengths)[start2chrom]
    GRanges(ans_seqnames, ans_ranges, seqinfo=ans_seqinfo)
}

