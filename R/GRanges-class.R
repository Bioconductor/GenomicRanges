### =========================================================================
### GRanges objects
### -------------------------------------------------------------------------
###

setClass("GRanges",
    contains="GenomicRanges",
    representation(
        seqnames="Rle",
        ranges="IRanges",
        strand="Rle",
        elementMetadata="DataFrame",
        seqinfo="Seqinfo"
    ),
    prototype(
        seqnames=Rle(factor()),
        strand=Rle(strand())
    )
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity.
###

.valid.GRanges.ranges <- function(x)
{
    if (!is.null(x@ranges@elementMetadata))
        return("slot 'ranges' cannot have metadata columns")
    NULL
}

.valid.GRanges.mcols <- function(x)
{
    x_mcols <- x@elementMetadata
    if (!is.null(rownames(x_mcols)))
        return("'mcols(x)' cannot have row names")
    NULL
}

.valid.GRanges <- function(x)
{
    c(.valid.GRanges.ranges(x), .valid.GRanges.mcols(x))
}

setValidity2("GRanges", .valid.GRanges)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor.
###

### Hidden constructor shared with other GRanges-like objects.
newGRanges <- function(class,
                       seqnames=Rle(), ranges=NULL, strand=NULL,
                       mcols=DataFrame(),
                       seqlengths=NULL, seqinfo=NULL)
{
    if (is.null(ranges)) {
        if (length(seqnames) != 0L) {
            ans <- as(seqnames, class)
            ans_seqnames <- seqnames(ans)
            ans_ranges <- ranges(ans)
            if (is.null(strand)) {
                ans_strand <- strand(ans)
            } else {
                ans_strand <- strand
            }
            if (ncol(mcols) == 0L) {
                ans_mcols <- mcols(ans)
            } else {
                ans_mcols <- mcols
            }
            if (is.null(seqlengths)) {
                ans_seqlengths <- seqlengths(ans)
            } else {
                ans_seqlengths <- seqlengths
            }
            if (is.null(seqinfo)) {
                ans_seqinfo <- seqinfo(ans)
            } else {
                ans_seqinfo <- seqinfo
            }
            ans <- newGRanges(class,
                              ans_seqnames, ans_ranges, ans_strand,
                              ans_mcols,
                              ans_seqlengths, ans_seqinfo)
            return(ans)
        }
        ranges <- IRanges()
    } else if (class(ranges) != "IRanges") {
        ranges <- as(ranges, "IRanges")
    }

    if (!is(seqnames, "Rle"))
        seqnames <- Rle(seqnames)
    if (!is.factor(runValue(seqnames))) 
        runValue(seqnames) <- factor(runValue(seqnames),
                                     levels=unique(runValue(seqnames)))

    if (is.null(strand)) {
        strand <- Rle(strand("*"), length(seqnames))
    } else {
        if (!is(strand, "Rle"))
            strand <- Rle(strand)
        if (!is.factor(runValue(strand)) ||
            !identical(levels(runValue(strand)), levels(strand())))
            runValue(strand) <- strand(runValue(strand))
        if (S4Vectors:::anyMissing(runValue(strand))) {
            warning("missing values in strand converted to \"*\"")
            runValue(strand)[is.na(runValue(strand))] <- "*"
        }
    }

    lx <- max(length(seqnames), length(ranges), length(strand))
    if (lx > 1) {
        if (length(seqnames) == 1)
            seqnames <- rep(seqnames, lx)
        if (length(ranges) == 1)
            ranges <- rep(ranges, lx)
        if (length(strand) == 1)
            strand <- rep(strand, lx)
    }

    if (is.null(seqlengths))
        seqlengths <- setNames(rep(NA_integer_, length(levels(seqnames))),
                               levels(seqnames))

    if (is.null(seqinfo))
        seqinfo <- Seqinfo(names(seqlengths), seqlengths)

    ## in case we have seqlengths for unrepresented sequences
    runValue(seqnames) <- factor(runValue(seqnames), levels=seqnames(seqinfo))

    if (!is(mcols, "DataFrame"))  # should never happen when calling GRanges()
        stop("'mcols' must be a DataFrame object")
    if (ncol(mcols) == 0L) {
        mcols <- mcols(ranges)
        if (is.null(mcols))
            mcols <- new("DataFrame", nrows=length(seqnames))
    }
    if (!is.null(mcols(ranges)))
        mcols(ranges) <- NULL
    if (!is.null(rownames(mcols))) {
        if (is.null(names(ranges)))
            names(ranges) <- rownames(mcols)
        rownames(mcols) <- NULL
    }
    new(class, seqnames = seqnames, ranges = ranges, strand = strand,
        seqinfo = seqinfo, elementMetadata = mcols)
}

GRanges <- function(seqnames=Rle(), ranges=NULL, strand=NULL,
                    ...,
                    seqlengths=NULL, seqinfo=NULL)
{
    newGRanges("GRanges", seqnames=seqnames, ranges=ranges, strand=strand,
                          mcols=DataFrame(..., check.names=FALSE),
                          seqlengths=seqlengths, seqinfo=seqinfo)
}

setMethod("updateObject", "GRanges",
    function(object, ..., verbose=FALSE)
    {
        if (verbose)
            message("updateObject(object = 'GRanges')")
        if (is(try(object@seqinfo, silent=TRUE), "try-error")) {
            object <- new(class(object),
                          seqnames = object@seqnames,
                          ranges = object@ranges,
                          strand = object@strand,
                          elementMetadata = object@elementMetadata,
                          metadata = object@metadata,
                          seqinfo = Seqinfo(seqnames = names(object@seqlengths),
                                            seqlengths = object@seqlengths))
            return(object)
        }
        if (is(try(validObject(object@seqinfo, complete=TRUE), silent=TRUE),
               "try-error")) {
            object@seqinfo <- updateObject(object@seqinfo)
            return(object)
        }
        object
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercion.
###

.from_character_to_GRanges <- function(from)
{
    stopifnot(is.character(from))
    if (anyNA(from))
        stop(wmsg("converting a character vector to a GRanges object ",
                  "does not support NAs"))
    error_msg <- wmsg(
        "The character vector to convert to a GRanges object must contain ",
        "strings of the form \"chr1:2501-2800\" or \"chr1:2501-2800:+\" ",
        "(\"..\" being also supported as a separator between the start and ",
        "end positions). Strand can be \"+\", \"-\", \"*\", or missing."
    )
    split0 <- CharacterList(strsplit(from, ":", fixed=TRUE))
    split0_eltlens <- elementLengths(split0)
    if (S4Vectors:::anyMissingOrOutside(split0_eltlens, 2L, 3L))
        stop(error_msg)
    ans_strand <- as.character(ptail(split0, n=-2L))
    ans_strand[is.na(ans_strand)] <- "*"
    split1 <- phead(split0, n=2L)
    ans_seqnames <- as.character(phead(split1, n=1L))
    ranges <- as.character(ptail(split1, n=-1L))
    ## We want to split on the first occurence of  "-" that is preceeded by
    ## a digit (ignoring and removing the spaces in between if any).
    ranges <- sub("([[:digit:]])[[:space:]]*-", "\\1..", ranges)
    split2 <- CharacterList(strsplit(ranges, "..", fixed=TRUE))
    split2_eltlens <- elementLengths(split2)
    if (!all(split2_eltlens == 2L))
        stop(error_msg)
    ans_start <- as.integer(phead(split2, n=1L))
    ans_end <- as.integer(ptail(split2, n=1L))
    ans_ranges <- IRanges(ans_start, ans_end, names=names(from))
    GRanges(ans_seqnames, ans_ranges, ans_strand)
}

setAs("character", "GRanges", .from_character_to_GRanges)
setAs("character", "GenomicRanges", .from_character_to_GRanges)

.from_factor_to_GRanges <- function(from)
{
    from <- setNames(as.character(from), names(from))
    .from_character_to_GRanges(from)
}

setAs("factor", "GRanges", .from_factor_to_GRanges)
setAs("factor", "GenomicRanges", .from_factor_to_GRanges)

### Does NOT propagate the ranges names and metadata columns i.e. always
### returns an unnamed GRanges object with no metadata columns.
setAs("RangesList", "GRanges",
      function(from)
      {
        if (!length(from))
          return(GRanges())
        from <- as(from, "CompressedIRangesList")
        ranges <- unlist(from, use.names=FALSE)
        ranges <- IRanges(start=start(ranges), width=width(ranges))
        ## From now, ranges is guaranteed to be an IRanges *instance*.
        if (is.null(space(from))) {
          stop("Cannot create GRanges when 'space(from)' is NULL")
        }
        gr <- GRanges(seqnames = space(from),
                      ranges = ranges,
                      strand = Rle("*", length(ranges)))
        seqinfo(gr) <- seqinfo(from)
        metadata(gr) <- metadata(from)
        gr
      })

setAs("RangedData", "GRanges",
    function(from)
    {
        ans_ranges <- unlist(ranges(from), use.names=FALSE)
        ans_mcols <- unlist(values(from), use.names=FALSE)
        rownames(ans_mcols) <- NULL
        whichStrand <- match("strand", colnames(ans_mcols))
        if (is.na(whichStrand)) {
            ans_strand <- Rle(strand("*"), length(ans_ranges))
        } else {
            ans_strand <- Rle(strand(from))
            ans_mcols <- ans_mcols[-whichStrand]
        }
        ans <- GRanges(seqnames=space(from),
                       ranges=ans_ranges,
                       strand=ans_strand,
                       ans_mcols,
                       seqinfo=seqinfo(from))
        metadata(ans) <- metadata(from)
        ans
    }
)

setAs("Seqinfo", "GRanges", .fromSeqinfoToGRanges)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Slot getters and setters.
###

setMethod("seqnames", "GRanges", function(x) x@seqnames)
setMethod("ranges", "GRanges",
    function(x, use.mcols=FALSE)
    {
        if (!isTRUEorFALSE(use.mcols))
            stop("'use.mcols' must be TRUE or FALSE")
        ans <- x@ranges
        if (use.mcols)
            mcols(ans) <- mcols(x)
        ans
    }
)
setMethod("strand", "GRanges", function(x) x@strand)
setMethod("seqinfo", "GRanges", function(x) x@seqinfo)

