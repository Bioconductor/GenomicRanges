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
### "update" method
###

### Having update() redirect to BiocGenerics:::replaceSlots() on GRanges
### objects makes all the methods for GenomicRanges objects defined in
### R/GenomicRanges-class.R work on GRanges objects.
setMethod("update", "GRanges",
    function(object, ...) BiocGenerics:::replaceSlots(object, ...)
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity
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
### Constructor
###

.set_strand_mcols_seqinfo <- function(x, strand=NULL, mcols=NULL,
                                         seqlengths=NULL, seqinfo=NULL)
{
    if (is.null(strand)) {
        x_strand <- strand(x)
    } else {
        x_strand <- strand
    }
    if (length(mcols) == 0L) {
        x_mcols <- mcols(x)
    } else {
        x_mcols <- mcols
    }
    if (is.null(seqlengths)) {
        x_seqlengths <- seqlengths(x)
    } else {
        x_seqlengths <- seqlengths
    }
    if (is.null(seqinfo)) {
        x_seqinfo <- seqinfo(x)
    } else {
        x_seqinfo <- seqinfo
    }
    new_GRanges(class(x), seqnames(x), ranges(x), x_strand,
                          x_mcols, x_seqlengths, x_seqinfo)
}

### Internal low-level constructor. Shared with other GRanges-like objects.
new_GRanges <- function(Class, seqnames=NULL, ranges=NULL, strand=NULL,
                               mcols=NULL, seqlengths=NULL, seqinfo=NULL)
{
    if (is.null(ranges)) {
        if (!is.null(seqnames)) {
            x <- as(seqnames, Class)
            return(.set_strand_mcols_seqinfo(x, strand, mcols,
                                                seqlengths, seqinfo))
        }
        ranges <- IRanges()
    } else {
        ranges <- as(ranges, "IRanges")
    }

    if (is.null(seqnames)) {
        seqnames <- Rle()
    } else {
        if (!is(seqnames, "Rle"))
            seqnames <- Rle(seqnames)
        if (!is.factor(runValue(seqnames))) 
            runValue(seqnames) <- factor(runValue(seqnames),
                                         levels=unique(runValue(seqnames)))
    }

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

    ranges_mcols <- mcols(ranges)
    if (!is.null(ranges_mcols))
        mcols(ranges) <- NULL

    ## Normalize 'mcols'.
    if (length(mcols) == 0L) {
        mcols <- ranges_mcols
        if (is.null(mcols))
            mcols <- DataFrame()
    } else if (!is(mcols, "DataFrame")) {
        stop("'mcols' must be a DataFrame object")
    }
    if (nrow(mcols) == 0L && ncol(mcols) == 0L) {
        mcols <- S4Vectors:::make_zero_col_DataFrame(length(ranges))
    } else if (!is.null(rownames(mcols))) {
        if (is.null(names(ranges)))
            names(ranges) <- rownames(mcols)
        rownames(mcols) <- NULL
    }

    new(Class, seqnames=seqnames, ranges=ranges, strand=strand,
               elementMetadata=mcols, seqinfo=seqinfo)
}

GRanges <- function(seqnames=NULL, ranges=NULL, strand=NULL,
                    ..., seqlengths=NULL, seqinfo=NULL)
{
    mcols <- DataFrame(..., check.names=FALSE)
    new_GRanges("GRanges", seqnames=seqnames, ranges=ranges, strand=strand,
                           mcols=mcols, seqlengths=seqlengths, seqinfo=seqinfo)
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
### Coercion
###

setMethod("granges", "GenomicRanges",
    function(x, use.names=TRUE, use.mcols=FALSE)
    {
        if (!isTRUEorFALSE(use.mcols))
            stop("'use.mcols' must be TRUE or FALSE")
        ans <- GRanges(seqnames(x),
                       ranges(x, use.names=use.names),
                       strand(x),
                       seqinfo=seqinfo(x))
        if (use.mcols)
            mcols(ans) <- cbind(extraColumnSlotsAsDF(x), mcols(x))
        ans
    }
)

setAs("GenomicRanges", "GRanges",
    function(from) granges(from, use.mcols=TRUE)
)

.from_character_to_GRanges <- function(from)
{
    stopifnot(is.character(from))
    if (anyNA(from))
        stop(wmsg("converting a character vector to a GRanges object ",
                  "does not support NAs"))
    error_msg <- wmsg(
        "The character vector to convert to a GRanges object must contain ",
        "strings of the form \"chr:start-end\" or \"chr:start-end:strand\", ",
        "with end >= start - 1, or \"chr:pos\" or \"chr:pos:strand\". ",
        "For example: \"chr1:2501-2900\", \"chr1:2501-2900:+\", or ",
        "\"chr1:740\". Note that \"..\" is a valid alternate start/end ",
        "separator. Strand can be \"+\", \"-\", \"*\", or missing."
    )
    split0 <- CharacterList(strsplit(from, ":", fixed=TRUE))
    split0_eltNROWS <- elementNROWS(split0)
    if (S4Vectors:::anyMissingOrOutside(split0_eltNROWS, 2L, 3L))
        stop(error_msg)
    ans_strand <- as.character(ptail(split0, n=-2L))
    ans_strand[is.na(ans_strand)] <- "*"
    split1 <- phead(split0, n=2L)
    ans_seqnames <- as.character(phead(split1, n=1L))
    ranges <- ptail(split1, n=-1L)
    ranges <- setNames(as.character(ranges), names(ranges))
    ans_ranges <- try(as(ranges, "IRanges"), silent=TRUE)
    if (is(ans_ranges, "try-error"))
        stop(error_msg)
    GRanges(ans_seqnames, ans_ranges, ans_strand)
}
setAs("character", "GRanges", .from_character_to_GRanges)

.from_factor_to_GRanges <- function(from)
{
    from <- setNames(as.character(from), names(from))
    .from_character_to_GRanges(from)
}
setAs("factor", "GRanges", .from_factor_to_GRanges)

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

.from_Seqinfo_to_GRanges <- function(from)
{
    if (anyNA(seqlengths(from)))
        stop(wmsg("cannot create a GRanges object ",
                  "from a Seqinfo object with NA seqlengths"))
    GRanges(seqnames(from),
            IRanges(rep.int(1L, length(from)),
                    width=seqlengths(from),
                    names=seqnames(from)),
            seqinfo=from)
}

setAs("Seqinfo", "GRanges", .from_Seqinfo_to_GRanges)

setAs("Seqinfo", "RangesList",
    function(from) as(as(from, "GRanges"), "RangesList")
)

setAs("ANY", "GenomicRanges", function(from) as(from, "GRanges"))


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Slot getters and setters
###

setMethod("seqnames", "GRanges", function(x) x@seqnames)
setMethod("ranges", "GRanges",
    function(x, use.names=TRUE, use.mcols=FALSE)
    {
        if (!isTRUEorFALSE(use.names))
            stop("'use.names' must be TRUE or FALSE")
        if (!isTRUEorFALSE(use.mcols))
            stop("'use.mcols' must be TRUE or FALSE")
        ans <- x@ranges
        if (!use.names)
            names(ans) <- NULL
        if (use.mcols)
            mcols(ans) <- mcols(x)
        ans
    }
)
setMethod("strand", "GRanges", function(x) x@strand)
setMethod("seqinfo", "GRanges", function(x) x@seqinfo)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Subsetting
###

.DollarNames.GRanges <- .DollarNames.GenomicRanges

setMethod("extractROWS", "GRanges",
    function(x, i)
    {
        if (missing(i) || !is(i, "Ranges"))
            i <- normalizeSingleBracketSubscript(i, x)
        ans_seqnames <- extractROWS(seqnames(x), i)
        ans_ranges <- extractROWS(ranges(x), i)
        ans_strand <- extractROWS(strand(x), i)
        ans_mcols <- extractROWS(mcols(x), i)
        ans_ecs <- lapply(extraColumnSlots(x), extractROWS, i)
        BiocGenerics:::replaceSlots(x, seqnames=ans_seqnames,
                                       ranges=ans_ranges,
                                       strand=ans_strand,
                                       elementMetadata=ans_mcols,
                                       .slotList=ans_ecs,
                                       check=FALSE)
    }
)

setMethod("replaceROWS", "GRanges",
    function(x, i, value)
    {
        if (missing(i) || !is(i, "Ranges"))
            i <- normalizeSingleBracketSubscript(i, x)
        seqinfo(x) <- merge(seqinfo(x), seqinfo(value))
        ans_seqnames <- replaceROWS(seqnames(x), i, seqnames(value))
        ans_ranges <- replaceROWS(ranges(x), i, ranges(value))
        ans_strand <- replaceROWS(strand(x), i, strand(value))
        ans_mcols <- replaceROWS(mcols(x), i, mcols(value))
        ans_ecs_names <- extraColumnSlotNames(x)
        ans_necs <- length(ans_ecs_names)
        if (ans_necs == 0L) {
            ans_ecs <- NULL
        } else {
            value_ecs_names <- extraColumnSlotNames(value)
            if (!identical(value_ecs_names[seq_len(ans_necs)],
                           ans_ecs_names))
                stop("'value' can have more extra column slots but not less")
            ans_ecs <- extraColumnSlotsAsDF(x)
            value_ecs <- extraColumnSlotsAsDF(value)
            ans_ecs <- replaceROWS(ans_ecs, i, value_ecs[seq_len(ans_necs)])
        }
        BiocGenerics:::replaceSlots(x, seqnames=ans_seqnames,
                                       ranges=ans_ranges,
                                       strand=ans_strand,
                                       elementMetadata=ans_mcols,
                                       .slotList=as.list(ans_ecs))
    }
)

### TODO: Refactor to use replaceROWS(). This will make the code much simpler
### and avoid a lot of duplication with the above "replaceROWS" method.
setReplaceMethod("[", "GRanges",
    function(x, i, j, ..., value)
    {
        if (!is(value, "GenomicRanges"))
            stop("replacement value must be a GenomicRanges object")
        seqinfo(x) <- merge(seqinfo(x), seqinfo(value))
        seqnames <- seqnames(x)
        ranges <- ranges(x)
        strand <- strand(x)
        ans_mcols <- mcols(x, FALSE)
        value_ecs <- extraColumnSlotsAsDF(value)
        x_ecs <- extraColumnSlotsAsDF(x)
        new_ecs <- value_ecs[!names(value_ecs) %in% names(x_ecs)]
        ecs_to_replace <- intersect(names(value_ecs), names(x_ecs))
        if (missing(i)) {
            seqnames[] <- seqnames(value)
            ranges[] <- ranges(value)
            strand[] <- strand(value)
            if (missing(j))
                ans_mcols[ , ] <- mcols(value, FALSE)
            else
                ans_mcols[ , j] <- mcols(value, FALSE)
            if (length(new_ecs) > 0L)
                ans_mcols[names(new_ecs)] <- new_ecs
            x_ecs[ecs_to_replace] <- value_ecs[ecs_to_replace]
        } else {
            i <- extractROWS(setNames(seq_along(x), names(x)), i)
            seqnames[i] <- seqnames(value)
            ranges[i] <- ranges(value)
            strand[i] <- strand(value)
            if (missing(j))
                ans_mcols[i, ] <- mcols(value, FALSE)
            else
                ans_mcols[i, j] <- mcols(value, FALSE)
            if (length(new_ecs) > 0L)
                ans_mcols[i, names(new_ecs)] <- DataFrame(new_ecs)
            if (length(ecs_to_replace) > 0L) {
              x_ecs[i, ecs_to_replace] <- value_ecs[ecs_to_replace]
            }
        }
        BiocGenerics:::replaceSlots(x, seqnames=seqnames,
                                       ranges=ranges,
                                       strand=strand,
                                       elementMetadata=ans_mcols,
                                       .slotList=as.list(x_ecs))
    }
)

