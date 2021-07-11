### =========================================================================
### GRanges objects
### -------------------------------------------------------------------------
###


setClassUnion("IRanges_OR_IPos", c("IRanges", "IPos"))

setClass("GRanges",
    contains="GenomicRanges",
    representation(
        seqnames="Rle",
        ranges="IRanges_OR_IPos",  # an IPos only for GPos
        strand="Rle",
        seqinfo="Seqinfo"
    ),
    prototype(
        seqnames=Rle(factor()),
        strand=Rle(strand())
    )
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### parallel_slot_names()
###

### Combine the new "parallel slots" with those of the parent class. Make
### sure to put the new parallel slots **first**. See R/Vector-class.R file
### in the S4Vectors package for what slots should or should not be considered
### "parallel".
setMethod("parallel_slot_names", "GRanges",
    #function(x) c("seqnames", "ranges", "strand", callNextMethod())

    ## TEMPORARY DEFINITION.
    ## TODO: Remove this temporary definition and add parallel_slot_names()
    ## methods to all packages that define extraColumnSlotNames() methods
    ## (e.g. VariantAnnotation, GenomicTuples, InteractionSet, SGSeq).
    function(x) c(extraColumnSlotNames(x), "seqnames", "ranges", "strand",
                  callNextMethod())
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
### updateObject()
###
### Internal representation of GRanges objects has changed in GenomicRanges
### 1.31.16 (Bioc 3.7).
###

.get_GRanges_version <- function(object)
{
    if (.hasSlot(object, "elementType")) "current" else "< 1.31.16"
}

setMethod("updateObject", "GRanges",
    function(object, ..., verbose=FALSE)
    {
        ## elementType slot.
        version <- .get_GRanges_version(object)
        if (version == "current") {
            if (verbose)
                message("[updateObject] Internal representation of ",
                        class(object), " object is current.\n",
                        "[updateObject] Nothing to update.")
        } else {
            if (verbose)
                message("[updateObject] ", class(object), " object ",
                        "uses internal representation from\n",
                        "[updateObject] GenomicRanges ", version, ". ",
                        "Updating it ... ", appendLF=FALSE)
            object@elementType <- new(class(object))@elementType
            if (verbose)
                message("OK")
        }

        ## ranges slot.
        object@ranges <- updateObject(object@ranges, ..., verbose=verbose)

        callNextMethod()
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### "update" method
###

### Having update() redirect to BiocGenerics:::replaceSlots() on GRanges
### objects makes all the methods for GenomicRanges objects defined in
### R/GenomicRanges-class.R work on GRanges objects.
setMethod("update", "GRanges",
    function(object, ...)
    {
        ## Fix old GRanges instances on-the-fly.
        object <- updateObject(object, check=FALSE)
        BiocGenerics:::replaceSlots(object, ...)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor
###

### Return a factor-Rle with no NAs.
.normarg_seqnames1 <- function(seqnames)
{
    if (is.null(seqnames))
        return(Rle(factor()))
    if (!is(seqnames, "Rle"))
        seqnames <- Rle(seqnames)
    run_vals <- runValue(seqnames)
    if (anyNA(run_vals))
        stop(wmsg("'seqnames' cannot contain NAs"))
    if (!is.factor(run_vals)) {
        if (!is.character(run_vals))
            run_vals <- as.character(run_vals)
        runValue(seqnames) <- factor(run_vals, levels=unique(run_vals))
    }
    seqnames
}

### 'seqnames' is assumed to be a factor-Rle with no NAs (which should
### be the case if it went thru .normarg_seqnames1()).
### 'seqinfo' is assumned to be a Seqinfo object.
.normarg_seqnames2 <- function(seqnames, seqinfo)
{
    ans_seqlevels <- seqlevels(seqinfo)
    run_vals <- runValue(seqnames)
    seqnames_levels <- levels(run_vals)
    is_used <- tabulate(run_vals, nbins=length(seqnames_levels)) != 0L
    seqnames_levels_in_use <- seqnames_levels[is_used]
    if (!all(seqnames_levels_in_use %in% ans_seqlevels))
        stop(wmsg("'seqnames' contains sequence names ",
                  "with no entries in 'seqinfo'"))
    if (!all(seqnames_levels %in% ans_seqlevels))
        warning(wmsg("levels in 'seqnames' with no entries ",
                     "in 'seqinfo' were dropped"))
    runValue(seqnames) <- factor(run_vals, levels=ans_seqlevels)
    seqnames
}

### Return a factor-Rle with levels +|-|* and no NAs.
.normarg_strand <- function(strand, seqnames)
{
    if (is.null(strand))
        return(Rle(strand("*"), length(seqnames)))
    if (!is(strand, "Rle"))
        strand <- Rle(strand)
    run_vals <- runValue(strand)
    if (anyNA(run_vals)) {
        warning(wmsg("missing values in 'strand' converted to \"*\""))
        run_vals[is.na(run_vals)] <- "*"
    }
    if (!is.factor(run_vals) || !identical(levels(run_vals), levels(strand())))
        run_vals <- strand(run_vals)
    runValue(strand) <- run_vals
    strand
}

### Internal low-level constructor. Used by high-level GRanges/GPos
### constructors. Not meant to be used directly by the end user.
### NOTE: 'ranges' is trusted! (should have been checked by the caller).
new_GRanges <- function(Class, seqnames=NULL, ranges=NULL, strand=NULL,
                               mcols=NULL, seqinfo=NULL)
{
    seqnames <- .normarg_seqnames1(seqnames)

    if (is.null(seqinfo)) {
        seqinfo <- Seqinfo(levels(seqnames))
    } else {
        seqinfo <- normarg_seqinfo1(seqinfo)
        seqnames <- .normarg_seqnames2(seqnames, seqinfo)
    }

    strand <- .normarg_strand(strand, seqnames)

    seqnames_len <- length(seqnames)
    ranges_len <- length(ranges)
    strand_len <- length(strand)
    ans_len <- max(seqnames_len, ranges_len, strand_len)
    if (ans_len != 0L) {
        stop_if_wrong_length <- function(what, ans_len)
            stop(wmsg(what, " must have the length of the object ",
                      "to construct (", ans_len, ") or length 1"))
        if (seqnames_len == 0L)
            stop_if_wrong_length("'seqnames'", ans_len)
        if (ranges_len == 0L) {
            what <- if (is(ranges, "IPos")) "'pos'" else "'ranges'"
            stop_if_wrong_length(what, ans_len)
        }
        if (strand_len == 0L)
            stop_if_wrong_length("'strand'", ans_len)
        if (ans_len > 1L) {
            if (seqnames_len == 1L)
                seqnames <- rep(seqnames, ans_len)
            else if (seqnames_len != ans_len)
                stop_if_wrong_length("'seqnames'", ans_len)
            if (ranges_len == 1L)
                ranges <- rep(ranges, ans_len)
            else if (ranges_len != ans_len) {
                what <- if (is(ranges, "IPos")) "'pos'" else "'ranges'"
                stop_if_wrong_length(what, ans_len)
            }
            if (strand_len == 1L)
                strand <- rep(strand, ans_len)
            else if (strand_len != ans_len)
                stop_if_wrong_length("'strand'", ans_len)
        }
    }

    ## From now on, 'seqnames', 'ranges', and 'strand' are guaranteed
    ## to be of length 'ans_len'.

    if (length(mcols) == 0L)
        mcols <- mcols(ranges, use.names=FALSE)
    mcols <- S4Vectors:::normarg_mcols(mcols, Class, ans_len)

    if (!is.null(mcols(ranges, use.names=FALSE)))
        mcols(ranges) <- NULL
    if (!is.null(rownames(mcols))) {
        if (is.null(names(ranges)))
            names(ranges) <- rownames(mcols)
        rownames(mcols) <- NULL
    }

    new2(Class, seqnames=seqnames, ranges=ranges, strand=strand,
                elementMetadata=mcols, seqinfo=seqinfo, check=FALSE)
}

### High-level GRanges constructor.
GRanges <- function(seqnames=NULL, ranges=NULL, strand=NULL,
                    ..., seqinfo=NULL, seqlengths=NULL)
{
    mcols <- DataFrame(..., check.names=FALSE)

    if (!is.null(ranges)) {
        ranges <- as(ranges, "IRanges")
    } else if (is.null(seqnames)) {
        ranges <- IRanges()
    } else {
        x <- as(seqnames, "GRanges")
        seqnames <- x@seqnames
        ranges <- x@ranges
        if (is.null(strand))
            strand <- x@strand
        if (length(mcols) == 0L)
            mcols <- mcols(x, use.names=FALSE)
        if (is.null(seqinfo))
            seqinfo <- seqinfo(x)
    }

    seqinfo <- normarg_seqinfo2(seqinfo, seqlengths)

    ans <- new_GRanges("GRanges", seqnames=seqnames,
                                  ranges=ranges, strand=strand,
                                  mcols=mcols, seqinfo=seqinfo)
    ## new_GRanges() doesn't check validity so we do it here. Note that 'ans'
    ## should be valid except for the silly INVALID.GR.COLNAMES restriction.
    ## If it wasn't for this restriction, we actually wouldn't need to
    ## validate 'ans'.
    validObject(ans)
    ans
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Accessors
###

setMethod("seqnames", "GRanges", function(x) x@seqnames)

setMethod("strand", "GRanges", function(x) x@strand)

setMethod("seqinfo", "GRanges", function(x) x@seqinfo)

### Range squeezer.
setMethod("ranges", "GRanges",
    function(x, use.names=TRUE, use.mcols=FALSE)
    {
        if (!isTRUEorFALSE(use.names))
            stop(wmsg("'use.names' must be TRUE or FALSE"))
        if (!isTRUEorFALSE(use.mcols))
            stop(wmsg("'use.mcols' must be TRUE or FALSE"))
        ## We call updateObject() in case 'x@ranges' is an old IPos object
        ## (see updateObject() method for IPos objects).
        ans <- updateObject(x@ranges, check=FALSE)
        if (!use.names)
            names(ans) <- NULL
        if (use.mcols)
            mcols(ans) <- mcols(x, use.names=FALSE)
        ans
    }
)

### Genomic range squeezer.
setMethod("granges", "GenomicRanges",
    function(x, use.names=TRUE, use.mcols=FALSE)
    {
        if (!isTRUEorFALSE(use.mcols))
            stop(wmsg("'use.mcols' must be TRUE or FALSE"))
        ans <- GRanges(seqnames(x),
                       ranges(x, use.names=use.names),
                       strand(x),
                       seqinfo=seqinfo(x))
        if (use.mcols)
            mcols(ans) <- cbind(extraColumnSlotsAsDF(x),
                                mcols(x, use.names=FALSE))
        ans
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercion
###

setAs("GenomicRanges", "GRanges",
    function(from) granges(from, use.mcols=TRUE)
)

.from_character_to_GRanges <- function(from)
{
    stopifnot(is.character(from))
    if (anyNA(from))
        stop(wmsg("converting a character vector to a GRanges object ",
                  "does not support NAs"))
    error_msg <- c(
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
        stop(wmsg(error_msg))
    ans_strand <- as.character(tails(split0, n=-2L))
    ans_strand[is.na(ans_strand)] <- "*"
    split1 <- heads(split0, n=2L)
    ans_seqnames <- as.character(heads(split1, n=1L))
    ranges <- tails(split1, n=-1L)
    ranges <- setNames(as.character(ranges), names(ranges))
    ans_ranges <- try(as(ranges, "IRanges"), silent=TRUE)
    if (is(ans_ranges, "try-error"))
        stop(wmsg(error_msg))
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
setAs("IntegerRangesList", "GRanges",
      function(from)
      {
        if (!length(from))
          return(GRanges())
        from <- as(from, "CompressedIRangesList")
        ranges <- unlist(from, use.names=FALSE)
        ranges <- IRanges(start=start(ranges), width=width(ranges))
        ## From now, ranges is guaranteed to be an IRanges *instance*.
        if (is.null(space(from))) {
          stop(wmsg("cannot create GRanges when 'space(from)' is NULL"))
        }
        gr <- GRanges(seqnames = space(from),
                      ranges = ranges,
                      strand = Rle("*", length(ranges)),
                      seqinfo = seqinfo(from))
        metadata(gr) <- metadata(from)
        gr
      })

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

setAs("Seqinfo", "IntegerRangesList",
    function(from) as(as(from, "GRanges"), "IntegerRangesList")
)

setAs("ANY", "GenomicRanges", function(from) as(from, "GRanges"))


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Subsetting
###

.DollarNames.GRanges <- .DollarNames.GenomicRanges

setMethod("replaceROWS", "GRanges",
    function(x, i, value)
    {
        i <- normalizeSingleBracketSubscript(i, x, as.NSBS=TRUE)
        seqinfo(x) <- merge(seqinfo(x), seqinfo(value))
        ans_seqnames <- replaceROWS(seqnames(x), i, seqnames(value))
        ans_ranges <- replaceROWS(ranges(x), i, ranges(value))
        ans_strand <- replaceROWS(strand(x), i, strand(value))
        ans_mcols <- replaceROWS(mcols(x, use.names=FALSE), i,
                                 mcols(value, use.names=FALSE))
        ans_ecs_names <- extraColumnSlotNames(x)
        ans_necs <- length(ans_ecs_names)
        if (ans_necs == 0L) {
            ans_ecs <- NULL
        } else {
            value_ecs_names <- extraColumnSlotNames(value)
            if (!identical(head(value_ecs_names, n=ans_necs),
                           ans_ecs_names))
                stop(wmsg("'value' can have more extra column slots ",
                          "but not less"))
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

