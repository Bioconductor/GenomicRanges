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
### parallelSlotNames()
###

### Combine the new parallel slots with those of the parent class. Make sure
### to put the new parallel slots *first*.
setMethod("parallelSlotNames", "GRanges",
    #function(x) c("seqnames", "ranges", "strand", callNextMethod())

    ## TEMPORARY DEFINITION.
    ## TODO: Remove this temporary definition and add "parallelSlotNames"
    ## methods to all packages that define "extraColumnSlotNames" methods
    ## (e.g. VariantAnnotation, GenomicTuples, InteractionSet, SGSeq).
    function(x) c(extraColumnSlotNames(x), "seqnames", "ranges", "strand",
                  callNextMethod())
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
        ### Fix old GRanges instances on-the-fly.
        object <- updateObject(object)
        BiocGenerics:::replaceSlots(object, ...)
    }
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
        version <- .get_GRanges_version(object)
        if (version == "current") {
            if (verbose)
                message("[updateObject] Internal representation of ",
                        class(object), " object is current.\n",
                        "[updateObject] Nothing to update.")
            return(object)
        }
        if (verbose)
            message("[updateObject] ", class(object), " object uses ",
                    "internal representation from GenomicRanges\n",
                    "[updateObject] ", version, ". Updating it ...")
        object@elementType <- new("GRanges")@elementType
        object
    }
)


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

### Genomic range squeezer.
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
    ans_strand <- as.character(tails(split0, n=-2L))
    ans_strand[is.na(ans_strand)] <- "*"
    split1 <- heads(split0, n=2L)
    ans_seqnames <- as.character(heads(split1, n=1L))
    ranges <- tails(split1, n=-1L)
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
        ans_mcols <- replaceROWS(mcols(x), i, mcols(value))
        ans_ecs_names <- extraColumnSlotNames(x)
        ans_necs <- length(ans_ecs_names)
        if (ans_necs == 0L) {
            ans_ecs <- NULL
        } else {
            value_ecs_names <- extraColumnSlotNames(value)
            if (!identical(head(value_ecs_names, n=ans_necs),
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

