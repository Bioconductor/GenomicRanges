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
        return("slot 'ranges' cannot have element metadata")
    NULL
}

.valid.GRanges.elementMetadata <- function(x)
{
    if (!is.null(rownames(x@elementMetadata)))
        return("slot 'elementMetadata' cannot contain row names")
    NULL
}

.valid.GRanges <- function(x)
{
    c(.valid.GRanges.ranges(x), .valid.GRanges.elementMetadata(x))
}

setValidity2("GRanges", .valid.GRanges)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor.
###

### TODO: Revisit this constructor to make it more user friendly.
### Also provide a way to supply the sequence circularity flags.
newGRanges <- ## hidden constructor shared with other GRanges-like objects
function(class, seqnames = Rle(), ranges = IRanges(),
         strand = Rle("*", length(seqnames)),
         ...,
         seqlengths =
         structure(rep(NA_integer_, length(levels(seqnames))),
                   names = levels(seqnames)))
{
    ## occurs first for generation of default seqlengths
    if (!is(seqnames, "Rle"))
        seqnames <- Rle(seqnames)
    if (!is.factor(runValue(seqnames))) 
        runValue(seqnames) <- factor(runValue(seqnames))

    ranges <- force(ranges)
    if (class(ranges) != "IRanges")
        ranges <- as(ranges, "IRanges")

    if (!is(strand, "Rle"))
        strand <- Rle(strand)
    if (!is.factor(runValue(strand)) ||
        !identical(levels(runValue(strand)), levels(strand())))
        runValue(strand) <- strand(runValue(strand))
    if (IRanges:::anyMissing(runValue(strand))) {
        warning("missing values in strand converted to \"*\"")
        runValue(strand)[is.na(runValue(strand))] <- "*"
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

    if (!is.integer(seqlengths))
        seqlengths <-
          structure(as.integer(seqlengths), names = names(seqlengths))
    seqinfo <- Seqinfo(seqnames=names(seqlengths), seqlengths=seqlengths)
    ## in case we have seqlengths for unrepresented sequences
    runValue(seqnames) <- factor(runValue(seqnames), names(seqlengths))

    if (!is.null(elementMetadata(ranges))) {
        warning("'ranges' has element metadata, dropping them")
        elementMetadata(ranges) <- NULL
    }
    elementMetadata <- DataFrame(...)
    if (ncol(elementMetadata) == 0L)
        elementMetadata <- new("DataFrame", nrows = length(seqnames))
    if (!is.null(rownames(elementMetadata))) {
        if (!is.null(names(ranges)))
            names(ranges) <- rownames(elementMetadata)
        rownames(elementMetadata) <- NULL
    }

    new(class, seqnames = seqnames, ranges = ranges, strand = strand,
        seqinfo = seqinfo, elementMetadata = elementMetadata)
}

GRanges <-
  function(seqnames = Rle(), ranges = IRanges(),
           strand = Rle("*", length(seqnames)),
           ...,
           seqlengths = # this default is more accurate than in newGRanges
           structure(rep(NA_integer_, length(unique(seqnames))),
                     names = levels(as.factor(runValue(as(seqnames, "Rle"))))))
{
  if (missing(seqlengths)) # avoid potentially expensive seqnames conversion
    newGRanges("GRanges", seqnames = seqnames, ranges = ranges, strand = strand,
               ...)
  else newGRanges("GRanges", seqnames = seqnames, ranges = ranges,
                  strand = strand, ..., seqlengths = seqlengths)
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

setAs("RangedData", "GRanges",
    function(from)
    {
        ranges <- unlist(ranges(from), use.names=FALSE)
        values <- unlist(values(from), use.names=FALSE)
        nms <- rownames(from)
        rownames(values) <- NULL
        whichStrand <- which(colnames(values) == "strand")
        if (length(whichStrand) > 0)
            values <- values[-whichStrand]
        gr <- GRanges(seqnames = space(from),
                      ranges = ranges,
                      strand = Rle(strand(from)),
                      values)
        seqinfo(gr) <- seqinfo(from)
        metadata(gr) <- metadata(from)
        gr
    }
)

### Does NOT propagate the ranges names and elementMetadata i.e. always
### returns an unnamed GRanges object with no elementMetadata.
setAs("RangesList", "GRanges",
      function(from)
      {
        if (!length(from))
          return(GRanges())
        ranges <- unlist(from, use.names=FALSE)
        ranges <- IRanges(start=start(ranges), width=width(ranges))
        ## From now, ranges is guaranteed to be an IRanges *instance*.
        gr <- GRanges(seqnames = space(from),
                      ranges = ranges,
                      strand = Rle("*", length(ranges)))
        seqinfo(gr) <- seqinfo(from)
        metadata(gr) <- metadata(from)
        gr
      })

setAs("RleList", "GRanges", function(from) {
  rd <- as(from, "RangedData")
  rd$strand <- "*"
  gr <- as(rd, "GRanges")
  seqlengths(gr) <- elementLengths(from)
  gr
})

setAs("RleViewsList", "GRanges", function(from) {
  as(as(from, "RangedData"), "GRanges")
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Slot getters and setters.
###

setMethod("seqnames", "GRanges", function(x) x@seqnames)
setMethod("ranges", "GRanges", function(x, ...) x@ranges)
setMethod("strand", "GRanges", function(x) x@strand)
setMethod("seqinfo", "GRanges", function(x) x@seqinfo)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Combining and Splitting
###

setMethod("split", "GRanges",
    function(x, f, drop = FALSE, ...)
    {
        if (missing(f)) {
            nms <- names(x)
            if (is.null(nms))
                f <- seq_len(length(x))
            else
                f <- factor(nms, levels = nms)
        }
        nrows <- nlevels(f)
        if (nrows == 0)
            nrows <- sum(!is.na(unique(f)))
        IRanges:::newCompressedList("GRangesList", x, splitFactor = f,
                                    drop = drop,
                                    elementMetadata =
                                    new("DataFrame", nrows = nrows))
    }
)

