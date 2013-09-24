### =========================================================================
### The GenomicRanges interface
### -------------------------------------------------------------------------
###

### TODO: The 'constraint' slot could be moved to the Vector class (or to the
### Annotated class) so any Vector object could be constrained.
setClass("GenomicRanges",
    contains="Vector",
    representation(
        "VIRTUAL"#,
        #No more constraint slot for now...
        #constraint="ConstraintORNULL"
    )
)

setClassUnion("GenomicRangesORmissing", c("GenomicRanges", "missing"))

### The code in this file will work out-of-the-box on 'x' as long as
### seqnames(x), ranges(x), strand(x), seqlengths(x), seqinfo(),
### update(x) and clone(x) are defined.


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### 2 non-exported low-level helper functions.
###
### For the 2 functions below, 'x' must be a GenomicRanges object.
### They both return a named integer vector where the names are guaranteed
### to be 'seqlevels(x)'.
###

minStartPerGRangesSequence <- function(x)
{
    cil <- splitAsList(start(x), seqnames(x))  # CompressedIntegerList object
    v <- Views(cil@unlistData, cil@partitioning)  # XIntegerViews object
    ans <- viewMins(v)
    ans[width(v) == 0L] <- NA_integer_
    names(ans) <- names(v)
    ans
}

maxEndPerGRangesSequence <- function(x)
{
    cil <- splitAsList(end(x), seqnames(x))  # CompressedIntegerList object
    v <- Views(cil@unlistData, cil@partitioning)  # XIntegerViews object
    ans <- viewMaxs(v)
    ans[width(v) == 0L] <- NA_integer_
    names(ans) <- names(v)
    ans
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Getters.
###

setMethod("length", "GenomicRanges", function(x) length(seqnames(x)))

setMethod("names", "GenomicRanges", function(x) names(ranges(x)))

#setMethod("constraint", "GenomicRanges", function(x) x@constraint)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Extra column slots (implemented by subclasses)
###

extraColumnSlots <- function(x) {
  sapply(extraColumnSlotNames(x), slot, object = x, simplify = FALSE)
}

extraColumnSlotsAsDF <- function(x) {
  ## low-level fast path; otherwise, would need to wrap some things with I()
  new("DataFrame", listData = extraColumnSlots(x), nrows = length(x))
}

setGeneric("extraColumnSlotNames",
           function(x) standardGeneric("extraColumnSlotNames"))

setMethod("extraColumnSlotNames", "ANY", function(x) character())

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity.
###
### TODO: Should we enforce ranges(x) to be an IRanges *instance* (i.e.
### class(ranges(x) == "IRanges")) instead of just an IRanges *object* (i.e.
### is(ranges(x), "IRanges"))? Right now I can create a GRanges object where
### the ranges are a Views object on a very long DNAString subject with
### something like 'GRanges("chr1", Views(subject, start=1:2, end=5))'.
### Sounds cool but there are also some potential complications with this...

.valid.GenomicRanges.length <- function(x)
{
    x_len <- length(x)
    checkCoreGetterReturnedLength <- function(getter) {
        if (NROW(get(getter)(x)) != x_len)
            paste0("NROW(", getter, "(x)) != length(x)")
    }
    pbs1 <- unlist(lapply(c("seqnames", "ranges", "strand", "mcols"),
                          checkCoreGetterReturnedLength))
    checkExtraColumnLength <- function(slotname) {
        if (NROW(slot(x, slotname)) != x_len)
            paste0("NROW(x@", slotname, ") != length(x)")
    }
    pbs2 <- unlist(lapply(extraColumnSlotNames(x), checkExtraColumnLength))
    c(pbs1, pbs2)
}

.valid.GenomicRanges.seqnames <- function(x)
{
    if (!is.factor(runValue(seqnames(x))))
        return("'seqnames' should be a 'factor' Rle")
    if (IRanges:::anyMissing(runValue(seqnames(x))))
        return("'seqnames' contains missing values")
    NULL
}

.valid.GenomicRanges.ranges <- function(x)
{
    if (class(ranges(x)) != "IRanges")
        return("'ranges(x)' must be an IRanges instance")
    NULL
}

.valid.GenomicRanges.strand <- function(x)
{
    if (!is.factor(runValue(strand(x))) ||
        !identical(levels(runValue(strand(x))), levels(strand())))
    {
        msg <- c("'strand' should be a 'factor' Rle with levels c(",
                 paste0('"', levels(strand()), '"', collapse=", "),
                 ")")
        return(paste(msg, collapse=""))
    }
    if (IRanges:::anyMissing(runValue(strand(x))))
        return("'strand' contains missing values")
    NULL
}

### NOTE: This list is also included in the man page for GRanges objects.
### Keep the 2 lists in sync!
### We don't put "genome" in that list in order to facilitate import of GFF3
### files as GRanges objects (see ?import.gff3 in rtracklayer).
INVALID.GR.COLNAMES <- c("seqnames", "ranges", "strand",
                         "seqlevels", "seqlengths", "isCircular",
                         #"genome",
                         "start", "end", "width", "element")

.valid.GenomicRanges.mcols <- function(x)
{    
    if (any(INVALID.GR.COLNAMES %in% colnames(mcols(x)))) {
        msg <- c("names of metadata columns cannot be one of ",
                 paste0("\"", INVALID.GR.COLNAMES, "\"", collapse=", "))
        return(paste(msg, collapse=" "))
    }
    NULL
}

### Also used by the validity method for GAlignments objects.
valid.GenomicRanges.seqinfo <- function(x)
{
    x_seqinfo <- seqinfo(x)
    if (!identical(seqlevels(x_seqinfo), levels(seqnames(x)))) {
        msg <- c("'seqlevels(seqinfo(x))' and 'levels(seqnames(x))'",
                 "are not identical")
        return(paste(msg, collapse=" "))
    }
    x_seqlengths <- seqlengths(x_seqinfo)
    seqs_with_known_length <- names(x_seqlengths)[!is.na(x_seqlengths)]
    if (length(seqs_with_known_length) == 0L)
        return(NULL)
    if (any(x_seqlengths < 0L, na.rm=TRUE))
        return("'seqlengths(x)' contains negative values")
    ## We check only the ranges that are on a non-circular sequence with
    ## a known length.
    x_isCircular <- isCircular(x_seqinfo)
    non_circ_seqs <- names(x_isCircular)[!(x_isCircular %in% TRUE)]
    ncswkl <- intersect(non_circ_seqs, seqs_with_known_length)
    if (length(ncswkl) == 0L)
        return(NULL)
    x_seqnames <- seqnames(x)
    runValue(x_seqnames) <- runValue(x_seqnames)[drop=TRUE]
    minStarts <- minStartPerGRangesSequence(x)
    maxEnds <- maxEndPerGRangesSequence(x)
    if (any(minStarts[ncswkl] < 1L, na.rm=TRUE)
     || any(maxEnds[ncswkl] >
            seqlengths(x)[ncswkl], na.rm=TRUE))
        warning("'ranges' contains values outside of sequence bounds. ",
                "See ?trim to subset ranges.")
    NULL
}

## For convenience, validate the extra column slots that are virtual
## classes. Since they are not directly constructed, any validity
## checks specific to the virtual class have probably not been called.
.valid.GenomicRanges.ecs <- function(x) {
  virtuals <- Filter(isVirtualClass, getSlots(class(x))[extraColumnSlotNames(x)])
  unlist(lapply(names(virtuals), function(nm) validObject(slot(x, nm))))
}

.valid.GenomicRanges <- function(x)
{
    c(.valid.GenomicRanges.length(x),
      .valid.GenomicRanges.seqnames(x),
      .valid.GenomicRanges.ranges(x),
      .valid.GenomicRanges.strand(x),
      .valid.GenomicRanges.mcols(x),
      valid.GenomicRanges.seqinfo(x),
      .valid.GenomicRanges.ecs(x))
      #checkConstraint(x, constraint(x)))
}

setValidity2("GenomicRanges", .valid.GenomicRanges)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercion.
###

setAs("GenomicRanges", "RangedData",
    function(from)
    {
        mcols <- mcols(from)
        ecs <- extraColumnSlotsAsDF(from)
        if (length(ecs))
          mcols <- cbind(mcols, ecs)
        rd <- RangedData(ranges(from), strand=strand(from),
                         mcols, space=seqnames(from))
        mcols(ranges(rd)) <- DataFrame(seqlengths=seqlengths(from),
                                       isCircular=isCircular(from),
                                       genome=genome(from))
        metadata(ranges(rd)) <- metadata(from)
        metadata(ranges(rd))$seqinfo <- seqinfo(from)
        rd
    }
)

setAs("GenomicRanges", "RangesList",
    function(from)
    {
        strand_mcols <- DataFrame(strand=strand(from), mcols(from))
        ecs <- extraColumnSlotsAsDF(from)
        if (length(ecs))
          strand_mcols <- cbind(strand_mcols, ecs)
        rngs <- ranges(from)
        mcols(rngs) <- strand_mcols
        rl <- split(rngs, seqnames(from))
        mcols(rl) <- DataFrame(seqlengths=seqlengths(from),
                               isCircular=isCircular(from),
                               genome=genome(from))
        metadata(rl) <- metadata(from)
        metadata(rl)$seqinfo <- seqinfo(from)
        rl
    }
)

setMethod("as.data.frame", "GenomicRanges",
    function(x, row.names=NULL, optional=FALSE, ...)
    {
        ranges <- ranges(x)
        if (missing(row.names))
            row.names <- names(x)
        if (!is.null(names(x)))
            names(x) <- NULL
        mcols_df <- as.data.frame(mcols(x), ...)
        if (length(extraColumnSlotNames(x)) > 0L)
            mcols_df <- cbind(as.data.frame(extraColumnSlotsAsDF(x), ...),
                              mcols_df)
        data.frame(seqnames=as.factor(seqnames(x)),
                   start=start(x),
                   end=end(x),
                   width=width(x),
                   strand=as.factor(strand(x)),
                   mcols_df,
                   row.names=row.names,
                   stringsAsFactors=FALSE)
    }
)

.fromSeqinfoToGRanges <- function(from)
{
    if (any(is.na(seqlengths(from))))
        stop("cannot create a GenomicRanges from a Seqinfo ",
             "with NA seqlengths")
    GRanges(seqnames(from),
            IRanges(rep(1L, length(from)),
                    width=seqlengths(from),
                    names=seqnames(from)),
            seqinfo=from)
}

setAs("Seqinfo", "GenomicRanges", .fromSeqinfoToGRanges)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Setters.
###

setReplaceMethod("names", "GenomicRanges",
    function(x, value)
    {
        names(ranges(x)) <- value
        x
    }
)

setReplaceMethod("seqnames", "GenomicRanges",
    function(x, value) 
    {
        if (!is(value, "Rle"))
            value <- Rle(value)
        if (!is.factor(runValue(value)))
            runValue(value) <- factor(runValue(value))
        if (!identical(levels(value), seqlevels(x)))
            stop("levels of supplied 'seqnames' must be ",
                 "identical to 'seqlevels(x)'")
        n <- length(x)
        k <- length(value)
        if (k != n) {
            if ((k == 0L) || (k > n) || (n %% k != 0L))
                stop(k, " elements in value to replace ", n, " elements")
            value <- rep(value, length.out=n)
        }
        update(x, seqnames=value)
    }
)

setReplaceMethod("ranges", "GenomicRanges",
    function(x, value) 
    {
        if (!is(value, "IRanges"))
            value <- as(value, "IRanges")
        n <- length(x)
        k <- length(value)
        if (k != n) {
            if ((k == 0L) || (k > n) || (n %% k != 0L))
                stop(k, " elements in value to replace ", n, "elements")
            value <- rep(value, length.out=n)
        }
        update(x, ranges=value, check=FALSE)
    }
)

normargGenomicRangesStrand <- function(strand, n)
{
    if (!is(strand, "Rle"))
        strand <- Rle(strand)
    if (!is.factor(runValue(strand))
     || !identical(levels(runValue(strand)), levels(strand())))
        runValue(strand) <- strand(runValue(strand))
    k <- length(strand)
    if (k != n) {
        if (k != 1L && (k == 0L || k > n || n %% k != 0L))
            stop("supplied 'strand' has ", k, " elements (", n, " expected)")
        strand <- rep(strand, length.out=n)
    }
    strand
}

setReplaceMethod("strand", "GenomicRanges",
    function(x, value) 
    {
        value <- normargGenomicRangesStrand(value, length(x))
        x <- update(x, strand=value, check=FALSE)
        msg <- .valid.GenomicRanges.strand(x)
        if (!is.null(msg))
            stop(msg)
        x
    }
)

setReplaceMethod("elementMetadata", "GenomicRanges",
    function(x, ..., value)
    {
        value <- normalizeMetadataColumnsReplacementValue(value, x)
        x <- update(x, elementMetadata=value, check=FALSE)
        msg <- .valid.GenomicRanges.mcols(x)
        if (!is.null(msg))
            stop(msg)
        x
    }
)

setReplaceMethod("seqinfo", "GenomicRanges",
    function(x, new2old=NULL, force=FALSE, value)
    {
        if (!is(value, "Seqinfo"))
            stop("the supplied 'seqinfo' must be a Seqinfo object")
        dangling_seqlevels <- getDanglingSeqlevels(x,
                                  new2old=new2old, force=force,
                                  seqlevels(value))
        if (length(dangling_seqlevels) != 0L)
            x <- x[!(seqnames(x) %in% dangling_seqlevels)]
        old_seqinfo <- seqinfo(x)
        new_seqnames <- makeNewSeqnames(x, new2old=new2old, seqlevels(value))
        x <- update(x, seqnames=new_seqnames, seqinfo=value, check=FALSE)
        geom_has_changed <- sequenceGeometryHasChanged(seqinfo(x), old_seqinfo,
                                                       new2old=new2old)
        if (any(geom_has_changed, na.rm=TRUE)) {
            msg <- valid.GenomicRanges.seqinfo(x)
            if (!is.null(msg))
                stop(msg)
        }
        x
    }
)

setMethod("score", "GenomicRanges", function(x) mcols(x)$score)

#setReplaceMethod("constraint", "GenomicRanges",
#    function(x, value)
#    {
#        if (isSingleString(value))
#            value <- new(value)
#        if (!is(value, "ConstraintORNULL"))
#            stop("the supplied 'constraint' must be a ",
#                 "Constraint object, a single string, or NULL")
#        x@constraint <- value
#        validObject(x)
#        x
#    }
#)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Updating and cloning.
###
### An object is either 'update'd in place (usually with a replacement
### method) or 'clone'd (copied), with specified slots/fields overridden.

### For an object with a pure S4 slot representation, these both map to
### initialize. Reference classes will want to override 'update'. Other
### external representations need further customization.

setGeneric("clone", function(x, ...) standardGeneric("clone"))  # not exported

setMethod("clone", "ANY",  # not exported
    function(x, ...)
    {
        if (nargs() > 1L)
            update(x, ...)
        else
            x
    }
)

setMethod("update", "GenomicRanges",
          function(object, ...)
          {
            BiocGenerics:::updateS4(object, ...)
          })

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Ranges methods.
###

setMethod("start", "GenomicRanges", function(x, ...) start(ranges(x)))
setMethod("end", "GenomicRanges", function(x, ...) end(ranges(x)))
setMethod("width", "GenomicRanges", function(x) width(ranges(x)))

setReplaceMethod("start", "GenomicRanges",
    function(x, check=TRUE, value)
    {
        if (!is.integer(value))
            value <- as.integer(value)
        ranges <- ranges(x)
        starts <- start(ranges)
        starts[] <- value
        ## TODO: Revisit this to handle circularity (maybe).
        if (!IRanges:::anyMissing(seqlengths(x))) {
            if (IRanges:::anyMissingOrOutside(starts, 1L)) {
                warning("trimmed start values to be positive")
                starts[starts < 1L] <- 1L
            }
        }
        start(ranges, check=check) <- starts
        update(x, ranges=ranges, check=FALSE)
    }
)

setReplaceMethod("end", "GenomicRanges",
    function(x, check=TRUE, value)
    {
        if (!is.integer(value))
            value <- as.integer(value)
        ranges <- ranges(x)
        ends <- end(ranges)
        ends[] <- value
        seqlengths <- seqlengths(x)
        ## TODO: Revisit this to handle circularity.
        if (!IRanges:::anyMissing(seqlengths)) {
            seqlengths <- seqlengths[seqlevels(x)]
            maxEnds <- seqlengths[as.integer(seqnames(x))]
            trim <- which(ends > maxEnds)
            if (length(trim) > 0L) {
                warning("trimmed end values to be <= seqlengths")
                ends[trim] <- maxEnds[trim]
            }
        }
        end(ranges, check=check) <- ends
        update(x, ranges=ranges, check=FALSE)
    }
)

setReplaceMethod("width", "GenomicRanges",
    function(x, check=TRUE, value)
    {
        if (!is.integer(value))
            value <- as.integer(value)
        if (!IRanges:::anyMissing(seqlengths(x))) {
            end(x) <- start(x) + (value - 1L)
        } else {
            ranges <- ranges(x)
            width(ranges, check=check) <- value
            x <- update(x, ranges=ranges, check=FALSE)
        }
        x
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Subsetting and combining.
###
        
setMethod(IRanges:::extractROWS, "GenomicRanges",
    function(x, i)
    {
        if (missing(i) || !is(i, "Ranges"))
            i <- IRanges:::normalizeSingleBracketSubscript(i, x)
        ans_seqnames <- IRanges:::extractROWS(seqnames(x), i)
        ans_ranges <- IRanges:::extractROWS(ranges(x), i)
        ans_strand <- IRanges:::extractROWS(strand(x), i)
        ans_mcols <- IRanges:::extractROWS(mcols(x), i)
        ans_ecs <- lapply(extraColumnSlots(x), IRanges:::extractROWS, i)
        clone(x, seqnames=ans_seqnames,
                 ranges=ans_ranges,
                 strand=ans_strand,
                 elementMetadata=ans_mcols,
                 .slotList=ans_ecs)
    }
)

### Needed only because we want to support x[i, j] subsetting.
setMethod("[", "GenomicRanges",
    function(x, i, j, ..., drop)
    {
        if (length(list(...)) > 0L)
            stop("invalid subsetting")
        ans <- IRanges:::extractROWS(x, i)
        if (missing(j))
            return(ans)
        ans_mcols <- mcols(ans)[ , j, drop=FALSE]
        clone(ans, elementMetadata=ans_mcols)
    }
)

setMethod(IRanges:::replaceROWS, "GenomicRanges",
    function(x, i, value)
    {
        if (missing(i) || !is(i, "Ranges"))
            i <- IRanges:::normalizeSingleBracketSubscript(i, x)
        seqinfo(x) <- merge(seqinfo(x), seqinfo(value))
        ans_seqnames <- IRanges:::replaceROWS(seqnames(x), i, seqnames(value))
        ans_ranges <- IRanges:::replaceROWS(ranges(x), i, ranges(value))
        ans_strand <- IRanges:::replaceROWS(strand(x), i, strand(value))
        ans_mcols <- IRanges:::replaceROWS(mcols(x), i, mcols(value))
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
            ans_ecs <- IRanges:::replaceROWS(ans_ecs, i,
                                             value_ecs[seq_len(ans_necs)])
        }
        update(x, seqnames=ans_seqnames,
                  ranges=ans_ranges,
                  strand=ans_strand,
                  elementMetadata=ans_mcols,
                  .slotList=as.list(ans_ecs))
    }
)

setReplaceMethod("[", "GenomicRanges",
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
            i <- IRanges:::extractROWS(setNames(seq_along(x), names(x)), i)
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
        update(x, seqnames=seqnames, ranges=ranges,
               strand=strand, elementMetadata=ans_mcols,
               .slotList=as.list(x_ecs))
    }
)

### Not exported. 'x' *must* be an unnamed list of length >= 1 (not checked).
.unlist_list_of_GenomicRanges <- function(x, ignore.mcols=FALSE)
{
    if (!isTRUEorFALSE(ignore.mcols))
        stop("'ignore.mcols' must be TRUE or FALSE")
    ans_class <- class(x[[1L]])
    ans_seqinfo <- do.call(merge, lapply(x, seqinfo))
    ans_seqnames <- do.call(c, lapply(x, seqnames))
    ans_ranges <- do.call(c, lapply(x, ranges))
    ans_strand <- do.call(c, lapply(x, strand))
    if (ignore.mcols) {
        ans_mcols <- new("DataFrame", nrows=length(ans_ranges))
    } else {
        ans_mcols <- do.call(rbind, lapply(x, mcols, FALSE))
    }
    if (length(extraColumnSlotNames(x[[1L]])) > 0L) {
        ans_ecs <- do.call(rbind, lapply(x, extraColumnSlotsAsDF))
        do.call(new, c(list(ans_class, seqnames=ans_seqnames, ranges=ans_ranges,
                            strand=ans_strand, elementMetadata=ans_mcols,
                            seqinfo=ans_seqinfo), as.list(ans_ecs)))
    } else {
        new(ans_class, seqnames=ans_seqnames, ranges=ans_ranges,
            strand=ans_strand, elementMetadata=ans_mcols, seqinfo=ans_seqinfo)
    }
}

setMethod("c", "GenomicRanges",
    function(x, ..., ignore.mcols=FALSE, recursive=FALSE)
    {
        if (!identical(recursive, FALSE))
            stop("'recursive' argument not supported")
        if (missing(x))
            args <- unname(list(...))
        else
            args <- unname(list(x, ...))
        .unlist_list_of_GenomicRanges(args, ignore.mcols=ignore.mcols)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### $ and $<- methods
###
### Provided as a convenience, for GenomicRanges *only*, and as the result
### of strong popular demand.
### Note that those methods are not consistent with the other $ and $<-
### methods in the IRanges/GenomicRanges infrastructure, and might confuse
### some users by making them believe that a GenomicRanges object can be
### manipulated as a data.frame-like object.
### Therefore we recommend using them only interactively, and we discourage
### their use in scripts or packages. For the latter, use 'mcols(x)$name'
### instead of 'x$name'.
###

.DollarNames.GenomicRanges <- function(x, pattern)
    grep(pattern, names(mcols(x)), value=TRUE)

setMethod("$", "GenomicRanges",
    function(x, name) mcols(x)[[name]]
)

setReplaceMethod("$", "GenomicRanges",
    function(x, name, value) {mcols(x)[[name]] <- value; x}
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### "show" method.
###

.makeNakedMatFromGenomicRanges <- function(x)
{
    lx <- length(x)
    nc <- ncol(mcols(x))
    ans <- cbind(seqnames=as.character(seqnames(x)),
                 ranges=IRanges:::showAsCell(ranges(x)),
                 strand=as.character(strand(x)))
    extraColumnNames <- extraColumnSlotNames(x)
    if (length(extraColumnNames) > 0L) {
        ans <- do.call(cbind,
                       c(list(ans), lapply(extraColumnSlots(x), showAsCell)))
    }
    if (nc > 0L) {
        tmp <- do.call(data.frame, c(lapply(mcols(x),
                                            IRanges:::showAsCell),
                                     list(check.names=FALSE)))
        ans <- cbind(ans, `|`=rep.int("|", lx), as.matrix(tmp))
    }
    ans
}

showGenomicRanges <- function(x, margin="",
                              with.classinfo=FALSE, print.seqlengths=FALSE)
{
    lx <- length(x)
    nc <- ncol(mcols(x))
    cat(class(x), " with ",
        lx, " ", ifelse(lx == 1L, "range", "ranges"),
        " and ",
        nc, " metadata ", ifelse(nc == 1L, "column", "columns"),
        ":\n", sep="")
    out <- IRanges:::makePrettyMatrixForCompactPrinting(x,
               .makeNakedMatFromGenomicRanges)
    if (with.classinfo) {
        .COL2CLASS <- c(
            seqnames="Rle",
            ranges="IRanges",
            strand="Rle"
        )
        extraColumnNames <- extraColumnSlotNames(x)
        .COL2CLASS <- c(.COL2CLASS, getSlots(class(x))[extraColumnNames])
        classinfo <- makeClassinfoRowForCompactPrinting(x, .COL2CLASS)
        ## A sanity check, but this should never happen!
        stopifnot(identical(colnames(classinfo), colnames(out)))
        out <- rbind(classinfo, out)
    }
    if (nrow(out) != 0L)
        rownames(out) <- paste0(margin, rownames(out))
    print(out, quote=FALSE, right=TRUE)
    if (print.seqlengths) {
        cat(margin, "---\n", sep="")
        showSeqlengths(x, margin=margin)
    }
}

setMethod("show", "GenomicRanges",
    function(object)
        showGenomicRanges(object, margin="  ",
                          with.classinfo=TRUE, print.seqlengths=TRUE)
)

setMethod("showAsCell", "GenomicRanges",
          function(object) {
            if (length(object) > 0L) {
              paste0(seqnames(object), ":", strand(object), ":",
                     showAsCell(ranges(object)))
            } else character(0)
          })
