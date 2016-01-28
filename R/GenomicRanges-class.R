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

setMethod("parallelVectorNames", "GenomicRanges", function(x) {
  colnames(as.data.frame(new(class(x))))
})


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### 2 low-level helper functions to deal with out-of-bound ranges.
###

### Returns index of out-of-bound ranges located on non-circular sequences
### whose length is not NA. Works on a GenomicRanges or GAlignments object.
get_out_of_bound_index <- function(x)
{
    if (length(x) == 0L)
        return(integer(0))
    x_seqnames_id <- as.integer(seqnames(x))
    x_seqlengths <- unname(seqlengths(x))
    seqlevel_is_circ <- unname(isCircular(x)) %in% TRUE
    seqlength_is_na <- is.na(x_seqlengths)
    seqlevel_has_bounds <- !(seqlevel_is_circ | seqlength_is_na)
    which(seqlevel_has_bounds[x_seqnames_id] &
          (start(x) < 1L | end(x) > x_seqlengths[x_seqnames_id]))
}

### Also works on a GenomicRanges or GAlignments object. Note that GAlignments
### objects are not trimmable so use 'suggest.trim=FALSE' on them.
make_out_of_bound_warning_msg <- function(x, idx, suggest.trim)
{
    where <- seqlevels(x)[unique(as.integer(seqnames(x))[idx])]
    if (length(where) == 1L) {
        on_what <- paste0("sequence ", where)
    } else if (length(where) == 2L) {
        on_what <- paste0("sequences ", where[1L], " and ", where[2L])
    } else {
        seqlevels_in1string <- paste0(head(where, n=-1L), collapse=", ")
        on_what <- paste0("sequences ", seqlevels_in1string,
                              ", and ", tail(where, n=1L))
    }
    msg <- c(class(x), " object contains ", length(idx), " out-of-bound ",
             "range", if (length(idx) >= 2L) "s" else "", " located on ",
             on_what, ". ",
             "Note that only ranges located on a non-circular ",
             "sequence whose length is not NA can be considered ",
             "out-of-bound (use seqlengths() and isCircular() to ",
             "get the lengths and circularity flags of the underlying ",
             "sequences).")
    if (suggest.trim)
        msg <- c(msg, " You can use trim() to trim these ranges. ",
                 "See ?`trim,GenomicRanges-method` for more information.")
    msg
}


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

### Used in GenomicAlignments.
.valid.GenomicRanges.seqnames <- function(x)
{
    x_seqnames <- seqnames(x)
    if (!is(x_seqnames, "Rle") || !is.factor(runValue(x_seqnames)))
        return("'seqnames(x)' must be a 'factor' Rle")
    if (!is.null(names(x_seqnames)))
        return("'seqnames(x)' must be a 'factor' Rle with no names")
    if (S4Vectors:::anyMissing(runValue(x_seqnames)))
        return("'seqnames(x)' contains missing values")
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
    if (S4Vectors:::anyMissing(runValue(strand(x))))
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
valid.GenomicRanges.seqinfo <- function(x, suggest.trim=FALSE)
{
    x_seqinfo <- seqinfo(x)
    if (!identical(seqlevels(x_seqinfo), levels(seqnames(x)))) {
        msg <- c("'seqlevels(seqinfo(x))' and 'levels(seqnames(x))'",
                 "are not identical")
        return(paste(msg, collapse=" "))
    }
    idx <- get_out_of_bound_index(x)
    if (length(idx) != 0L) {
        msg <- make_out_of_bound_warning_msg(x, idx, suggest.trim)
        warning(wmsg(msg))
    }
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
      valid.GenomicRanges.seqinfo(x, suggest.trim=TRUE),
      .valid.GenomicRanges.ecs(x))
      #checkConstraint(x, constraint(x)))
}

setValidity2("GenomicRanges", .valid.GenomicRanges)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercion.
###

### Propagate the names.
setMethod("as.character", "GenomicRanges",
    function(x, ignore.strand=FALSE)
    {
        if (!isTRUEorFALSE(ignore.strand))
            stop(wmsg("'ignore.strand' must be TRUE or FALSE"))
        if (length(x) == 0L)
            return(setNames(character(0), names(x)))
        ans <- paste0(seqnames(x), ":", start(x), "-", end(x))
        names(ans) <- names(x)
        if (ignore.strand)
            return(ans)
        x_strand <- strand(x)
        if (all(x_strand == "*"))
            return(ans)
        setNames(paste0(ans, ":", x_strand), names(x))
    }
)

### The as.factor() generic doesn't have the ... argument so this method
### cannot support the 'ignore.strand' argument.
setMethod("as.factor", "GenomicRanges",
    function(x)
        factor(as.character(x), levels=as.character(sort(unique(x))))
)

### TODO: Turn this into an S3/S4 combo for as.data.frame.GenomicRanges
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
setAs("Seqinfo", "RangesList",
    function(from) as(as(from, "GenomicRanges"), "RangesList")
)

setMethod("granges", "GenomicRanges",
    function(x, use.mcols=FALSE)
    {
        if (!isTRUEorFALSE(use.mcols))
            stop("'use.mcols' must be TRUE or FALSE")
        ans <- GRanges(seqnames(x), ranges(x), strand(x), seqinfo=seqinfo(x))
        if (use.mcols)
            mcols(ans) <- cbind(extraColumnSlotsAsDF(x), mcols(x))
        ans
    }
)


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

### Used in GenomicAlignments.
.normalize_seqnames_replacement_value <- function(value, x)
{
    if (!is(value, "Rle"))
        value <- Rle(value)
    if (!is.factor(runValue(value)))
        runValue(value) <- factor(runValue(value))
    if (!identical(levels(value), seqlevels(x)))
        stop("levels of supplied 'seqnames' must be ",
             "identical to 'seqlevels(x)'")
    S4Vectors:::V_recycle(value, x, x_what="value", skeleton_what="x")
}

setReplaceMethod("seqnames", "GenomicRanges",
    function(x, value)
    {
        value <- .normalize_seqnames_replacement_value(value, x)
        BiocGenerics:::replaceSlots(x, seqnames=value)
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
                stop(k, " elements in value to replace ", n, " elements")
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

setReplaceMethod("seqinfo", "GenomicRanges",
    function(x, new2old=NULL, force=FALSE, value)
    {
        if (!is(value, "Seqinfo"))
            stop("the supplied 'seqinfo' must be a Seqinfo object")
        dangling_seqlevels <- GenomeInfoDb:::getDanglingSeqlevels(x,
                                  new2old=new2old, force=force,
                                  seqlevels(value))
        if (length(dangling_seqlevels) != 0L)
            x <- x[!(seqnames(x) %in% dangling_seqlevels)]
        old_seqinfo <- seqinfo(x)
        new_seqnames <- GenomeInfoDb:::makeNewSeqnames(x,
                                  new2old=new2old, seqlevels(value))
        x <- update(x, seqnames=new_seqnames, seqinfo=value, check=FALSE)
        geom_has_changed <- GenomeInfoDb:::sequenceGeometryHasChanged(
                                  seqinfo(x), old_seqinfo, new2old=new2old)
        if (any(geom_has_changed, na.rm=TRUE)) {
            msg <- valid.GenomicRanges.seqinfo(x, suggest.trim=TRUE)
            if (!is.null(msg))
                stop(msg)
        }
        x
    }
)

setMethod("score", "GenomicRanges", function(x) mcols(x)$score)
setReplaceMethod("score", "GenomicRanges", function(x, value) {
  mcols(x)$score <- value
  x
})

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
            BiocGenerics:::replaceSlots(object, ...)
          })

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Ranges methods.
###

setMethod("start", "GenomicRanges", function(x, ...) start(ranges(x)))
setMethod("end", "GenomicRanges", function(x, ...) end(ranges(x)))
setMethod("width", "GenomicRanges", function(x) width(ranges(x)))

setReplaceMethod("start", "GenomicRanges",
    function(x, ..., value)
    {
        new_ranges <- `start<-`(ranges(x), ..., value=value)
        update(x, ranges=new_ranges, ...)
    }
)

setReplaceMethod("end", "GenomicRanges",
    function(x, ..., value)
    {
        new_ranges <- `end<-`(ranges(x), ..., value=value)
        update(x, ranges=new_ranges, ...)
    }
)

setReplaceMethod("width", "GenomicRanges",
    function(x, ..., value)
    {
        new_ranges <- `width<-`(ranges(x), ..., value=value)
        update(x, ranges=new_ranges, ...)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Subsetting.
###

setMethod("extractROWS", "GenomicRanges",
    function(x, i)
    {
        if (missing(i) || !is(i, "Ranges"))
            i <- normalizeSingleBracketSubscript(i, x)
        ans_seqnames <- extractROWS(seqnames(x), i)
        ans_ranges <- extractROWS(ranges(x), i)
        ans_strand <- extractROWS(strand(x), i)
        ans_mcols <- extractROWS(mcols(x), i)
        ans_ecs <- lapply(extraColumnSlots(x), extractROWS, i)
        clone(x, seqnames=ans_seqnames,
                 ranges=ans_ranges,
                 strand=ans_strand,
                 elementMetadata=ans_mcols,
                 .slotList=ans_ecs,
                 check=FALSE)
    }
)

### Needed only because we want to support x[i, j] subsetting.
setMethod("[", "GenomicRanges",
    function(x, i, j, ..., drop)
    {
        if (length(list(...)) > 0L)
            stop("invalid subsetting")
        if (!missing(i))
            x <- extractROWS(x, i)
        if (missing(j))
            return(x)
        new_mcols <- mcols(x)[ , j, drop=FALSE]
        clone(x, elementMetadata=new_mcols, check=FALSE)
    }
)

setMethod("replaceROWS", "GenomicRanges",
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
        update(x, seqnames=ans_seqnames,
                  ranges=ans_ranges,
                  strand=ans_strand,
                  elementMetadata=ans_mcols,
                  .slotList=as.list(ans_ecs))
    }
)

### TODO: Refactor to use replaceROWS(). This will make the code much simpler
### and avoid a lot of duplication with the above "replaceROWS" method.
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
        update(x, seqnames=seqnames, ranges=ranges,
               strand=strand, elementMetadata=ans_mcols,
               .slotList=as.list(x_ecs))
    }
)

### Subset a named list-like object *by* a GenomicRanges subscript.
### The returned object 'ans' is as follow:
###   (a) 'ans' is parallel to 'gr'.
###   (b) 'names(ans)' is identical to 'as.character(seqnames(gr))'.
###   (c) 'elementNROWS(ans)' is the same as 'width(gr)'.
###   (d) 'class(ans)' is 'relistToClass(x[[1]])' e.g. CompressedRleList if
###       'x' is an RleList object, or DNAStringSet is 'x' is a DNAStringSet
###       object.
.subset_by_GenomicRanges <- function(x, gr)
{
    if (!(is.list(x) || is(x, "List")))
        stop(wmsg("'x' must be a list-like object when subsetting ",
                  "by a GenomicRanges subscript"))
    x_names <- names(x)
    if (is.null(x_names))
        stop(wmsg("'x' must have names when subsetting ",
                  "by a GenomicRanges subscript"))
    if (anyDuplicated(x_names))
        stop(wmsg("'x' must have unique names when subsetting ",
                  "by a GenomicRanges subscript"))
    irl <- split(ranges(gr), seqnames(gr), drop=TRUE)
    seqlevels_in_use <- names(irl)
    seqlevels2names <- match(seqlevels_in_use, x_names)
    if (any(is.na(seqlevels2names)))
        stop(wmsg("when subsetting by a GenomicRanges subscript, the names ",
                  "of the object to subset must contain the seqnames of the ",
                  "subscript"))

    ## Handle empty case.
    if (length(gr) == 0L) {
        if (length(x) != 0L) {
            x1 <- x[[1L]]
        } else if (is(x, "CompressedList")) {
            x1 <- unlist(x, use.names=FALSE)
        } else {
            x1 <- new(elementType(x))
        }
        unlisted_ans <- x1[0]
        ans_partitioning <- PartitioningByEnd(names=character(0))
        return(relist(unlisted_ans, ans_partitioning))
    }

    tmp <- lapply(seq_along(seqlevels_in_use),
                  function(i) {
                      seqlevel <- seqlevels_in_use[i]
                      name <- x_names[seqlevels2names[i]]
                      extractList(x[[name]], irl[[seqlevel]])
                  })

    ## Unsplit 'tmp'.
    ans <- do.call(c, tmp)
    ans_len <- length(gr)
    idx <- unlist(split(seq_len(ans_len), seqnames(gr), drop=TRUE))
    revidx <- integer(ans_len)
    revidx[idx] <- seq_len(ans_len)
    names(ans) <- names(idx)
    ans <- ans[revidx]
    ans
}

setMethod("[", c("List", "GenomicRanges"),
    function(x, i, j, ..., drop=TRUE)
    {
        if (!missing(j) || length(list(...)) > 0L)
            stop("invalid subsetting")
        .subset_by_GenomicRanges(x, i)
    }
)

setMethod("[", c("list", "GenomicRanges"),
    function(x, i, j, ..., drop=TRUE)
    {
        if (!missing(j) || length(list(...)) > 0L)
            stop("invalid subsetting")
        .subset_by_GenomicRanges(x, i)
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
### Combining.
###

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
### "show" method.
###

.makeNakedMatFromGenomicRanges <- function(x)
{
    x_len <- length(x)
    x_mcols <- mcols(x)
    x_nmc <- if (is.null(x_mcols)) 0L else ncol(x_mcols)
    ans <- cbind(seqnames=as.character(seqnames(x)),
                 ranges=showAsCell(ranges(x)),
                 strand=as.character(strand(x)))
    extraColumnNames <- extraColumnSlotNames(x)
    if (length(extraColumnNames) > 0L) {
        ans <- do.call(cbind,
                       c(list(ans), lapply(extraColumnSlots(x), showAsCell)))
    }
    if (x_nmc > 0L) {
        tmp <- do.call(data.frame, c(lapply(x_mcols, showAsCell),
                                     list(check.names=FALSE)))
        ans <- cbind(ans, `|`=rep.int("|", x_len), as.matrix(tmp))
    }
    ans
}

showGenomicRanges <- function(x, margin="",
                              print.classinfo=FALSE, print.seqinfo=FALSE)
{
    x_class <- class(x)
    x_len <- length(x)
    x_mcols <- mcols(x)
    x_nmc <- if (is.null(x_mcols)) 0L else ncol(x_mcols)
    cat(x_class, " object with ",
        x_len, " ", ifelse(x_len == 1L, "range", "ranges"),
        " and ",
        x_nmc, " metadata ", ifelse(x_nmc == 1L, "column", "columns"),
        ":\n", sep="")
    ## S4Vectors:::makePrettyMatrixForCompactPrinting() assumes that 'x' is
    ## subsettable but not all GenomicRanges objects are (and if they are,
    ## subsetting them could be costly). However GRanges objects are assumed
    ## to be subsettable so if 'x' is not one then we turn it into one (this
    ## coercion is expected to work on any GenomicRanges object).
    if (!is(x, "GRanges"))
        x <- as(x, "GRanges")
    out <- S4Vectors:::makePrettyMatrixForCompactPrinting(x,
               .makeNakedMatFromGenomicRanges)
    if (print.classinfo) {
        .COL2CLASS <- c(
            seqnames="Rle",
            ranges="IRanges",
            strand="Rle"
        )
        extraColumnNames <- extraColumnSlotNames(x)
        .COL2CLASS <- c(.COL2CLASS, getSlots(class(x))[extraColumnNames])
        classinfo <-
            S4Vectors:::makeClassinfoRowForCompactPrinting(x, .COL2CLASS)
        ## A sanity check, but this should never happen!
        stopifnot(identical(colnames(classinfo), colnames(out)))
        out <- rbind(classinfo, out)
    }
    if (nrow(out) != 0L)
        rownames(out) <- paste0(margin, rownames(out))
    ## We set 'max' to 'length(out)' to avoid the getOption("max.print")
    ## limit that would typically be reached when 'showHeadLines' global
    ## option is set to Inf.
    print(out, quote=FALSE, right=TRUE, max=length(out))
    if (print.seqinfo) {
        cat(margin, "-------\n", sep="")
        cat(margin, "seqinfo: ", summary(seqinfo(x)), "\n", sep="")
    }
}

setMethod("show", "GenomicRanges",
    function(object)
        showGenomicRanges(object, margin="  ",
                          print.classinfo=TRUE, print.seqinfo=TRUE)
)

setMethod("showAsCell", "GenomicRanges",
          function(object) {
            if (length(object) > 0L) {
              paste0(seqnames(object), ":", strand(object), ":",
                     showAsCell(ranges(object)))
            } else character(0)
          })
