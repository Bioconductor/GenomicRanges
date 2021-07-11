### =========================================================================
### The GenomicRanges interface
### -------------------------------------------------------------------------
###

### TODO: The 'constraint' slot could be moved to the Vector class (or to the
### Annotated class) so any Vector object could be constrained.
setClass("GenomicRanges",
    contains="Ranges",
    representation(
        "VIRTUAL",
        #No more constraint slot for now...
        #constraint="Constraint_OR_NULL",
        elementMetadata="DataFrame"
    ),
    prototype(
        elementMetadata=new("DFrame")
    )
)

setClassUnion("GenomicRanges_OR_missing", c("GenomicRanges", "missing"))

### The code in this file will work out-of-the-box on 'x' as long as
### seqnames(x), ranges(x), strand(x), seqlengths(x), seqinfo(),
### and update(x) are defined.

setClass("GenomicPos",
    contains=c("GenomicRanges", "Pos"),
    representation("VIRTUAL")
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Getters.
###

### The method for Ranges does 'length(start(x))' which is very inefficient
### on a StitchedGPos object so we overwrite it with a method that should be
### fast on all GenomicRanges derivatives (including StitchedGPos objects).
setMethod("length", "GenomicRanges", function(x) length(ranges(x)))

setMethod("names", "GenomicRanges", function(x) names(ranges(x)))

#setMethod("constraint", "GenomicRanges", function(x) x@constraint)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Extra column slots (added by GRanges subclasses)
###
### The "extra column slots" are parallel slots added by GRanges subclasses
### in addition to GRanges parallel slots "seqnames", "ranges", "strand",
### and "elementMetadata".
###
### TODO: Packages defining "extraColumnSlotNames" methods (e.g.
### VariantAnnotation, GenomicTuples, InteractionSet, SGSeq) should define
### parallel_slot_names() methods instead. They'll get more things working
### out-of-the-box by doing so (e.g. c()). Then the extraColumnSlotNames()
### generic will become a GenomicRanges private business and could be made
### a regular function.
###

setGeneric("extraColumnSlotNames",
    function(x) standardGeneric("extraColumnSlotNames")
)

setMethod("extraColumnSlotNames", "ANY", function(x) character())

#setMethod("extraColumnSlotNames", "GRanges",
#    function(x)
#    {
#        GRanges_pslotnames <- parallel_slot_names(new("GRanges"))
#        setdiff(parallel_slot_names(x), GRanges_pslotnames)
#    }
#)

extraColumnSlots <- function(x)
    sapply(extraColumnSlotNames(x), slot, object=x, simplify=FALSE)

extraColumnSlotsAsDF <- function(x)
    S4Vectors:::new_DataFrame(extraColumnSlots(x), nrows=length(x))


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
             "Note that ranges located on a sequence whose length is unknown ",
             "(NA) or on a circular sequence are not considered out-of-bound ",
             "(use seqlengths() and isCircular() to get the lengths and ",
             "circularity flags of the underlying sequences).")
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
    if (!(class(ranges(x)) %in% c("IRanges", "StitchedIPos", "UnstitchedIPos")))
        return(paste0("'ranges(x)' must be an IRanges, StitchedIPos, ",
                      "or UnstitchedIPos instance"))
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
### TODO: Get rid of this restriction. Not sure why we ever had it. Doesn't
### make much sense to me.
INVALID.GR.COLNAMES <- c("seqnames", "ranges", "strand",
                         "seqlevels", "seqlengths", "isCircular",
                         #"genome",
                         "start", "end", "width", "element")

.valid.GenomicRanges.mcols <- function(x)
{
    if (any(INVALID.GR.COLNAMES %in% colnames(mcols(x, use.names=FALSE)))) {
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

.valid.GenomicRanges <- function(x)
{
    c(.valid.GenomicRanges.seqnames(x),
      .valid.GenomicRanges.ranges(x),
      .valid.GenomicRanges.strand(x),
      .valid.GenomicRanges.mcols(x),
      valid.GenomicRanges.seqinfo(x, suggest.trim=TRUE))
      #checkConstraint(x, constraint(x)))
}

setValidity2("GenomicRanges", .valid.GenomicRanges)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercion
###

### Propagate the names.
setMethod("as.character", "GenomicRanges",
    function(x, ignore.strand=FALSE)
    {
        if (!isTRUEorFALSE(ignore.strand))
            stop(wmsg("'ignore.strand' must be TRUE or FALSE"))
        if (length(x) == 0L)
            return(setNames(character(0), names(x)))
        ans <- paste0(seqnames(x), ":", as.character(ranges(x)))
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

setAs("GenomicRanges", "Grouping", function(from) {
    to <- as(grouping(seqnames(from), strand(from), start(from), end(from)),
             "ManyToOneGrouping")
    mcols(to)$granges <- granges(from)[end(PartitioningByEnd(to))]
    to
})

setMethod("as.data.frame", "GenomicRanges",
    function(x, row.names=NULL, optional=FALSE, ...)
    {
        if (missing(row.names))
            row.names <- names(x)
        if (!is.null(names(x)))
            names(x) <- NULL
        mcols_df <- as.data.frame(mcols(x, use.names=FALSE), ...)
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

.from_GenomicRanges_to_CompressedIRangesList <- function(from)
{
    strand_mcols <- DataFrame(strand=strand(from),
                              mcols(from, use.names=FALSE))
    ecs <- extraColumnSlotsAsDF(from)
    if (length(ecs))
      strand_mcols <- cbind(strand_mcols, ecs)
    rngs <- ranges(from)
    mcols(rngs) <- strand_mcols
    ans <- split(rngs, seqnames(from))
    mcols(ans) <- DataFrame(seqlengths=seqlengths(from),
                           isCircular=isCircular(from),
                           genome=genome(from))
    metadata(ans) <- metadata(from)
    metadata(ans)$seqinfo <- seqinfo(from)
    ans
}

setAs("GenomicRanges", "CompressedIRangesList",
    .from_GenomicRanges_to_CompressedIRangesList
)
setAs("GenomicRanges", "IRangesList",
    .from_GenomicRanges_to_CompressedIRangesList
)
setAs("GenomicRanges", "IntegerRangesList",
    .from_GenomicRanges_to_CompressedIRangesList
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Setters.
###

setReplaceMethod("names", "GenomicRanges",
    function(x, value)
    {
        x_ranges <- setNames(ranges(x), value)
        update(x, ranges=x_ranges, check=FALSE)
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
        update(x, seqnames=value)
    }
)

setReplaceMethod("ranges", "GenomicRanges",
    function(x, value)
    {
        if (class(value) != "IRanges")
            value <- as(value, "IRanges")
        mcols(value) <- NULL
        value <- S4Vectors:::V_recycle(value, x,
                                       x_what="value", skeleton_what="x")
        update(x, ranges=value)
    }
)

setReplaceMethod("strand", "GenomicRanges",
    function(x, value)
    {
        value <- normalize_strand_replacement_value(value, x)
        x <- update(x, strand=value, check=FALSE)
        msg <- .valid.GenomicRanges.strand(x)
        if (!is.null(msg))
            stop(msg)
        x
    }
)

### Does NOT suppoprt pruning mode "fine". Pruning modes "coarse" and "tidy"
### are equivalent on a GenomicRanges object.
set_GenomicRanges_seqinfo <-
    function(x, new2old=NULL,
             pruning.mode=c("error", "coarse", "fine", "tidy"),
             value)
{
    if (!is(value, "Seqinfo"))
        stop("the supplied 'seqinfo' must be a Seqinfo object")
    pruning.mode <- match.arg(pruning.mode)
    if (pruning.mode == "fine")
        stop(wmsg("\"fine\" pruning mode is not supported on ",
                  class(x), " objects"))
    dangling_seqlevels <- GenomeInfoDb:::getDanglingSeqlevels(x,
                              new2old=new2old,
                              pruning.mode=pruning.mode,
                              seqlevels(value))
    if (length(dangling_seqlevels) != 0L) {
        ## Prune 'x'.
        idx <- !(seqnames(x) %in% dangling_seqlevels)
        x <- x[idx]
    }
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
setReplaceMethod("seqinfo", "GenomicRanges", set_GenomicRanges_seqinfo)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The score() accessor
###
### Kind of a silly accessor. And why is it defined at the GenomicRanges
### level and not at the Vector level? Or at least at the Ranges and
### RangesList levels so it works on IRanges and IRangesList objects too.
###

setMethod("score", "GenomicRanges",
    function(x) mcols(x, use.names=FALSE)$score
)

setReplaceMethod("score", "GenomicRanges",
    function(x, value)
    {
        ## Fix old GRanges instances on-the-fly.
        x <- updateObject(x)
        mcols(x)$score <- value
        x
    }
)

#setReplaceMethod("constraint", "GenomicRanges",
#    function(x, value)
#    {
#        if (isSingleString(value))
#            value <- new(value)
#        if (!is(value, "Constraint_OR_NULL"))
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


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Range accessors
###

setMethod("start", "GenomicRanges", function(x, ...) start(ranges(x)))
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

setMethod("[", c("list_OR_List", "GenomicRanges"),
    function(x, i, j, ..., drop=TRUE)
    {
        if (!missing(j) || length(list(...)) > 0L)
            stop("invalid subsetting")
        .subset_by_GenomicRanges(x, i)
    }
)

### Avoid infinite recursion that we would otherwise get:
###   GRanges("chr1:1-10")[[1]]
###   # Error: C stack usage  7969428 is too close to the limit
setMethod("getListElement", "GenomicRanges",
    function(x, i, exact=TRUE)
    {
        ## A temporary situation
        stop(wmsg(class(x), " objects don't support [[, as.list(), ",
                  "lapply(), or unlist() at the moment"))
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

.DollarNames.GenomicRanges <- function(x, pattern = "")
    grep(pattern, names(mcols(x, use.names=FALSE)), value=TRUE)

setMethod("$", "GenomicRanges",
    function(x, name) mcols(x, use.names=FALSE)[[name]]
)

setReplaceMethod("$", "GenomicRanges",
    function(x, name, value)
    {
        ## Fix old GRanges instances on-the-fly.
        x <- updateObject(x)
        mcols(x)[[name]] <- value
        x
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Display
###

### S3/S4 combo for summary.GenomicRanges
summary.GenomicRanges <- summary.IPosRanges
setMethod("summary", "GenomicRanges", summary.GenomicRanges)

.from_GenomicRanges_to_naked_character_matrix_for_display <- function(x)
{
    m <- cbind(seqnames=showAsCell(seqnames(x)),
               ranges=showAsCell(ranges(x)),
               strand=showAsCell(strand(x)))
    extra_col_names <- extraColumnSlotNames(x)
    if (length(extra_col_names) != 0L) {
        extra_cols <- lapply(extraColumnSlots(x), showAsCell)
        m <- do.call(cbind, c(list(m), extra_cols))
    }
    cbind_mcols_for_display(m, x)
}
setMethod("makeNakedCharacterMatrixForDisplay", "GenomicRanges",
    .from_GenomicRanges_to_naked_character_matrix_for_display
)

### If 'x' is a GRanges object, 'coerce.internally.to.GRanges' has no effect.
### If it's a GenomicRanges object that is not a GRanges object, then
### show_GenomicRanges() will coerce it to a GRanges object unless
### 'coerce.internally.to.GRanges' is set to FALSE. Use this if coercing 'x'
### to GRanges is not supported or is too expensive but only if 'x' supports
### head() and tail().
show_GenomicRanges <- function(x, margin="",
                               print.classinfo=FALSE, print.seqinfo=FALSE,
                               coerce.internally.to.GRanges=TRUE)
{
    cat(margin, summary(x), ":\n", sep="")
    ## makePrettyMatrixForCompactPrinting() assumes that head() and tail()
    ## work on 'xx'.
    if (coerce.internally.to.GRanges) {
        xx <- as(x, "GRanges", strict=FALSE)
    } else {
        xx <- x
    }
    out <- makePrettyMatrixForCompactPrinting(xx)
    if (print.classinfo) {
        .COL2CLASS <- c(
            seqnames="Rle",
            ranges="IRanges",
            strand="Rle"
        )
        extra_col_names <- extraColumnSlotNames(x)
        .COL2CLASS <- c(.COL2CLASS, getSlots(class(x))[extra_col_names])
        classinfo <- makeClassinfoRowForCompactPrinting(x, .COL2CLASS)
        ## A sanity check, but this should never happen!
        stopifnot(identical(colnames(classinfo), colnames(out)))
        out <- rbind(classinfo, out)
    }
    if (nrow(out) != 0L)
        rownames(out) <- paste0(margin, "  ", rownames(out))
    ## We set 'max' to 'length(out)' to avoid the getOption("max.print")
    ## limit that would typically be reached when 'showHeadLines' global
    ## option is set to Inf.
    print(out, quote=FALSE, right=TRUE, max=length(out))
    if (print.seqinfo) {
        cat(margin, "  -------\n", sep="")
        cat(margin, "  seqinfo: ", summary(seqinfo(x)), "\n", sep="")
    }
}

setMethod("show", "GenomicRanges",
    function(object)
        show_GenomicRanges(object, print.classinfo=TRUE, print.seqinfo=TRUE)
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Concatenation
###

### 'x' can be any list-like object and its list elements can be anything
### that supports seqinfo().
### Used in the seqinfo() method for GenomicRangesList objects.
combine_seqinfo_from_GenomicRanges_objects <- function(x)
{
    x_len <- length(x)
    ## This will typically happen in the context of the seqinfo() method
    ## for GenomicRangesList objects.
    if (x_len == 0L)
        stop(wmsg("seqinfo() is not supported on a zero length ",
                  class(x), " object"))
    do.call(merge, lapply(seq_len(x_len), function(i) seqinfo(x[[i]])))
}

### NOT exported but used in the GenomicAlignments package.
concatenate_GenomicRanges_objects <-
    function(x, objects=list(), use.names=TRUE, ignore.mcols=FALSE, check=TRUE)
{
    ## Fix old GRanges and GAlignments instance on-the-fly.
    x <- updateObject(x, check=FALSE)

    objects <- S4Vectors:::prepare_objects_to_bind(x, objects)
    all_objects <- c(list(x), objects)

    ## Combine seqinfo.
    seqinfo(x) <- combine_seqinfo_from_GenomicRanges_objects(all_objects)

    ## Call method for Vector objects to concatenate all the parallel slots.
    callNextMethod()
}

setMethod("bindROWS", "GenomicRanges", concatenate_GenomicRanges_objects)

