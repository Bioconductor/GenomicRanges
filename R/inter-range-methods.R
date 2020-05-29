### =========================================================================
### Inter-range methods
### -------------------------------------------------------------------------
###

### 'revmap_unlisted' and 'revmap_partitioning': IntegerList and Partitioning
### objects representing the "revmap object", which is *conceptually* an
### IntegerListList object (of the same length as 'revmap_partitioning').
### *Conceptually* because, well, we don't actually have such container...
### 'old2new': IntegerList of the same length as the "revmap object".
.translate_revmap <- function(revmap_unlisted, revmap_partitioning, old2new)
{
    ## 'times' has the length of the "revmap object".
    times <- sum(relist(width(PartitioningByEnd(revmap_unlisted)),
                        revmap_partitioning))
    ## 'offset' has the length of 'revmap_unlisted@unlistData'.
    offset <- rep.int(start(PartitioningByEnd(old2new)) - 1L, times)
    revmap_flat <- revmap_unlisted@unlistData
    revmap_unlisted@unlistData <- old2new@unlistData[revmap_flat + offset]
    revmap_unlisted
}

### 'rgl2' and 'rgl' must be List objects (typically IntegerRangesList or
### GRangesList) of the same length, both with a "revmap" inner metadata
### column.
.fix_inner_revmap_mcol <- function(rgl2, rgl)
{
    unlisted_rgl2 <- unlist(rgl2, use.names=FALSE)
    unlisted_revmap2 <- mcols(unlisted_rgl2, use.names=FALSE)$revmap
    revmap <- relist(mcols(unlist(rgl, use.names=FALSE), use.names=FALSE)$revmap, rgl)
    mcols(unlisted_rgl2)$revmap <- .translate_revmap(unlisted_revmap2, rgl2,
                                                     revmap)
    rgl2 <- relist(unlisted_rgl2, rgl2)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Deconstruction/reconstruction of a GenomicRanges object into/from an
### IRangesList object.
###
### For internal use only (not exported).
###

.make_spaceid <- function(x, ignore.strand=FALSE, drop=FALSE)
{
    spaceid <- seqnames(x)
    runValue(spaceid) <- 3L * as.integer(runValue(spaceid))
    if (!ignore.strand) {
        strandid <- strand(x)
        runValue(strandid) <- as.integer(runValue(strandid)) - 3L
        spaceid <- spaceid + strandid
    }
    if (!drop) {
        levels <- as.character(seq_len(3L * length(seqlevels(x))))
        runValue(spaceid) <- structure(runValue(spaceid),
                                       levels=levels,
                                       class="factor")
    }
    spaceid
}

### Work on any GenomicRanges derivative and return a CompressedIRangesList
### instance.
deconstructGRintoRGL <- function(x, with.revmap=FALSE,
                                    ignore.strand=FALSE, drop=FALSE)
{
    if (!isTRUEorFALSE(with.revmap))
        stop("'with.revmap' must be TRUE or FALSE")
    if (!isTRUEorFALSE(ignore.strand))
        stop("'ignore.strand' must be TRUE or FALSE")
    x_ranges <- unname(ranges(x))
    if (with.revmap)
        mcols(x_ranges) <- DataFrame(revmap=seq_along(x_ranges))
    x_spaceid <- .make_spaceid(x, ignore.strand=ignore.strand, drop=drop)
    split(x_ranges, x_spaceid)
}

### Return a GRanges instance.
reconstructGRfromRGL <- function(rgl, x)
{
    ## Prepare 'ans_ranges'.
    ans_ranges <- unlist(rgl, use.names=FALSE)

    ## Prepare 'ans_seqnames' and 'ans_strand'.
    rgl_eltNROWS <- elementNROWS(rgl)
    spaceid <- as.integer(names(rgl)) - 1L
    ans_seqnames <- Rle(structure(spaceid %/% 3L + 1L,
                                  levels=seqlevels(x),
                                  class="factor"),
                        rgl_eltNROWS)
    ans_strand <- Rle(structure(spaceid %% 3L + 1L,
                                levels=levels(strand()),
                                class="factor"),
                      rgl_eltNROWS)

    ## Prepare 'ans_mcols'.
    ans_mcols <- mcols(ans_ranges, use.names=FALSE)
    if (is.null(ans_mcols)) {
        ans_mcols <- make_zero_col_DFrame(length(ans_ranges))
    } else {
        mcols(ans_ranges) <- NULL
    }

    ## Prepare 'ans_seqinfo'.
    ans_seqinfo <- seqinfo(x)

    ## To be as fast as possible, we don't use internal low-level constructor
    ## new_GRanges() and we don't check the new object.
    new2("GRanges", seqnames=ans_seqnames,
                    ranges=ans_ranges,
                    strand=ans_strand,
                    elementMetadata=ans_mcols,
                    seqinfo=ans_seqinfo,
                    check=FALSE)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Deconstruction/reconstruction of a GRangesList object into/from a GRanges
### object.
###
### For internal use only (not exported).
###

### 'f' is assumed to be an integer vector of non-negative values with no NAs.
.pad_with_zeros <- function(f)
{
    stopifnot(is.integer(f))
    if (length(f) == 0L)
        return(character())
    f_max <- max(f)
    if (f_max <= 9L)
        return(as.character(f))
    nd <- as.integer(log10(f_max)) + 1L
    sprintf(paste0("%0", nd, "d"), f)
}

.paste12 <- function(f1, f2)
{
    paste(.pad_with_zeros(f1), .pad_with_zeros(f2), sep="|")
}

### Unlist GRangesList object 'x' into a GRanges object but the differences
### with the "unlist" method for GRangesList objects are:
###   - The sequence names of the returned GRanges object are modified by
###     embedding the "grouping by top-level element" information in them.
###   - The seqinfo is modified accordingly.
deconstructGRLintoGR <- function(x, expand.levels=FALSE)
{
    ans <- x@unlistData
    f1 <- rep.int(seq_along(x), elementNROWS(x))
    f2 <- as.integer(seqnames(ans))
    f12 <- .paste12(f1, f2)

    ## Compute 'ans_seqinfo'.
    if (expand.levels) {
        x_nlev <- length(seqlevels(x))
        i1 <- rep(seq_len(length(x)), each=x_nlev)
        i2 <- rep.int(seq_len(x_nlev), length(x))
    } else {
        oo <- orderIntegerPairs(f1, f2)
        of1 <- f1[oo]
        of2 <- f2[oo]
        ## TODO: Support 'method="presorted"' in duplicatedIntegerPairs() for
        ## when the 2 input vectors are already sorted.
        notdups <- !duplicatedIntegerPairs(of1, of2)
        i1 <- of1[notdups]
        i2 <- of2[notdups]
    }
    x_seqinfo <- seqinfo(x)
    ans_seqlevels <- .paste12(i1, i2)
    ans_seqlengths <- unname(seqlengths(x_seqinfo))[i2]
    ans_isCircular <- unname(isCircular(x_seqinfo))[i2]
    ans_seqinfo <- Seqinfo(ans_seqlevels, ans_seqlengths, ans_isCircular)

    ## The 2 following modifications must be seen as a single atomic
    ## operation since doing the 1st without doing the 2nd would leave 'ans'
    ## in a broken state.
    ans@seqnames <- Rle(factor(f12, ans_seqlevels))
    ans@seqinfo <- ans_seqinfo
    ans
}

### The "inverse" transform of deconstructGRLintoGR().
### More precisely, reconstructGRLfromGR() transforms GRanges object 'gr'
### with sequence names in the "f1|f2" format (as produced by
### deconstructGRLintoGR() above) back into a GRangesList object with the
### same length & names & metadata columns & seqinfo as 'x'.
### The fundamental property of this deconstruction/reconstruction mechanism
### is that, for any GRangesList object 'x':
###
###   reconstructGRLfromGR(deconstructGRLintoGR(x), x) is identical to x
###
reconstructGRLfromGR <- function(gr, x, with.revmap=FALSE)
{
    snames <- strsplit(as.character(seqnames(gr)), "|", fixed=TRUE)
    m12 <- matrix(as.integer(unlist(snames)), ncol=2, byrow=TRUE)

    ## Restore the real sequence names.
    f2 <- m12[ , 2L]
    x_seqlevels <- seqlevels(x)
    ## The 2 following modifications must be seen as a single atomic
    ## operation since doing the 1st without doing the 2nd would leave 'ans'
    ## in a broken state.
    gr@seqnames <- Rle(factor(x_seqlevels[f2], x_seqlevels))
    gr@seqinfo <- seqinfo(x)

    ## Split.
    f1 <- m12[ , 1L]
    ans <- split(gr, factor(f1, levels=seq_len(length(x))))

    ## Decorate 'ans'.
    metadata(ans) <- metadata(x)
    names(ans) <- names(x)
    if (with.revmap) {
        unlisted_ans <- unlist(ans, use.names=FALSE)
        mcols(unlisted_ans)$revmap <-
            IRanges:::global2local_revmap(mcols(unlisted_ans, use.names=FALSE)$revmap, ans, x)
        ans <- relist(unlisted_ans, ans)
    }
    mcols(ans) <- mcols(x, use.names=FALSE)
    ans
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### range()
###

### Always return a GRanges *instance* whatever GenomicRanges derivative the
### input is, so does NOT act like an endomorphism in general. 
setMethod("range", "GenomicRanges",
    function(x, ..., with.revmap=FALSE, ignore.strand=FALSE, na.rm=FALSE)
    {
        if (!identical(na.rm, FALSE))
            warning("'na.rm' argument is ignored")
        args <- unname(list(x, ...))
        args <- lapply(args, granges)
        gr <- do.call(c, args)

        rgl <- deconstructGRintoRGL(gr, with.revmap=with.revmap,
                                    ignore.strand=ignore.strand, drop=TRUE)
        rgl2 <- callGeneric(rgl, with.revmap=with.revmap)
        if (with.revmap)
            rgl2 <- .fix_inner_revmap_mcol(rgl2, rgl)
        reconstructGRfromRGL(rgl2, gr)
    }
)

### Overwrite above method with optimized method for StitchedGPos objects.
### Like the above method, return a GRanges instance.
setMethod("range", "StitchedGPos",
    function(x, ..., with.revmap=FALSE, ignore.strand=FALSE, na.rm=FALSE)
        callGeneric(stitch_StitchedGPos(x), ...,
                    with.revmap=with.revmap, ignore.strand=ignore.strand,
                    na.rm=na.rm)
)

setMethod("range", "GRangesList",
    function(x, ..., with.revmap=FALSE, ignore.strand=FALSE, na.rm=FALSE)
    {
        if (class(x) == "GRangesList") {
            #warning(wmsg(OLD_GRANGESLIST_INSTANCE_MSG))
            x <- updateObject(x, check=FALSE)
        }
        if (is(x, "CompressedList")) {
            gr <- deconstructGRLintoGR(x)
            ## "range" method for GRanges objects is fast.
            gr2 <- callGeneric(gr, ..., with.revmap=with.revmap,
                               ignore.strand=ignore.strand, na.rm=na.rm)
            ans <- reconstructGRLfromGR(gr2, x, with.revmap=with.revmap)
            return(ans)
        }
        endoapply(x, range, ..., with.revmap=with.revmap,
                     ignore.strand=ignore.strand, na.rm=na.rm)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### reduce()
###

### Always return a GRanges *instance* whatever GenomicRanges derivative the
### input is, so does NOT act like an endomorphism in general. 
setMethod("reduce", "GenomicRanges",
    function(x, drop.empty.ranges=FALSE, min.gapwidth=1L,
                with.revmap=FALSE,
                with.inframe.attrib=FALSE, ignore.strand=FALSE)
    {
        if (!identical(with.inframe.attrib, FALSE))
            stop("'with.inframe.attrib' argument not supported ",
                 "when reducing a GenomicRanges object")
        rgl <- deconstructGRintoRGL(x, with.revmap=with.revmap,
                                       ignore.strand=ignore.strand, drop=TRUE)
        rgl2 <- callGeneric(rgl, drop.empty.ranges=drop.empty.ranges,
                                 min.gapwidth=min.gapwidth,
                                 with.revmap=with.revmap)
        if (with.revmap)
            rgl2 <- .fix_inner_revmap_mcol(rgl2, rgl)
        reconstructGRfromRGL(rgl2, x)
    }
)

setMethod("reduce", "GRangesList",
    function(x, drop.empty.ranges=FALSE, min.gapwidth=1L,
                with.revmap=FALSE,
                with.inframe.attrib=FALSE, ignore.strand=FALSE)
    {
        if (class(x) == "GRangesList") {
            #warning(wmsg(OLD_GRANGESLIST_INSTANCE_MSG))
            x <- updateObject(x, check=FALSE)
        }
        if (!identical(with.inframe.attrib, FALSE)) 
            stop("'with.inframe.attrib' argument is not supported ", 
                 "when reducing a GRangesList object")
        if (is(x, "CompressedList")) {
            gr <- deconstructGRLintoGR(x)
            gr2 <- callGeneric(gr, drop.empty.ranges=drop.empty.ranges,
                                   min.gapwidth=min.gapwidth,
                                   with.revmap=with.revmap,
                                   ignore.strand=ignore.strand)
            ans <- reconstructGRLfromGR(gr2, x, with.revmap=with.revmap)
            return(ans)
        }
        endoapply(x, reduce, drop.empty.ranges=drop.empty.ranges,
                             min.gapwidth=min.gapwidth,
                             with.revmap=with.revmap,
                             ignore.strand=ignore.strand)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### gaps()
###

### Always return a GRanges *instance* whatever GenomicRanges derivative the
### input is, so does NOT act like an endomorphism in general. 
setMethod("gaps", "GenomicRanges",
    function(x, start=1L, end=seqlengths(x))
    {
        seqlevels <- seqlevels(x)
        if (!is.null(names(start)))
            start <- start[seqlevels]
        if (!is.null(names(end)))
            end <- end[seqlevels]
        start <- S4Vectors:::recycleVector(start, length(seqlevels))
        start <- rep(start, each=3L)
        end <- S4Vectors:::recycleVector(end, length(seqlevels))
        end <- rep(end, each=3L)
        rgl <- deconstructGRintoRGL(x)
        rgl2 <- callGeneric(rgl, start=start, end=end)
        reconstructGRfromRGL(rgl2, x)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### disjoin()
###

### Always return a GRanges *instance* whatever GenomicRanges derivative the
### input is, so does NOT act like an endomorphism in general. 
setMethod("disjoin", "GenomicRanges",
    function(x, with.revmap=FALSE, ignore.strand=FALSE)
    {
        rgl <- deconstructGRintoRGL(x, with.revmap=with.revmap, 
                                       ignore.strand=ignore.strand, drop=TRUE)
        rgl2 <- callGeneric(rgl, with.revmap=with.revmap)
        if (with.revmap)
            rgl2 <- .fix_inner_revmap_mcol(rgl2, rgl)
        reconstructGRfromRGL(rgl2, x)
    }
)

setMethod("disjoin", "GRangesList",
    function(x, with.revmap=FALSE, ignore.strand=FALSE)
    {
        if (class(x) == "GRangesList") {
            #warning(wmsg(OLD_GRANGESLIST_INSTANCE_MSG))
            x <- updateObject(x, check=FALSE)
        }
        if (is(x, "CompressedList")) {
            gr <- deconstructGRLintoGR(x)
            gr2 <- callGeneric(gr, with.revmap=with.revmap,
                                   ignore.strand=ignore.strand)
            ans <- reconstructGRLfromGR(gr2, x, with.revmap=with.revmap)
            return(ans)
        }
        endoapply(x, disjoin, with.revmap=with.revmap,
                     ignore.strand=ignore.strand)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### isDisjoint()
###

setMethod("isDisjoint", "GenomicRanges",
    function(x, ignore.strand=FALSE)
    {
        rgl <- deconstructGRintoRGL(x, ignore.strand=ignore.strand, drop=TRUE)
        all(callGeneric(rgl))
    }
)

### Overwrite above method with optimized method for StitchedGPos objects.
setMethod("isDisjoint", "StitchedGPos",
    function(x, ignore.strand=FALSE)
        callGeneric(stitch_StitchedGPos(x), ignore.strand)
)

setMethod("isDisjoint", "GRangesList",
    function(x, ignore.strand=FALSE)
    {
        if (class(x) == "GRangesList") {
            #warning(wmsg(OLD_GRANGESLIST_INSTANCE_MSG))
            x <- updateObject(x, check=FALSE)
        }
        if (is(x, "CompressedList")) {
            gr <- deconstructGRLintoGR(x, expand.levels=TRUE)
            rgl <- deconstructGRintoRGL(gr, ignore.strand=ignore.strand)
            ans <- callGeneric(rgl)
            ans <- colSums(matrix(!ans, ncol=length(x))) == 0L
            names(ans) <- names(x)
            return(ans)
        }
        endoapply(x, isDisjoint, ignore.strand=ignore.strand)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### disjointBins()
###

setMethod("disjointBins", "GenomicRanges",
    function(x, ignore.strand=FALSE)
    {
        x_spaceid <- .make_spaceid(x, ignore.strand=ignore.strand, drop=TRUE)
        rgl <- split(unname(ranges(x)), x_spaceid)
        ans <- callGeneric(rgl)
        unsplit(ans, x_spaceid)
    }
)

