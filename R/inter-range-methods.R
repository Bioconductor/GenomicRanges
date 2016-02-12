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

### 'rgl2' and 'rgl' must be List objects (typically RangesList or GRangesList)
### of the same length, both with a "revmap" inner metadata column.
.fix_revmap_inner_mcol <- function(rgl2, rgl)
{
    unlisted_rgl2 <- unlist(rgl2, use.names=FALSE)
    unlisted_revmap2 <- mcols(unlisted_rgl2)$revmap
    revmap <- relist(mcols(unlist(rgl, use.names=FALSE))$revmap, rgl)
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
    ans_seqnames <- Rle(structure(spaceid %/% 3L + 1L, levels=seqlevels(x),
                                  class="factor"), rgl_eltNROWS)
    ans_strand <- Rle(structure(spaceid %% 3L + 1L, levels=levels(strand()),
                                class="factor"), rgl_eltNROWS)

    ## Prepare 'ans_mcols'.
    ans_mcols <- mcols(ans_ranges)
    if (is.null(ans_mcols)) {
        ans_mcols <- new("DataFrame", nrows=length(ans_ranges))
    } else {
        mcols(ans_ranges) <- NULL
    }

    ## Prepare 'ans_seqinfo'.
    ans_seqinfo <- seqinfo(x)

    ## To be as fast as possible, we don't use internal constructor
    ## newGRanges() and we don't check the new object.
    new2("GRanges", seqnames=ans_seqnames,
                    ranges=ans_ranges,
                    strand=ans_strand,
                    elementMetadata=ans_mcols,
                    seqinfo=ans_seqinfo,
                    check=FALSE)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### range()
###

### Always return a GRanges instance (whatever GenomicRanges derivative 'x'
### is), so is NOT an endomorphism. 
setMethod("range", "GenomicRanges",
    function(x, ..., ignore.strand=FALSE, na.rm=FALSE)
    {
        if (!identical(na.rm, FALSE))
            warning("'na.rm' argument is ignored")
        args <- unname(list(x, ...))
        args <- lapply(args, granges)
        x <- do.call(c, args)

        rgl <- deconstructGRintoRGL(x, ignore.strand=ignore.strand, drop=TRUE)
        rgl2 <- callGeneric(rgl)
        reconstructGRfromRGL(rgl2, x)
    }
)

setMethod("range", "GRangesList",
    function(x, ..., ignore.strand=FALSE, na.rm=FALSE)
    {
        gr <- deconstructGRLintoGR(x)
        ## "range" method for GRanges objects is fast.
        gr2 <- callGeneric(gr, ..., ignore.strand=ignore.strand, na.rm=na.rm)
        reconstructGRLfromGR(gr2, x)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### reduce()
###

### Always return a GRanges instance (whatever GenomicRanges derivative 'x'
### is), so is NOT an endomorphism. 
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
            rgl2 <- .fix_revmap_inner_mcol(rgl2, rgl)
        reconstructGRfromRGL(rgl2, x)
    }
)

### TODO: Support the 'with.revmap' argument.
setMethod("reduce", "GRangesList",
    function(x, drop.empty.ranges=FALSE, min.gapwidth=1L,
             with.inframe.attrib=FALSE, ignore.strand=FALSE)
    {
        if (!identical(with.inframe.attrib, FALSE)) 
            stop("'with.inframe.attrib' argument is not supported ", 
                 "when reducing a GRangesList object")
        gr <- deconstructGRLintoGR(x)
        ## "reduce" method for GRanges objects is fast.
        gr2 <- callGeneric(gr, drop.empty.ranges=drop.empty.ranges,
                               min.gapwidth=min.gapwidth,
                               ignore.strand=ignore.strand)
        reconstructGRLfromGR(gr2, x)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### gaps()
###

### Always return a GRanges instance (whatever GenomicRanges derivative 'x'
### is), so is NOT an endomorphism. 
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

### Always return a GRanges instance (whatever GenomicRanges derivative 'x'
### is), so is NOT an endomorphism. 
setMethod("disjoin", "GenomicRanges",
    function(x, ignore.strand=FALSE)
    {
        rgl <- deconstructGRintoRGL(x, ignore.strand=ignore.strand, drop=TRUE)
        rgl2 <- callGeneric(rgl)
        reconstructGRfromRGL(rgl2, x)
    }
)

setMethod("disjoin", "GRangesList",
    function(x, ignore.strand=FALSE)
    {
        gr <- deconstructGRLintoGR(x)
        gr2 <- callGeneric(gr, ignore.strand=ignore.strand)
        reconstructGRLfromGR(gr2, x)
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

setMethod("isDisjoint", "GPos",
    function(x, ignore.strand=FALSE) callGeneric(x@pos_runs, ignore.strand)
)

setMethod("isDisjoint", "GRangesList",
    function(x, ignore.strand=FALSE)
    {
        gr <- deconstructGRLintoGR(x, expand.levels=TRUE)
        rgl <- deconstructGRintoRGL(gr, ignore.strand=ignore.strand)
        ans <- callGeneric(rgl)
        ans <- colSums(matrix(!ans, ncol=length(x))) == 0L
        names(ans) <- names(x)
        ans
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

