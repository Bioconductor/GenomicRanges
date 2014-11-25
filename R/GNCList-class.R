### =========================================================================
### GNCList objects
### -------------------------------------------------------------------------
###
### GNCList is a container for storing a preprocessed GenomicRanges object
### that can be used for fast findOverlaps().
###

setClass("GNCList",
    contains="GenomicRanges",
    representation(
        nclists="list",
        granges="GRanges"
    )
)

.get_circle_length <- function(x)
{
    circle_length <- seqlengths(x)
    circle_length[!(isCircular(x) %in% TRUE)] <- NA_integer_
    circle_length
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Accessors
###

setMethod("granges", "GNCList",
    function(x, use.mcols=FALSE)
    {
        if (!isTRUEorFALSE(use.mcols))
            stop("'use.mcols' must be TRUE or FALSE")
        ans <- x@granges
        if (use.mcols)
            mcols(ans) <- mcols(x)
        ans
    }
)

setMethod("length", "GNCList", function(x) length(granges(x)))
setMethod("names", "GNCList", function(x) names(granges(x)))
setMethod("seqnames", "GNCList", function(x) seqnames(granges(x)))
setMethod("start", "GNCList", function(x, ...) start(granges(x)))
setMethod("end", "GNCList", function(x, ...) end(granges(x)))
setMethod("width", "GNCList", function(x) width(granges(x)))
setMethod("ranges", "GNCList", function(x, use.mcols=FALSE) ranges(granges(x)))
setMethod("strand", "GNCList", function(x) strand(granges(x)))
setMethod("seqinfo", "GNCList", function(x) seqinfo(granges(x)))

setAs("GNCList", "GRanges", function(from) granges(from))
setAs("GNCList", "NCLists",
    function(from)
    {
        from_granges <- granges(from)
        ans_rglist <- split(ranges(from_granges), seqnames(from_granges))
        new2("NCLists", nclists=from@nclists, rglist=ans_rglist, check=FALSE)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor
###

GNCList <- function(x)
{
    if (!is(x, "GenomicRanges"))
        stop("'x' must be a GenomicRanges object")
    if (!is(x, "GRanges"))
        x <- as(x, "GRanges")
    mcols(x) <- NULL
    x_groups <- split(seq_len(length(x)) - 1L, seqnames(x))
    x_nclists <- IRanges:::NCList_by_group(ranges(x),
                                           x_groups,
                                           .get_circle_length(x))
    new2("GNCList", nclists=x_nclists, granges=x, check=FALSE)
}

setAs("GenomicRanges", "GNCList", function(from) GNCList(from))


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### findOverlaps_GNCList()
###

### Rearrange the seqlevels in 'x' so that its first N seqlevels are exactly
### 'seqlevels' (in the same order).
.align_seqlevels <- function(x, seqlevels)
{
    seqlevels(x) <- union(seqlevels, seqlevels(x))
    x
}

### NOT exported.
findOverlaps_GNCList <- function(query, subject, min.score=1L,
                                 type=c("any", "start", "end",
                                        "within", "extend", "equal"),
                                 select=c("all", "first", "last", "arbitrary"),
                                 ignore.strand=FALSE)
{
    if (!(is(query, "GenomicRanges") && is(subject, "GenomicRanges")))
        stop("'query' and 'subject' must be GenomicRanges objects")
    if (!isSingleNumber(min.score))
        stop("'min.score' must be a single integer")
    if (!is.integer(min.score))
        min.score <- as.integer(min.score)
    type <- match.arg(type)
    select <- match.arg(select)
    if (!isTRUEorFALSE(ignore.strand))
        stop("'ignore.strand' must be TRUE or FALSE")

    q_len <- length(query)
    s_len <- length(subject)
    q_seqlevels <- seqlevels(query)
    s_seqlevels <- seqlevels(subject)
    if (is(subject, "GNCList")) {
        si <- merge(seqinfo(subject), seqinfo(query))
        query <- .align_seqlevels(query, s_seqlevels)
        nclists <- subject@nclists
        nclist_is_q <- rep.int(FALSE, length(nclists))
    } else if (is(query, "GNCList")) {
        si <- merge(seqinfo(query), seqinfo(subject))
        subject <- .align_seqlevels(subject, q_seqlevels)
        nclists <- query@nclists
        nclist_is_q <- rep.int(TRUE, length(nclists))
    } else {
        ## We'll do "on-the-fly preprocessing".
        if (length(s_seqlevels) <= length(q_seqlevels)) {
            si <- merge(seqinfo(subject), seqinfo(query))
            query <- .align_seqlevels(query, s_seqlevels)
            nclists <- vector(mode="list", length=length(s_seqlevels))
        } else {
            si <- merge(seqinfo(query), seqinfo(subject))
            subject <- .align_seqlevels(subject, q_seqlevels)
            nclists <- vector(mode="list", length=length(q_seqlevels))
        }
        nclist_is_q <- rep.int(NA, length(nclists))
    }
    if (ignore.strand) {
        q_space <- s_space <- NULL
    } else {
        q_space <- as.integer(strand(query)) - 3L
        s_space <- as.integer(strand(subject)) - 3L
    }
    q_groups <- split(seq_len(q_len) - 1L, seqnames(query))
    s_groups <- split(seq_len(s_len) - 1L, seqnames(subject))
    circle_length <- .get_circle_length(si)
    IRanges:::NCList_find_overlaps_by_group_and_combine(
                          ranges(query), q_space, q_groups,
                          ranges(subject), s_space, s_groups,
                          nclists, nclist_is_q,
                          min.score, type, select, circle_length)
}

