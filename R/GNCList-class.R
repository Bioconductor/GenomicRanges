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
    function(x, use.names=TRUE, use.mcols=FALSE)
    {
        if (!isTRUEorFALSE(use.names))
            stop("'use.names' must be TRUE or FALSE")
        if (!isTRUEorFALSE(use.mcols))
            stop("'use.mcols' must be TRUE or FALSE")
        ans <- x@granges
        if (!use.names)
            names(ans) <- NULL
        if (use.mcols)
            mcols(ans) <- mcols(x, use.names=FALSE)
        ans
    }
)

setMethod("length", "GNCList", function(x) length(granges(x)))
setMethod("names", "GNCList", function(x) names(granges(x)))
setMethod("seqnames", "GNCList", function(x) seqnames(granges(x)))
setMethod("start", "GNCList", function(x, ...) start(granges(x)))
setMethod("end", "GNCList", function(x, ...) end(granges(x)))
setMethod("width", "GNCList", function(x) width(granges(x)))
setMethod("ranges", "GNCList",
    function(x, use.names=TRUE, use.mcols=FALSE)
        ranges(granges(x, use.names=use.names, use.mcols=use.mcols),
               use.names=TRUE, use.mcols=use.mcols)
)
setMethod("strand", "GNCList", function(x) strand(granges(x)))
setMethod("seqinfo", "GNCList", function(x) seqinfo(granges(x)))

setAs("GNCList", "GRanges", function(from) granges(from, use.mcols=TRUE))


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor
###

.extract_groups_from_GenomicRanges <- function(x)
    splitAsList(seq_along(x) - 1L, seqnames(x))

GNCList <- function(x)
{
    if (!is(x, "GenomicRanges"))
        stop("'x' must be a GenomicRanges object")
    if (!is(x, "GRanges"))
        x <- as(x, "GRanges")
    ans_mcols <- mcols(x, use.names=FALSE)
    mcols(x) <- NULL
    x_groups <- .extract_groups_from_GenomicRanges(x)
    x_ranges <- IRanges:::.shift_ranges_in_groups_to_first_circle(ranges(x),
                                   x_groups, .get_circle_length(x))
    ranges(x) <- x_ranges
    x_nclists <- IRanges:::.nclists(x_ranges, x_groups)
    new2("GNCList", nclists=x_nclists,
                    granges=x,
                    elementMetadata=ans_mcols,
                    check=FALSE)
}

setAs("GenomicRanges", "GNCList", function(from) GNCList(from))


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Subsetting
###

setMethod("extractROWS", "GNCList",
    function(x, i) as(callGeneric(as(x, "GRanges"), i), class(x))
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### findOverlaps_GNCList()
###

### NOT exported.
findOverlaps_GNCList <- function(query, subject,
             maxgap=-1L, minoverlap=0L,
             type=c("any", "start", "end", "within", "extend", "equal"),
             select=c("all", "first", "last", "arbitrary", "count"),
             ignore.strand=FALSE)
{
    if (!(is(query, "GenomicRanges") && is(subject, "GenomicRanges")))
        stop("'query' and 'subject' must be GenomicRanges objects")
    type <- match.arg(type)
    select <- match.arg(select)
    if (!isTRUEorFALSE(ignore.strand))
        stop("'ignore.strand' must be TRUE or FALSE")

    si <- merge(seqinfo(query), seqinfo(subject))
    q_seqlevels <- seqlevels(query)
    s_seqlevels <- seqlevels(subject)
    common_seqlevels <- intersect(q_seqlevels, s_seqlevels)
    NG <- length(common_seqlevels)
    q_group_idx <- match(common_seqlevels, q_seqlevels)  # of length NG
    s_group_idx <- match(common_seqlevels, s_seqlevels)  # of length NG

    ## Extract 'q_groups' and 's_groups' (both of length NG).
    q_groups <- .extract_groups_from_GenomicRanges(query)[q_group_idx]
    s_groups <- .extract_groups_from_GenomicRanges(subject)[s_group_idx]

    ## Extract 'nclists' and 'nclist_is_q' (both of length NG).
    if (is(subject, "GNCList")) {
        nclists <- subject@nclists[s_group_idx]
        nclist_is_q <- rep.int(FALSE, NG)
    } else if (is(query, "GNCList")) {
        nclists <- query@nclists[q_group_idx]
        nclist_is_q <- rep.int(TRUE, NG)
    } else {
        ## We'll do "on-the-fly preprocessing".
        nclists <- vector(mode="list", length=NG)
        nclist_is_q <- rep.int(NA, NG)
    }

    ## Extract 'circle_length' (of length NG).
    circle_length <- .get_circle_length(si)[q_group_idx]

    ## Extract 'q_space' and 's_space'.
    if (ignore.strand) {
        q_space <- s_space <- NULL
    } else {
        q_space <- as.integer(strand(query)) - 3L
        s_space <- as.integer(strand(subject)) - 3L
    }

    ## GO!
    IRanges:::find_overlaps_in_groups_NCList(
                          ranges(query), q_space, q_groups,
                          ranges(subject), s_space, s_groups,
                          nclists, nclist_is_q,
                          maxgap, minoverlap, type, select, circle_length)
}

