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


### Return the "split factor" of GenomicRanges object 'x'.
.split_factor <- function(x)
{
    #x_nseqlevels <- length(seqlevels(x))
    #levels_prefix <- rep(seq_len(x_nseqlevels), each=3L)
    #levels_suffix <- rep.int(levels(strand()), x_nseqlevels)
    #factor(paste0(as.integer(seqnames(x)), as.character(strand(x))),
    #       levels=paste0(levels_prefix, levels_suffix))
    seqnames(x)
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
        ans_rglist <- split(ranges(from_granges), .split_factor(from_granges))
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
    rglist <- split(ranges(x), .split_factor(x))
    ans_nclists <- NCLists(rglist)@nclists
    new2("GNCList", nclists=ans_nclists, granges=x, check=FALSE)
}

setAs("GenomicRanges", "GNCList", function(from) GNCList(from))


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### combine_and_select_hits()
###
### A low-level (non exported) utility function for post-processing a list
### of Hits objects. Used by findOverlaps_GNCList() below and "findOverlaps"
### method for GenomicRanges objects.
###

### Rearrange the seqlevels in 'x' so that its first N seqlevels are exactly
### 'seqlevels' (in the same order).
.align_seqlevels <- function(x, seqlevels)
{
    seqlevels(x) <- union(seqlevels, seqlevels(x))
    x
}

.remap_and_combine_Hits_objects <- function(objects,
                                            q_split_factor, s_split_factor)
{
    query_len <- length(q_split_factor)
    subject_len <- length(s_split_factor)
    if (length(objects) == 0L) {
        q_hits <- s_hits <- integer(0)
    } else {
        ## Compute 'q_hits'.
        query_maps <- split(seq_len(query_len), q_split_factor)
        q_hits <-
            lapply(seq_along(objects),
                   function(i) query_maps[[i]][queryHits(objects[[i]])])
        q_hits <- unlist(q_hits)
        ## Compute 's_hits'.
        subject_maps <- split(seq_len(subject_len), s_split_factor)
        s_hits <-
            lapply(seq_along(objects),
                   function(i) subject_maps[[i]][subjectHits(objects[[i]])])
        s_hits <- unlist(s_hits)
        ## Order 'q_hits' and 's_hits'.
        oo <- S4Vectors:::orderIntegerPairs(q_hits, s_hits)
        q_hits <- q_hits[oo]
        s_hits <- s_hits[oo]
    }
    new2("Hits", queryHits=q_hits, subjectHits=s_hits,
                 queryLength=query_len, subjectLength=subject_len,
                 check=FALSE)
}

.strand_is_incompatible <- function(qh_strand, sh_strand)
{
    qh_strand == "+" & sh_strand == "-" | qh_strand == "-" & sh_strand == "+"
}

.drop_overlaps_with_incompatible_strand <- function(hits, q_strand, s_strand)
{
    qh_strand <- q_strand[queryHits(hits)]
    sh_strand <- s_strand[subjectHits(hits)]
    drop_idx <- which(.strand_is_incompatible(qh_strand, sh_strand))
    if (length(drop_idx) != 0L)
        hits <- hits[-drop_idx]
    hits
}

combine_and_select_hits <- function(list_of_Hits,
                                    q_split_factor, s_split_factor,
                                    q_strand, s_strand,
                                    select, ignore.strand)
{
    hits <- .remap_and_combine_Hits_objects(list_of_Hits,
                                            q_split_factor, s_split_factor)
    if (!ignore.strand)
        hits <- .drop_overlaps_with_incompatible_strand(hits,
                                                        q_strand, s_strand)
    selectHits(hits, select=select)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### findOverlaps_GNCList()
###

### NOT exported.
findOverlaps_GNCList <- function(query, subject, min.score=1L,
                                 type=c("any", "start", "end",
                                        "within", "extend", "equal"),
                                 select=c("all", "first", "last", "arbitrary"),
                                 ignore.strand=FALSE)
{
    if (!(is(query, "GNCList") || is(subject, "GNCList")))
        stop("'query' or 'subject' must be a GNCList object")
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

    if (is(subject, "GNCList")) {
       si <- merge(seqinfo(subject), seqinfo(query))
       query <- .align_seqlevels(query, seqlevels(subject))
       q_split_factor <- .split_factor(query)
       s_split_factor <- .split_factor(subject)
       q_rglist <- split(ranges(query), q_split_factor)
       s_rglist <- as(subject, "NCLists")
    } else {
       si <- merge(seqinfo(query), seqinfo(subject))
       subject <- .align_seqlevels(subject, seqlevels(query))
       q_split_factor <- .split_factor(query)
       s_split_factor <- .split_factor(subject)
       q_rglist <- as(query, "NCLists")
       s_rglist <- split(ranges(subject), s_split_factor)
    }
    circle_length <- seqlengths(si)
    circle_length[!(isCircular(si) %in% TRUE)] <- NA_integer_
    circle_length <- head(circle_length, n=min(length(q_rglist),
                                               length(s_rglist)))
    list_of_Hits <- IRanges:::findOverlaps_NCLists(q_rglist, s_rglist,
                                  min.score=min.score, type=type,
                                  select="all",
                                  circle.length=circle_length)
    combine_and_select_hits(list_of_Hits,
                            q_split_factor, s_split_factor,
                            strand(query), strand(subject),
                            select, ignore.strand)
}

