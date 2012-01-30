### =========================================================================
### encodeOverlaps methods and related utilities
### -------------------------------------------------------------------------


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### "encodeOverlaps" methods.
###

if (FALSE) {
### Not sure we need to bother with methods that do 1-to-1 range overlap
### encodings. What would be the use cases?
setMethod("encodeOverlaps", c("GRanges", "GRanges"),
    function(query, subject)
    {
        seqinfo <- merge(seqinfo(query), seqinfo(subject))
        seqlevels(query) <- seqlevels(subject) <- seqlevels(seqinfo)
        query.space <- as.integer(seqnames(query)) * 3L +
                       as.integer(strand(query))
        subject.space <- as.integer(seqnames(subject)) * 3L +
                         as.integer(strand(subject))
        same_space <- query.space == subject.space
        ans <- character(length(same_space))
        ans[!same_space] <- "X"
        ans[same_space] <- encodeOverlaps(ranges(query)[same_space],
                                          ranges(subject)[same_space])
        ans
    }
)
}

.get_GRangesList_space <- function(x)
{
        x_seqnames <- seqnames(x)
        x_strand <- strand(x)
        seqnames <- as.integer(unlist(x_seqnames, use.names=FALSE))
        strand <- as.integer(unlist(x_strand, use.names=FALSE))
        as.list(relist(seqnames * 3L + strand, x_seqnames))
}

setMethod("encodeOverlaps", c("GRangesList", "GRangesList"),
    function(query, subject)
    {
        seqinfo <- merge(seqinfo(query), seqinfo(subject))
        seqlevels(query) <- seqlevels(subject) <- seqlevels(seqinfo)
        RangesList_encodeOverlaps(as.list(start(query)),
                                  as.list(width(query)),
                                  as.list(start(subject)),
                                  as.list(width(subject)),
                                  query.space=.get_GRangesList_space(query),
                                  subject.space=.get_GRangesList_space(subject))
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### isCompatibleWithSplicing().
###

setGeneric("isCompatibleWithSplicing",
    function(x) standardGeneric("isCompatibleWithSplicing")
)

.check_ngap_max <- function(x)
{
    ngap <- as.integer(sub(":.*", "", x)) - 1L
    if (max(ngap) > 3L)
        stop("reads with more than 3 gaps are not supported yet, sorry")
}

.get_CompatibleWithSplicing_regex <- function()
{
    part1 <- "[fgij]"                     # reads with 1 range (0 gap)
    part2 <- "[jg].:.[gf]"                # reads with 2 ranges (1 gap)
    part3 <- "[jg]..:.g.:..[gf]"          # reads with 3 ranges (2 gaps)
    part4 <- "[jg]...:.g..:..g.:...[gf]"  # reads with 4 ranges (3 gaps)
    paste0(":(", paste(part1, part2, part3, part4, sep="|"), "):")
}

.isCompatibleWithSplicing <- function(x)
{
    .check_ngap_max(x)
    grepl(.get_CompatibleWithSplicing_regex(), x)
}

.whichCompatibleWithSplicing <- function(x)
{
    .check_ngap_max(x)
    grep(.get_CompatibleWithSplicing_regex(), x)
}

setMethod("isCompatibleWithSplicing", "character", .isCompatibleWithSplicing)

setMethod("isCompatibleWithSplicing", "factor",
    function(x)
    {
        if (length(x) == 0L)
            return(logical(0))
        idx <- .whichCompatibleWithSplicing(levels(x))
        as.integer(x) %in% idx
    }
)

setMethod("isCompatibleWithSplicing", "OverlapEncodings",
    function(x)
        isCompatibleWithSplicing(encoding(x))
)

