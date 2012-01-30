### =========================================================================
### encodeOverlaps methods and related utilities
### -------------------------------------------------------------------------


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
### Related utilities.
###

readIsCompatibleWithSplicing <- function(x)
{
    if (!is(x, "OverlapEncodings"))
        stop("'x' must be an OverlapEncodings object")
    if (length(x) == 0L)
        return(logical(0))
    ngap <- as.integer(sub(":.*", "", levels(encoding(x)))) - 1L
    if (max(ngap) > 2)
        stop("don't support reads with more than 2 gaps for now, sorry")
    ## Regex to use for reads with 1 range (no gaps):
    pattern1 <- ":[fgij]:"
    ## Regex to use for reads with 2 ranges (1 gap):
    pattern2 <- ":[jg].:.[gf]:"
    ## Regex to use for reads with 3 ranges (2 gaps):
    pattern3 <- ":[jg]..:.g.:..[gf]:"
    pattern123 <- ":([fgij]|[jg].:.[gf]|[jg]..:.g.:..[gf]):"
    as.integer(encoding(x)) %in% grep(pattern123, levels(encoding(x)))
}

