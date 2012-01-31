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

.check_ngap_max <- function(x)
{
    ngap <- as.integer(sub(":.*", "", x)) - 1L
    if (max(ngap) > 3L)
        stop("reads with more than 3 gaps are not supported yet, sorry")
}

setGeneric("isCompatibleWithSplicing",
    function(x) standardGeneric("isCompatibleWithSplicing")
)

.get_CompatibleWithSplicing_regex <- function()
{
    subregex1 <- "[fgij]"                     # reads with 1 range (0 gap)
    subregex2 <- "[jg].:.[gf]"                # reads with 2 ranges (1 gap)
    subregex3 <- "[jg]..:.g.:..[gf]"          # reads with 3 ranges (2 gaps)
    subregex4 <- "[jg]...:.g..:..g.:...[gf]"  # reads with 4 ranges (3 gaps)
    paste0(":(",
           paste(subregex1, subregex2, subregex3, subregex4, sep="|"),
           "):")
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
    function(x) isCompatibleWithSplicing(encoding(x))
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### isCompatibleWithSkippedExons().
###

setGeneric("isCompatibleWithSkippedExons", signature="x",
    function(x, max.skipped.exons=NA)
        standardGeneric("isCompatibleWithSkippedExons")
)

.get_CompatibleWithSkippedExons_regex <- function(max.skipped.exons)
{
    if (!identical(max.skipped.exons, NA))
        stop("only 'max.skipped.exons=NA' is supported for now, sorry")
    ## Reads with 2 ranges (1 gap):
    subregex2 <- "[jg].:(..:)+.[gf]"
    ## Reads with 3 ranges (2 gaps):
    subregex3 <- "[jg]..:((...:)+.g.:(...:)*|(...:)*.g.:(...:)+)..[gf]"
    ## Reads with 4 ranges (3 gaps):
    subregex4 <- "[jg]...:((....:)+.g..:(....:)*..g.:(....:)*|(....:)*.g..:(....:)+..g.:(....:)*|(....:)*.g..:(....:)*..g.:(....:)+)...[gf]"
    paste0(":(",
           paste(subregex2, subregex3, subregex4, sep="|"),
           "):")
}

.isCompatibleWithSkippedExons <- function(x, max.skipped.exons=NA)
{
    .check_ngap_max(x)
    grepl(.get_CompatibleWithSkippedExons_regex(max.skipped.exons), x)
}

.whichCompatibleWithSkippedExons <- function(x, max.skipped.exons=NA)
{
    .check_ngap_max(x)
    grep(.get_CompatibleWithSkippedExons_regex(max.skipped.exons), x)
}

setMethod("isCompatibleWithSkippedExons", "character",
    .isCompatibleWithSkippedExons
)

setMethod("isCompatibleWithSkippedExons", "factor",
    function(x, max.skipped.exons=NA)
    {
        if (length(x) == 0L)
            return(logical(0))
        idx <- .whichCompatibleWithSkippedExons(levels(x),
                        max.skipped.exons=max.skipped.exons)
        as.integer(x) %in% idx
    }
)

setMethod("isCompatibleWithSkippedExons", "OverlapEncodings",
    function(x, max.skipped.exons=NA)
        isCompatibleWithSkippedExons(encoding(x),
                        max.skipped.exons=max.skipped.exons)
)

setGeneric("extractSkippedExonRanks",
    function(x) standardGeneric("extractSkippedExonRanks")
)

setMethod("extractSkippedExonRanks", "character",
    function(x)
    {
        xx <- strsplit(x, ":", fixed=TRUE)
        lapply(xx,
               function(s0) {
                   ngap <- as.integer(s0[1L]) - 1L
                   ans <- integer(0)
                   if (ngap <= 0L)
                       return(ans)
                   s1 <- s0[-1L]
                   if (length(s1) < ngap + 2L)
                       return(ans)
                   if (ngap == 1L) {
                       if (grepl("^[jg].$", s1[1L])
                        && grepl("^.[gf]$", s1[length(s1)])) {
                           ans <- 2L:(length(s1)-1L)
                       }
                   } else if (ngap == 2L) {
                       if (grepl("^[jg]..$", s1[1L])
                        && grepl("^..[gf]$", s1[length(s1)])) {
                           s2 <- s1[-c(1L, length(s1))]
                           i <- grep("^.g.$", s2)
                           if (length(i) == 0L)
                               return(ans)
                           if (length(i) != 1L)
                               stop("unexpected/unsupported overlap ",
                                    "situation: length(i) > 1L")
                           ans <- (2:(length(s1)-1L))[-i]
                       }
                   } else if (ngap == 3L) {
                       if (grepl("^[jg]...$", s1[1L])
                        && grepl("^...[gf]$", s1[length(s1)])) {
                           s2 <- s1[-c(1L, length(s1))]
                           i <- grep("^.g..$", s2)
                           j <- grep("^..g.$", s2)
                           if (length(i) == 0L || length(j) == 0L)
                               return(ans)
                           if (length(i) != 1L || length(j) != 1L)
                               stop("unexpected/unsupported overlap ",
                                    "situation: length(i) or length(j) > 1L")
                           if (j <= i)
                               return(ans)
                           ans <- (2:(length(s1)-1L))[-c(i,j)]
                       }
                   } else {
                       stop("reads with more than 3 gaps ",
                            "are not supported yet, sorry")
                   }
                   ans
               })
    }
)

setMethod("extractSkippedExonRanks", "factor",
    function(x)
    {
        if (length(x) == 0L)
            return(list())
        skipped_ranks <- extractSkippedExonRanks(levels(x))
        skipped_ranks[as.integer(x)]
    }
)

setMethod("extractSkippedExonRanks", "OverlapEncodings",
    function(x)
    {
        skipped_ranks <- extractSkippedExonRanks(encoding(x))
        tmp <- unlist(skipped_ranks, use.names=FALSE)
        tmp <- tmp + rep.int(Loffset(x), elementLengths(skipped_ranks))
        flevels <- seq_len(length(skipped_ranks))
        f <- factor(rep.int(flevels, elementLengths(skipped_ranks)),
                    levels=flevels)
        unname(split(tmp, f))
    }
)

