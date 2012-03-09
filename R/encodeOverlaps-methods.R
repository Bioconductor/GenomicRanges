### =========================================================================
### encodeOverlaps methods and related utilities
### -------------------------------------------------------------------------


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### "encodeOverlaps" methods.
###

if (FALSE) {
### Not sure we need to bother with methods that do 1-to-1 range overlap
### encodings. What would be the use cases?
setMethod("encodeOverlaps", c("GRanges", "GRanges", "missing"),
    function(query, subject, hits=NULL)
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

GRangesList_encodeOverlaps <- function(query, subject, Lquery.lengths=NULL)
{
    seqinfo <- merge(seqinfo(query), seqinfo(subject))
    seqlevels(query) <- seqlevels(subject) <- seqlevels(seqinfo)
    RangesList_encodeOverlaps(as.list(start(query)),
                              as.list(width(query)),
                              as.list(start(subject)),
                              as.list(width(subject)),
                              query.space=.get_GRangesList_space(query),
                              subject.space=.get_GRangesList_space(subject),
                              Lquery.lengths=Lquery.lengths)
}

setMethod("encodeOverlaps", c("GRangesList", "GRangesList", "missing"),
    function(query, subject, hits=NULL)
        GRangesList_encodeOverlaps(query, subject)
)

setMethod("encodeOverlaps", c("GappedAlignments", "GRangesList", "missing"),
    function(query, subject, hits=NULL)
        GRangesList_encodeOverlaps(as(query, "GRangesList"), subject)
)

setMethod("encodeOverlaps", c("GappedAlignmentPairs", "GRangesList", "missing"),
    function(query, subject, hits=NULL)
    {
        Lquery.lengths <- 1L + ngap(left(query))
        GRangesList_encodeOverlaps(as(query, "GRangesList"), subject,
                                      Lquery.lengths=Lquery.lengths)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### isCompatibleWithSplicing().
###

.check_ngap_max <- function(x)
{
    ngap <- as.integer(unlist(strsplit(sub(":.*", "", x), "--", fixed=TRUE),
                              use.names=FALSE)) - 1L
    if (max(ngap) > 3L)
        stop("reads with more than 3 gaps are not supported yet, sorry")
}

setGeneric("isCompatibleWithSplicing",
    function(x) standardGeneric("isCompatibleWithSplicing")
)

.get_CompatibleWithSplicing_regex <- function()
{
    subregex1 <- "[fgij]"
    subregex2 <- c("[jg].", ".[gf]")
    subregex3 <- c("[jg]..", ".g.", "..[gf]")
    subregex4 <- c("[jg]...", ".g..", "..g.", "...[gf]")

    Ssubregex1 <- paste0(subregex1, collapse=":")
    Ssubregex2 <- paste0(subregex2, collapse=":")
    Ssubregex3 <- paste0(subregex3, collapse=":")
    Ssubregex4 <- paste0(subregex4, collapse=":")
    Ssubregex <- paste(Ssubregex1, Ssubregex2, Ssubregex3, Ssubregex4, sep="|")
    Ssubregex <- paste0(":(", Ssubregex, "):")

    Rencoding <- "-[^-:]*" 
    Lsubregex1 <- paste0(":", subregex1, "-", collapse=Rencoding)
    Lsubregex2 <- paste0(":", subregex2, "-", collapse=Rencoding)
    Lsubregex3 <- paste0(":", subregex3, "-", collapse=Rencoding)
    Lsubregex4 <- paste0(":", subregex4, "-", collapse=Rencoding)
    Lsubregex <- paste(Lsubregex1, Lsubregex2, Lsubregex3, Lsubregex4, sep="|")
    Lsubregex <- paste0("(", Lsubregex, ")")

    Lencoding <- "[^-:]*-" 
    Rsubregex1 <- paste0("-", subregex1, ":", collapse=Lencoding)
    Rsubregex2 <- paste0("-", subregex2, ":", collapse=Lencoding)
    Rsubregex3 <- paste0("-", subregex3, ":", collapse=Lencoding)
    Rsubregex4 <- paste0("-", subregex4, ":", collapse=Lencoding)
    Rsubregex <- paste(Rsubregex1, Rsubregex2, Rsubregex3, Rsubregex4, sep="|")
    Rsubregex <- paste0("(", Rsubregex, ")")

    LRsubregex <- paste0(Lsubregex, ".*", Rsubregex)

    paste0("(", Ssubregex, "|", LRsubregex, ")")
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
### FIXME: This currently does NOT work with paired-end encodings!
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

