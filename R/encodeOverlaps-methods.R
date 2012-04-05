### =========================================================================
### encodeOverlaps methods and related utilities
### -------------------------------------------------------------------------


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### 2 non-exported utilities for safe extraction of group seqnames and strand
### of a GRangesList object.
###

.groupSeqnames <- function(x, errmsg)
{
    if (!is(x, "GRangesList"))
        stop("'x' must be a GRangesList object")
    group_seqnames <- runValue(seqnames(x))
    elt_lens <- elementLengths(group_seqnames)
    if (!all(elt_lens == 1L))
        stop(errmsg)
    unlist(group_seqnames, use.names=FALSE)
}

.groupStrand <- function(x, errmsg)
{
    if (!is(x, "GRangesList"))
        stop("'x' must be a GRangesList object")
    group_strand <- runValue(strand(x))
    elt_lens <- elementLengths(group_strand)
    if (!all(elt_lens == 1L))
        stop(errmsg)
    unlist(group_strand, use.names=FALSE)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### 2 non-exported utilities for inverting the strand of an object.
###
### TODO: We should probably have an invertStrand() generic with methods for
### GRanges, GRangesList, GappedAlignments, GappedAlignmentPairs, and possibly
### more, instead of this.

### Works on GRanges and GappedAlignments objects. More generally, it should
### work on any object that has: (1) a strand() getter that returns a
### 'factor'-Rle, and (2) a strand() setter.
invertRleStrand <- function(x)
{
    x_strand <- strand(x)
    runValue(x_strand) <- strand(runValue(x_strand) == "+")
    strand(x) <- x_strand
    x
}

invertRleListStrand <- function(x)
{
    x@unlistData <- invertRleStrand(x@unlistData)
    x
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### flipQuery()
###

flipQuery <- function(x, i)
{
    if (!is(x, "GRangesList"))
        stop("'x' must be a GRangesList object")
    if (missing(i))
        i <- seq_len(length(x))
    else
        i <- IRanges:::normalizeSingleBracketSubscript(i, x)
    xi <- x[i]
    x[i] <- invertRleListStrand(revElements(xi))
    xi_query.break <- elementMetadata(xi)$query.break
    if (!is.null(xi_query.break)) {
        revxi_query.break <- elementLengths(xi) - xi_query.break
        elementMetadata(x)$query.break[i] <- revxi_query.break
    }
    x
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Should we use a generic + methods for this?
###

.get_GRanges_spaces <- function(x)
{
        ans <- as.integer(seqnames(x))
        x_strand <- as.integer(strand(x))
        ans <- ans * 3L + x_strand
        is_minus <- which(x_strand == as.integer(strand("-")))
        ans[is_minus] <- - ans[is_minus]
        ans
}

.get_GRangesList_spaces <- function(x)
{
        unlisted_ans <- .get_GRanges_spaces(x@unlistData)
        as.list(relist(unlisted_ans, x))
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### "encodeOverlaps" methods.
###

.flipQueryIfWrongStrand <- function(query, subject)
{
    if (!is(query, "GRangesList") || !is(subject, "GRangesList"))
        stop("'query' and 'subject' must be GRangesList objects")

    ## Checking the query and subject seqnames.
    errmsg <- c("some alignments in 'query' have ranges on ",
                "more than 1 reference sequence (fusion reads?)")
    query_seqnames <- .groupSeqnames(query, errmsg)
    errmsg <- c("some transcripts in 'subject' mix exons from ",
                "different chromosomes (trans-splicing?)")
    subject_seqnames <- .groupSeqnames(subject, errmsg)
    ## Should never happen if 'encodeOverlaps(query, subject, hits)'
    ## was called with 'hits' being the result of a call to
    ## 'findOverlaps(query, subject)'.
    if (!all(query_seqnames == subject_seqnames))
        stop("cannot use 'flip.query.if.wrong.strand=TRUE' to ",
             "encode overlaps across chromosomes")

    ## Checking the query and subject strand.
    errmsg <- c("some alignments in 'query' have ranges on ",
                "both strands")
    query_strand <- .groupStrand(query, errmsg)
    errmsg <- c("some transcripts in 'subject' mix exons from ",
                "both strands (trans-splicing?)")
    subject_strand <- .groupStrand(subject, errmsg)

    ## Flip queries to put them on same strand as subjects.
    flip_idx <- which(query_strand != subject_strand)
    flipQuery(query, flip_idx)
}

.GRangesList_encodeOverlaps <- function(query, subject,
                                        flip.query.if.wrong.strand=FALSE)
{
    if (!isTRUEorFALSE(flip.query.if.wrong.strand))
        stop("'flip.query.if.wrong.strand' must be TRUE or FALSE")
    seqinfo <- merge(seqinfo(query), seqinfo(subject))
    seqlevels(query) <- seqlevels(subject) <- seqlevels(seqinfo)
    if (flip.query.if.wrong.strand)
        query <- .flipQueryIfWrongStrand(query, subject)
    query.breaks <- elementMetadata(query)$query.break
    RangesList_encodeOverlaps(as.list(start(query)),
                              as.list(width(query)),
                              as.list(start(subject)),
                              as.list(width(subject)),
                              query.spaces=.get_GRangesList_spaces(query),
                              subject.spaces=.get_GRangesList_spaces(subject),
                              query.breaks=query.breaks)
}

setMethod("encodeOverlaps", c("GRangesList", "GRangesList", "missing"),
    function(query, subject, hits=NULL, flip.query.if.wrong.strand=FALSE)
        .GRangesList_encodeOverlaps(query, subject,
                flip.query.if.wrong.strand=flip.query.if.wrong.strand)
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### selectEncodingWithCompatibleStrand().
###

selectEncodingWithCompatibleStrand <- function(x, y,
                                               query.strand, subject.strand,
                                               hits=NULL)
{
    if (!is(x, "OverlapEncodings"))
        stop("'x' must be an OverlapEncodings object")
    if (!is(y, "OverlapEncodings"))
        stop("'y' must be an OverlapEncodings object")
    if (!is.null(hits)) {
        if (!is(hits, "Hits"))
            stop("'hits' must be a Hits object or NULL")
        query.strand <- query.strand[queryHits(hits)]
        subject.strand <- subject.strand[subjectHits(hits)]
    }
    ans <- x
    names(ans) <- NULL
    elementMetadata(ans) <- NULL
    idx <- which(query.strand != subject.strand)
    ans@Loffset[idx] <- y@Loffset[idx]
    ans@Roffset[idx] <- y@Roffset[idx]
    ans_encoding <- as.character(ans@encoding)
    ans_encoding[idx] <- as.character(y@encoding[idx])
    ans@encoding <- as.factor(ans_encoding)
    ans
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### isCompatibleWithSplicing().
###

.extract_ngap_from_encoding <- function(x)
{
    as.integer(unlist(strsplit(sub(":.*", "", x), "--", fixed=TRUE),
                      use.names=FALSE)) - 1L
}

.check_ngap_max <- function(x)
{
    ngap <- .extract_ngap_from_encoding(x)
    if (length(ngap) != 0L && max(ngap) > 3L)
        stop("reads with more than 3 gaps are not supported yet, sorry")
}

setGeneric("isCompatibleWithSplicing",
    function(x) standardGeneric("isCompatibleWithSplicing")
)

.REGEX_BLOCKS1 <- "[fgij]"
.REGEX_BLOCKS2 <- c("[jg][^:-]",
                "[^:-][gf]")
.REGEX_BLOCKS3 <- c("[jg][^:-][^:-]",
                "[^:-]g[^:-]",
                "[^:-][^:-][gf]")
.REGEX_BLOCKS4 <- c("[jg][^:-][^:-][^:-]",
                "[^:-]g[^:-][^:-]",
                "[^:-][^:-]g[^:-]",
                "[^:-][^:-][^:-][gf]")

.get_CompatibleWithSplicing_regex <- function()
{
    ## Sub-regex for single-end reads.
    Ssubregex1 <- .REGEX_BLOCKS1
    Ssubregex2 <- paste0(.REGEX_BLOCKS2, collapse=":")
    Ssubregex3 <- paste0(.REGEX_BLOCKS3, collapse=":")
    Ssubregex4 <- paste0(.REGEX_BLOCKS4, collapse=":")
    Ssubregex <- paste(Ssubregex1, Ssubregex2, Ssubregex3, Ssubregex4, sep="|")
    Ssubregex <- paste0(":(", Ssubregex, "):")

    ## Sub-regex for paired-end reads.
    Rencoding <- "-[^:-]*" 
    Lsubregex1 <- paste0(":", .REGEX_BLOCKS1, "-")
    Lsubregex2 <- paste0(":", .REGEX_BLOCKS2, "-", collapse=Rencoding)
    Lsubregex3 <- paste0(":", .REGEX_BLOCKS3, "-", collapse=Rencoding)
    Lsubregex4 <- paste0(":", .REGEX_BLOCKS4, "-", collapse=Rencoding)
    Lsubregex <- paste(Lsubregex1, Lsubregex2, Lsubregex3, Lsubregex4, sep="|")
    Lsubregex <- paste0("(", Lsubregex, ")")

    Lencoding <- "[^:-]*-" 
    Rsubregex1 <- paste0("-", .REGEX_BLOCKS1, ":")
    Rsubregex2 <- paste0("-", .REGEX_BLOCKS2, ":", collapse=Lencoding)
    Rsubregex3 <- paste0("-", .REGEX_BLOCKS3, ":", collapse=Lencoding)
    Rsubregex4 <- paste0("-", .REGEX_BLOCKS4, ":", collapse=Lencoding)
    Rsubregex <- paste(Rsubregex1, Rsubregex2, Rsubregex3, Rsubregex4, sep="|")
    Rsubregex <- paste0("(", Rsubregex, ")")

    LRsubregex <- paste0(Lsubregex, ".*", Rsubregex)

    ## Final regex.
    paste0("(", Ssubregex, "|", LRsubregex, ")")
}

.isCompatibleWithSplicing <- function(x)
{
    if (!is.character(x))
        stop("'x' must be a character vector")
    .check_ngap_max(x)
    grepl(.get_CompatibleWithSplicing_regex(), x)
}

.whichCompatibleWithSplicing <- function(x)
{
    if (!is.character(x))
        stop("'x' must be a character vector")
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

.get_CompatibleWithSkippedExons_regex <- function(max.skipped.exons=NA)
{
    if (!identical(max.skipped.exons, NA))
        stop("only 'max.skipped.exons=NA' is supported for now, sorry")

    ## Sub-regex for single-end reads.
    Ssubregex1 <- .REGEX_BLOCKS1
    Ssubregex2 <- paste0(.REGEX_BLOCKS2, collapse=":(..:)*")
    Ssubregex3 <- paste0(.REGEX_BLOCKS3, collapse=":(...:)*")
    Ssubregex4 <- paste0(.REGEX_BLOCKS4, collapse=":(....:)*")
    Ssubregex <- paste(Ssubregex1, Ssubregex2, Ssubregex3, Ssubregex4, sep="|")
    Ssubregex <- paste0(":(", Ssubregex, "):")

    ## Sub-regex for paired-end reads.
    Lsubregex1 <- paste0(":", .REGEX_BLOCKS1, "-")
    Lsubregex2 <- paste0(":", .REGEX_BLOCKS2, "-", collapse=".*")
    Lsubregex3 <- paste0(":", .REGEX_BLOCKS3, "-", collapse=".*")
    Lsubregex4 <- paste0(":", .REGEX_BLOCKS4, "-", collapse=".*")
    Lsubregex <- paste(Lsubregex1, Lsubregex2, Lsubregex3, Lsubregex4, sep="|")
    Lsubregex <- paste0("(", Lsubregex, ")")

    Rsubregex1 <- paste0("-", .REGEX_BLOCKS1, ":")
    Rsubregex2 <- paste0("-", .REGEX_BLOCKS2, ":", collapse=".*")
    Rsubregex3 <- paste0("-", .REGEX_BLOCKS3, ":", collapse=".*")
    Rsubregex4 <- paste0("-", .REGEX_BLOCKS4, ":", collapse=".*")
    Rsubregex <- paste(Rsubregex1, Rsubregex2, Rsubregex3, Rsubregex4, sep="|")
    Rsubregex <- paste0("(", Rsubregex, ")")

    LRsubregex <- paste0(Lsubregex, ".*", Rsubregex)

    ## Final regex.
    paste0("(", Ssubregex, "|", LRsubregex, ")")
}

.isCompatibleWithSkippedExons <- function(x, max.skipped.exons=NA)
{
    if (!is.character(x))
        stop("'x' must be a character vector")
    .check_ngap_max(x)
    grepl(.get_CompatibleWithSkippedExons_regex(max.skipped.exons), x) &
    !grepl(.get_CompatibleWithSplicing_regex(), x)
}

.whichCompatibleWithSkippedExons <- function(x, max.skipped.exons=NA)
{
    if (!is.character(x))
        stop("'x' must be a character vector")
    .check_ngap_max(x)
    setdiff(grep(.get_CompatibleWithSkippedExons_regex(max.skipped.exons), x),
            grep(.get_CompatibleWithSplicing_regex(), x))
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


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### extractSpannedExonRanks().
###

setGeneric("extractSpannedExonRanks",
    function(x) standardGeneric("extractSpannedExonRanks")
)

.extractSteppedExonRanksFromEncodingBlocks <- function(encoding_blocks,
                                                       regex_blocks)
{
    regexprs <- paste0("^", regex_blocks, "$")
    ii <- lapply(regexprs, grep, encoding_blocks)
    ii_elt_lens <- elementLengths(ii)
    if (any(ii_elt_lens == 0L))
        return(integer(0))
    if (any(ii_elt_lens != 1L))
        stop("cannot unambiguously extract stepped exon ranks from ",
             "encoding \"", paste0(encoding_blocks, collapse=":"), "\"")
    ans <- unlist(ii, use.names=FALSE)
    diff_ans <- diff(ans)
    if (any(diff_ans <= 0L))
        return(integer(0))
    ans
}

### 'encoding' must be a single encoding.
.extractSteppedExonRanks <- function(encoding)
{
    encoding_blocks <- strsplit(encoding, ":", fixed=TRUE)[[1]]
    ngap <- .extract_ngap_from_encoding(encoding_blocks[1L])
    if (length(ngap) == 2L)
        stop("extractSpannedExonRanks() and extractSkippedExonRanks() ",
             "don't yet support overlap encodings of paired-end reads, sorry")
    encoding_blocks <- encoding_blocks[-1L]
    if (ngap == 0L) {
        regex_blocks <- .REGEX_BLOCKS1
    } else if (ngap == 1L) {
        regex_blocks <- .REGEX_BLOCKS2
    } else if (ngap == 2L) {
        regex_blocks <- .REGEX_BLOCKS3
    } else if (ngap == 3L) {
        regex_blocks <- .REGEX_BLOCKS4
    } else {
        stop("reads with more than 3 gaps are not supported yet, sorry")
    }
    .extractSteppedExonRanksFromEncodingBlocks(encoding_blocks, regex_blocks)
}

setMethod("extractSpannedExonRanks", "character",
    function(x)
    {
        .extractRanks <- function(encoding) {
            ranks <- .extractSteppedExonRanks(encoding)
            if (length(ranks) == 0L)
                return(c(NA_integer_, NA_integer_))
            c(ranks[1L], ranks[length(ranks)])
        }
        ranks <- lapply(x, .extractRanks)
        if (length(ranks) == 0L) {
            startRank <- endRank <- integer(0)
        } else {
            ranks <- unlist(ranks, use.names=FALSE)
            startRank <- ranks[c(TRUE, FALSE)]
            endRank <- ranks[c(FALSE, TRUE)]
        }
        data.frame(startRank=startRank, endRank=endRank,
                   check.names=FALSE, stringsAsFactors=FALSE)
    }
)

setMethod("extractSpannedExonRanks", "factor",
    function(x)
    {
        if (length(x) == 0L)
            return(list())
        ranks <- extractSpannedExonRanks(levels(x))
        ans <- ranks[as.integer(x), , drop=FALSE]
        rownames(ans) <- NULL
        ans
    }
)

setMethod("extractSpannedExonRanks", "OverlapEncodings",
    function(x)
    {
        ranks <- extractSpannedExonRanks(encoding(x))
        ranks$startRank <- ranks$startRank + Loffset(x)
        ranks$endRank <- ranks$endRank + Loffset(x)
        ranks
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### extractSkippedExonRanks().
###

setGeneric("extractSkippedExonRanks",
    function(x) standardGeneric("extractSkippedExonRanks")
)

setMethod("extractSkippedExonRanks", "character",
    function(x)
    {
        .extractRanks <- function(encoding) {
            ranks <- .extractSteppedExonRanks(encoding)
            if (length(ranks) == 0L)
                return(ranks)
            setdiff(ranks[1L]:ranks[length(ranks)], ranks)
        }
        lapply(x, .extractRanks)
    }
)

setMethod("extractSkippedExonRanks", "factor",
    function(x)
    {
        if (length(x) == 0L)
            return(list())
        ranks <- extractSkippedExonRanks(levels(x))
        ranks[as.integer(x)]
    }
)

setMethod("extractSkippedExonRanks", "OverlapEncodings",
    function(x)
    {
        ranks <- extractSkippedExonRanks(encoding(x))
        tmp <- unlist(ranks, use.names=FALSE)
        tmp <- tmp + rep.int(Loffset(x), elementLengths(ranks))
        flevels <- seq_len(length(ranks))
        f <- factor(rep.int(flevels, elementLengths(ranks)),
                    levels=flevels)
        unname(split(tmp, f))
    }
)

