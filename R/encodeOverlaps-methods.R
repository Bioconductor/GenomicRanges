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

.isWrongStrand <- function(query, subject)
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

    query_strand != subject_strand
}

.GRangesList_encodeOverlaps <- function(query, subject,
                                        flip.query.if.wrong.strand=FALSE)
{
    if (!isTRUEorFALSE(flip.query.if.wrong.strand))
        stop("'flip.query.if.wrong.strand' must be TRUE or FALSE")
    seqinfo <- merge(seqinfo(query), seqinfo(subject))
    seqlevels(query) <- seqlevels(subject) <- seqlevels(seqinfo)
    if (flip.query.if.wrong.strand) {
        is_wrong_strand <- .isWrongStrand(query, subject)
        query <- flipQuery(query, is_wrong_strand)
    }
    query.breaks <- elementMetadata(query)$query.break
    ans <- RangesList_encodeOverlaps(
                           as.list(start(query)),
                           as.list(width(query)),
                           as.list(start(subject)),
                           as.list(width(subject)),
                           query.spaces=.get_GRangesList_spaces(query),
                           subject.spaces=.get_GRangesList_spaces(subject),
                           query.breaks=query.breaks)
    if (flip.query.if.wrong.strand)
        ans@flippedQuery <- is_wrong_strand
    ans
}

setMethod("encodeOverlaps", c("GRangesList", "GRangesList", "missing"),
    function(query, subject, hits=NULL, flip.query.if.wrong.strand=FALSE)
        .GRangesList_encodeOverlaps(query, subject,
                flip.query.if.wrong.strand=flip.query.if.wrong.strand)
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### selectEncodingWithCompatibleStrand().
###

selectEncodingWithCompatibleStrand <- function(ovencA, ovencB,
                                               query.strand, subject.strand,
                                               hits=NULL)
{
    if (!is(ovencA, "OverlapEncodings"))
        stop("'ovencA' must be an OverlapEncodings object")
    if (!is(ovencB, "OverlapEncodings"))
        stop("'ovencB' must be an OverlapEncodings object")
    if (!is.null(hits)) {
        if (!is(hits, "Hits"))
            stop("'hits' must be a Hits object or NULL")
        query.strand <- query.strand[queryHits(hits)]
        subject.strand <- subject.strand[subjectHits(hits)]
    }
    ans <- ovencA
    names(ans) <- NULL
    elementMetadata(ans) <- NULL
    is_wrong_strand <- query.strand != subject.strand
    idx <- which(is_wrong_strand)
    ans@Loffset[idx] <- ovencB@Loffset[idx]
    ans@Roffset[idx] <- ovencB@Roffset[idx]
    ans_encoding <- as.character(ans@encoding)
    ans_encoding[idx] <- as.character(ovencB@encoding[idx])
    ans@encoding <- as.factor(ans_encoding)
    ans@flippedQuery[is_wrong_strand] <- TRUE
    ans
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### isCompatibleWithSplicing().
###

.build_encoding_patterns <- function(ngap)
{
    ## Each "atom" must match exactly 1 code in the encoding.
    ATOM0 <- "[fgij]"
    if (ngap == 0L)
        return(ATOM0)
    LEFT_ATOM <- "[jg]"
    RIGHT_ATOM <- "[gf]"
    WILDCARD_ATOM <- "[^:-]"
    sapply(seq_len(ngap + 1L),
           function(i) {
               if (i == 1L) {
                   atoms <- c(LEFT_ATOM, rep.int(WILDCARD_ATOM, ngap))
               } else if (i == ngap + 1L) {
                   atoms <- c(rep.int(WILDCARD_ATOM, ngap), RIGHT_ATOM)
               } else {
                   atoms <- c(rep.int(WILDCARD_ATOM, i-1L),
                              "g",
                              rep.int(WILDCARD_ATOM, ngap-i+1L))
               }
               paste0(atoms, collapse="")
           })
}

## Putting an arbitrary limit on the number of gaps in a read.
.NGAP_MAX <- 3L

.extract_ngap_from_encoding <- function(x)
{
    as.integer(unlist(strsplit(sub(":.*", "", x), "--", fixed=TRUE),
                      use.names=FALSE)) - 1L
}

.check_ngap_max <- function(x, max.ngap)
{
    ngap <- .extract_ngap_from_encoding(x)
    if (length(ngap) != 0L && max(ngap) > max.ngap)
        stop("reads with more than ", max.ngap,
             " gaps are not supported, sorry")
}

setGeneric("isCompatibleWithSplicing",
    function(x) standardGeneric("isCompatibleWithSplicing")
)

.build_CompatibleWithSplicing_pattern <- function(max.ngap)
{
    ## Subpattern for single-end reads.
    Ssubpattern <- sapply(0:max.ngap,
                     function(ngap)
                       paste0(.build_encoding_patterns(ngap),
                              collapse=":"))
    Ssubpattern <- paste0(":(", paste0(Ssubpattern, collapse="|"), "):")

    ## Subpattern for paired-end reads.
    Lsubpattern <- sapply(0:max.ngap,
                     function(ngap)
                       paste0(":", .build_encoding_patterns(ngap), "-",
                              collapse="-[^:-]*"))
    Lsubpattern <- paste0("(", paste0(Lsubpattern, collapse="|"), ")")

    Rsubpattern <- sapply(0:max.ngap,
                     function(ngap)
                       paste0("-", .build_encoding_patterns(ngap), ":",
                              collapse="[^:-]*-"))
    Rsubpattern <- paste0("(", paste0(Rsubpattern, collapse="|"), ")")

    LRsubpattern <- paste0(Lsubpattern, ".*", Rsubpattern)

    ## Final pattern.
    paste0("(", Ssubpattern, "|", LRsubpattern, ")")
}

.isCompatibleWithSplicing <- function(x)
{
    if (!is.character(x))
        stop("'x' must be a character vector")
    .check_ngap_max(x, .NGAP_MAX)
    grepl(.build_CompatibleWithSplicing_pattern(.NGAP_MAX), x)
}

.whichCompatibleWithSplicing <- function(x)
{
    if (!is.character(x))
        stop("'x' must be a character vector")
    .check_ngap_max(x, .NGAP_MAX)
    grep(.build_CompatibleWithSplicing_pattern(.NGAP_MAX), x)
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

.build_CompatibleWithSkippedExons_pattern <- function(max.ngap,
                                                      max.skipped.exons=NA)
{
    if (!identical(max.skipped.exons, NA))
        stop("only 'max.skipped.exons=NA' is supported for now, sorry")

    ## Subpattern for single-end reads.
    skipped_exons_subpatterns <- c(":(.:)*", ":(..:)*",
                                   ":(...:)*", ":(....:)*")
    Ssubpattern <- sapply(0:max.ngap,
                     function(ngap)
                       paste0(.build_encoding_patterns(ngap),
                              collapse=skipped_exons_subpatterns[ngap+1L]))
    Ssubpattern <- paste0(":(", paste0(Ssubpattern, collapse="|"), "):")

    ## Subpattern for paired-end reads.
    Lsubpattern <- sapply(0:max.ngap,
                     function(ngap)
                       paste0(":", .build_encoding_patterns(ngap), "-",
                              collapse=".*"))
    Lsubpattern <- paste0("(", paste0(Lsubpattern, collapse="|"), ")")

    Rsubpattern <- sapply(0:max.ngap,
                     function(ngap)
                       paste0("-", .build_encoding_patterns(ngap), ":",
                              collapse=".*"))
    Rsubpattern <- paste0("(", paste0(Rsubpattern, collapse="|"), ")")

    LRsubpattern <- paste0(Lsubpattern, ".*", Rsubpattern)

    ## Final pattern.
    paste0("(", Ssubpattern, "|", LRsubpattern, ")")
}

.isCompatibleWithSkippedExons <- function(x, max.skipped.exons=NA)
{
    if (!is.character(x))
        stop("'x' must be a character vector")
    .check_ngap_max(x, .NGAP_MAX)
    pattern1 <- .build_CompatibleWithSkippedExons_pattern(.NGAP_MAX,
                                                          max.skipped.exons)
    pattern2 <- .build_CompatibleWithSplicing_pattern(.NGAP_MAX)
    grepl(pattern1, x) & !grepl(pattern2, x)
}

.whichCompatibleWithSkippedExons <- function(x, max.skipped.exons=NA)
{
    if (!is.character(x))
        stop("'x' must be a character vector")
    .check_ngap_max(x, .NGAP_MAX)
    pattern1 <- .build_CompatibleWithSkippedExons_pattern(.NGAP_MAX,
                                                          max.skipped.exons)
    pattern2 <- .build_CompatibleWithSplicing_pattern(.NGAP_MAX)
    setdiff(grep(pattern1, x), grep(pattern2, x))
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
### extractSteppedExonRanks().
###

.extractSteppedExonRanksFromEncodingBlocks <- function(encoding_blocks,
                                                       encoding_patterns)
{
    patterns <- paste0("^", encoding_patterns, "$")
    ii <- lapply(patterns, grep, encoding_blocks)
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
### Returns a sorted integer vector, unnamed and strictly sorted if single-end
### read, named and not necessarily strictly sorted if paired-end read (last
### exon stepped by the left end can be the same as first exon stepped by right
### end).
.extractSteppedExonRanks <- function(encoding, for.query.right.end=FALSE)
{
    if (!isTRUEorFALSE(for.query.right.end))
        stop("'for.query.right.end' must be TRUE or FALSE")
    encoding_blocks <- strsplit(encoding, ":", fixed=TRUE)[[1L]]
    ngap <- .extract_ngap_from_encoding(encoding_blocks[1L])
    encoding_blocks <- encoding_blocks[-1L]
    if (length(ngap) == 1L) {
        ## Single-end read.
        if (for.query.right.end)
            stop("cannot use 'for.query.right.end=TRUE' ",
                 "on single-end encoding: ", encoding)
        encoding_patterns <- .build_encoding_patterns(ngap)
        return(.extractSteppedExonRanksFromEncodingBlocks(encoding_blocks,
                                                          encoding_patterns))
    }
    if (length(ngap) != 2L)  # should never happen
        stop(encoding, ": invalid encoding")
    ## Paired-end read.
    encoding_blocks <- strsplit(encoding_blocks, "--", fixed=TRUE)
    if (!all(elementLengths(encoding_blocks) == 2L))  # should never happen
        stop(encoding, ": invalid encoding")
    encoding_blocks <- matrix(unlist(encoding_blocks, use.names=FALSE), nrow=2L)
    Lencoding_patterns <- .build_encoding_patterns(ngap[1L])
    Lranks <- .extractSteppedExonRanksFromEncodingBlocks(encoding_blocks[1L, ],
                                                         Lencoding_patterns)
    Rencoding_patterns <- .build_encoding_patterns(ngap[2L])
    Rranks <- .extractSteppedExonRanksFromEncodingBlocks(encoding_blocks[2L, ],
                                                         Rencoding_patterns)
    if (length(Lranks) == 0L || length(Rranks) == 0L ||
        Lranks[length(Lranks)] > Rranks[1L]) {
        ranks <- integer(0)
        names(ranks) <- character(0)
        return(ranks)
    }
    if (for.query.right.end)
        return(Rranks)  # unnamed! (like for a single-end read)
    names(Rranks) <- rep.int("R", length(Rranks))
    names(Lranks) <- rep.int("L", length(Lranks))
    c(Lranks, Rranks)
}

setGeneric("extractSteppedExonRanks",
    function(x, for.query.right.end=FALSE)
        standardGeneric("extractSteppedExonRanks")
)

setMethod("extractSteppedExonRanks", "character",
    function(x, for.query.right.end=FALSE)
    {
        lapply(x, .extractSteppedExonRanks, for.query.right.end)
    }
)

setMethod("extractSteppedExonRanks", "factor",
    function(x, for.query.right.end=FALSE)
    {
        if (length(x) == 0L)
            return(list())
        ranks <- extractSteppedExonRanks(levels(x),
                     for.query.right.end=for.query.right.end)
        ranks[as.integer(x)]
    }
)

setMethod("extractSteppedExonRanks", "OverlapEncodings",
    function(x, for.query.right.end=FALSE)
    {
        ranks <- extractSteppedExonRanks(encoding(x),
                     for.query.right.end=for.query.right.end)
        ranks_elt_lens <- elementLengths(ranks)
        tmp <- unlist(unname(ranks), use.names=TRUE)  # we want the inner names
        tmp <- tmp + rep.int(Loffset(x), ranks_elt_lens)
        flevels <- seq_len(length(ranks))
        f <- factor(rep.int(flevels, ranks_elt_lens), levels=flevels)
        unname(split(tmp, f))
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### extractSpannedExonRanks().
###

setGeneric("extractSpannedExonRanks",
    function(x, for.query.right.end=FALSE)
        standardGeneric("extractSpannedExonRanks")
)

setMethod("extractSpannedExonRanks", "character",
    function(x, for.query.right.end=FALSE)
    {
        .extractRanks <- function(encoding) {
            ranks <- .extractSteppedExonRanks(encoding,
                          for.query.right.end=for.query.right.end)
            if (length(ranks) == 0L)
                return(c(NA_integer_, NA_integer_))
            c(ranks[1L], ranks[length(ranks)])
        }
        ranks <- lapply(x, .extractRanks)
        if (length(ranks) == 0L) {
            firstSpannedExonRank <- lastSpannedExonRank <- integer(0)
        } else {
            ranks <- unlist(ranks, use.names=FALSE)
            firstSpannedExonRank <- ranks[c(TRUE, FALSE)]
            lastSpannedExonRank <- ranks[c(FALSE, TRUE)]
        }
        data.frame(firstSpannedExonRank=firstSpannedExonRank,
                   lastSpannedExonRank=lastSpannedExonRank,
                   check.names=FALSE, stringsAsFactors=FALSE)
    }
)

setMethod("extractSpannedExonRanks", "factor",
    function(x, for.query.right.end=FALSE)
    {
        if (length(x) == 0L)
            return(list())
        ranks <- extractSpannedExonRanks(levels(x),
                     for.query.right.end=for.query.right.end)
        ans <- ranks[as.integer(x), , drop=FALSE]
        rownames(ans) <- NULL
        ans
    }
)

setMethod("extractSpannedExonRanks", "OverlapEncodings",
    function(x, for.query.right.end=FALSE)
    {
        ranks <- extractSpannedExonRanks(encoding(x),
                     for.query.right.end=for.query.right.end)
        ranks$firstSpannedExonRank <- ranks$firstSpannedExonRank + Loffset(x)
        ranks$lastSpannedExonRank <- ranks$lastSpannedExonRank + Loffset(x)
        ranks
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### extractSkippedExonRanks().
###

setGeneric("extractSkippedExonRanks",
    function(x, for.query.right.end=FALSE)
        standardGeneric("extractSkippedExonRanks")
)

setMethod("extractSkippedExonRanks", "character",
    function(x, for.query.right.end=FALSE)
    {
        .extractRanks <- function(encoding) {
            ranks <- .extractSteppedExonRanks(encoding,
                          for.query.right.end=for.query.right.end)
            if (length(ranks) == 0L)
                return(ranks)
            ranks_names <- names(ranks)
            if (is.null(ranks_names))  # single-end read
                return(setdiff(ranks[1L]:ranks[length(ranks)], ranks))
            ## Paired-end read.
            ranks <- split(unname(ranks), ranks_names)
            Lranks <- ranks$L
            Lranks <- setdiff(Lranks[1L]:Lranks[length(Lranks)], Lranks)
            Rranks <- ranks$R
            Rranks <- setdiff(Rranks[1L]:Rranks[length(Rranks)], Rranks)
            names(Lranks) <- rep.int("L", length(Lranks))
            names(Rranks) <- rep.int("R", length(Rranks))
            c(Lranks, Rranks)
        }
        lapply(x, .extractRanks)
    }
)

setMethod("extractSkippedExonRanks", "factor",
    function(x, for.query.right.end=FALSE)
    {
        if (length(x) == 0L)
            return(list())
        ranks <- extractSkippedExonRanks(levels(x),
                     for.query.right.end=for.query.right.end)
        ranks[as.integer(x)]
    }
)

setMethod("extractSkippedExonRanks", "OverlapEncodings",
    function(x, for.query.right.end=FALSE)
    {
        ranks <- extractSkippedExonRanks(encoding(x),
                     for.query.right.end=for.query.right.end)
        ranks_elt_lens <- elementLengths(ranks)
        tmp <- unlist(unname(ranks), use.names=TRUE)  # we want the inner names
        tmp <- tmp + rep.int(Loffset(x), ranks_elt_lens)
        flevels <- seq_len(length(ranks))
        f <- factor(rep.int(flevels, ranks_elt_lens), levels=flevels)
        unname(split(tmp, f))
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### extractQueryStartInTranscript().
###

### TODO: Maybe put this in IRanges and rename it setElementLengths, or even
### better introduce an "elementLengths<-" generic and make this the method
### for CompressedList objects?
.setElementLengths <- function(x, elt_lens)
{
    if (!is(x, "CompressedList"))
        stop("'x' must be a CompressedList object")
    if (!is.numeric(elt_lens) || length(elt_lens) != length(x))
        stop("'elt_lens' must be an integer vector of the same length as 'x'")
    if (!is.integer(elt_lens))
        elt_lens <- as.integer(elt_lens)
    if (IRanges:::anyMissingOrOutside(elt_lens, lower=0L))
        stop("'elt_lens' cannot contain NAs or negative values")
    x_elt_lens <- elementLengths(x)
    if (!all(elt_lens <= x_elt_lens))
        stop("'all(elt_lens <= elementLengths(x))' must be TRUE")
    offset <- cumsum(c(0L, x_elt_lens[-length(x_elt_lens)]))
    ii <- IRanges:::fancy_mseq(elt_lens, offset=offset)
    x@unlistData <- x@unlistData[ii]
    x@partitioning@end <- unname(cumsum(elt_lens))
    x
}

### Returns a data.frame with 1 row per overlap, and 3 integer columns:
###     1. startInTranscript
###     2. firstSpannedExonRank
###     3. startInFirstSpannedExon
### Rows for overlaps that are not "compatible" or "almost compatible"
### contain NAs.
extractQueryStartInTranscript <- function(query, subject,
                                          hits=NULL, ovenc=NULL,
                                          flip.query.if.wrong.strand=FALSE,
                                          for.query.right.end=FALSE)
{
    if (!is(query, "GRangesList") || !is(subject, "GRangesList"))
        stop("'query' and 'subject' must be GRangesList objects")
    seqinfo <- merge(seqinfo(query), seqinfo(subject))
    seqlevels(query) <- seqlevels(subject) <- seqlevels(seqinfo)
    if (is.null(hits)) {
        if (length(query) != length(subject))
            stop("'query' and 'subject' must have the same length")
    } else {
        if (!is(hits, "Hits"))
            stop("'hits' must be a Hits object or NULL")
        if (queryLength(hits) != length(query) ||
            subjectLength(hits) != length(subject))
            stop("'hits' is not compatible with 'query' and 'subject' ",
                 "('queryLength(hits)' and 'subjectLength(hits)' don't ",
                 "match the lengths of 'query' and 'subject')")
        query <- query[queryHits(hits)]
        subject <- subject[subjectHits(hits)]
    }
    if (is.null(ovenc)) {
        ovenc <- encodeOverlaps(query, subject,
                         flip.query.if.wrong.strand=flip.query.if.wrong.strand)
    } else {
        if (!is(ovenc, "OverlapEncodings"))
            stop("'ovenc' must be an OverlapEncodings object")
        if (length(ovenc) != length(query))
            stop("when not NULL, 'ovenc' must have the same length ",
                 "as 'hits', if specified, otherwiseaas 'query'")
    }
    if (!isTRUEorFALSE(for.query.right.end))
        stop("'for.query.right.end' must be TRUE or FALSE")

    query <- flipQuery(query, flippedQuery(ovenc))

    ## Extract start/end/strand of the first range
    ## in each top-level element of 'query'.
    qii1 <- start(query@partitioning)
    if (for.query.right.end) {
        query.break <- elementMetadata(query)$query.break
        if (is.null(query.break))
            stop("using 'for.query.right.end=TRUE' requires that ",
                 "'elementMetadata(query)' has a \"query.break\" column ",
                 "indicating for each paired-end read the position of the ",
                 "break between the ranges coming from one end and those ",
                 "coming from the other end")
        qii1 <- qii1 + query.break
    }
    query_start1 <- start(query@unlistData)[qii1]
    query_end1 <- end(query@unlistData)[qii1]
    query_strand1 <- as.factor(strand(query@unlistData))[qii1]

    ## Extract start/end/strand of the first spanned exon
    ## in each top-level element of 'subject'.
    exrank <- extractSpannedExonRanks(ovenc,
                  for.query.right.end=for.query.right.end)$firstSpannedExonRank
    sii1 <- start(subject@partitioning) + exrank - 1L
    subject_start1 <- start(subject@unlistData)[sii1]
    subject_end1 <- end(subject@unlistData)[sii1]
    subject_strand1 <- as.factor(strand(subject@unlistData))[sii1]

    ## A sanity check.
    if (any(!is.na(exrank) & (query_strand1 != subject_strand1))) {
        ## TODO: Error message needs to take into account whether 'hits'
        ## and/or 'ovenc' was supplied or not.
        stop("'ovenc' is incompatible with the supplied 'query' ",
             "and/or 'subject' and/or 'hits'")
    }

    ## Compute the "query start in first spanned exon".
    startInFirstSpannedExon <- rep.int(NA_integer_, length(query))
    is_on_plus <- query_strand1 == "+"
    idx <- which(!is.na(exrank) & is_on_plus)
    startInFirstSpannedExon[idx] <- query_start1[idx] - subject_start1[idx] + 1L
    idx <- which(!is.na(exrank) & !is_on_plus)
    startInFirstSpannedExon[idx] <- subject_end1[idx] - query_end1[idx] + 1L

    ## Truncate each transcript in 'subject' right before the first spanned
    ## exon and compute the cumulated width of the truncated object.
    subject2_elt_lens <- exrank - 1L
    subject2_elt_lens[is.na(exrank)] <- 0L
    subject2 <- .setElementLengths(subject, subject2_elt_lens)
    subject2_cumwidth <- unname(sum(width(subject2)))
    subject2_cumwidth[is.na(exrank)] <- NA_integer_

    ## Compute the "query start in transcript".
    startInTranscript <- subject2_cumwidth + startInFirstSpannedExon

    data.frame(startInTranscript=startInTranscript,
               firstSpannedExonRank=exrank,
               startInFirstSpannedExon=startInFirstSpannedExon,
               check.names=FALSE, stringsAsFactors=FALSE)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### High-level convenience wrappers: findCompatibleOverlaps() and
### countCompatibleOverlaps().
###

setGeneric("findCompatibleOverlaps",
    function(query, subject) standardGeneric("findCompatibleOverlaps")
)

.GappedAlignmentsORGappedAlignmentPairs.findCompatibleOverlaps <-
    function(query, subject)
{
    grl <- grglist(query, order.as.in.query=TRUE)
    ## TODO: Use 'type="within"' when it's supported for circular
    ## sequences like the mitochondrial chromosome.
    ov <- findOverlaps(grl, subject, ignore.strand=TRUE)
    ovenc <- encodeOverlaps(grl, subject, hits=ov,
                            flip.query.if.wrong.strand=TRUE)
    ov_is_compat <- isCompatibleWithSplicing(ovenc)
    ov[ov_is_compat]
}

setMethod("findCompatibleOverlaps", c("GappedAlignments", "GRangesList"),
    .GappedAlignmentsORGappedAlignmentPairs.findCompatibleOverlaps
)

setMethod("findCompatibleOverlaps", c("GappedAlignmentPairs", "GRangesList"),
    .GappedAlignmentsORGappedAlignmentPairs.findCompatibleOverlaps
)

countCompatibleOverlaps <- function(query, subject)
{
    compatov <- findCompatibleOverlaps(query, subject)
    tabulate(queryHits(compatov), nbins=queryLength(compatov))
}

