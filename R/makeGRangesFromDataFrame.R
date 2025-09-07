### =========================================================================
### makeGRangesFromDataFrame() and makeGPosFromDataFrame()
### -------------------------------------------------------------------------


### Vectorized.
.has_suffix <- function(x, suffix)
{
    stopifnot(is.character(x), isSingleString(suffix))
    x_nc <- nchar(x)
    substr(x, x_nc - nchar(suffix) + 1L, x_nc) == suffix
}

.collect_prefixes <- function(x, suffixes)
{
    stopifnot(is.character(x), is.character(suffixes))
    all_prefixes <- lapply(suffixes,
        function(suffix) {
            ok <- .has_suffix(x, suffix) & x != suffix
            x2 <- x[ok]
	    prefix_nc <- nchar(x2) - nchar(suffix)
            substr(x2, 1L, prefix_nc)
        })
    unique(unlist(all_prefixes, use.names=FALSE))
}

.normarg_field <- function(field, what)
{
    if (!is.character(field) || any(is.na(field)))
        stop("'", what, ".field' must be a character vector with no NAs")
    tolower(field)
}

.find_seqnames_col <- function(df_colnames, seqnames.field, startend_prefix="")
{
    colidx <- which(df_colnames %in% paste0(startend_prefix, seqnames.field))
    if (length(colidx) == 0L)
        colidx <- which(df_colnames %in% seqnames.field)
    if (length(colidx) == 0L)
        stop("cannnot find seqnames column")
    if (length(colidx) >= 2L)
        stop("cannnot determine seqnames column unambiguously")
    colidx
}

.find_strand_col <- function(df_colnames, strand.field, ignore.strand,
                             startend_prefix="")
{
    if (ignore.strand)
        return(NA_integer_)
    colidx <- which(df_colnames %in% paste0(startend_prefix, strand.field))
    if (length(colidx) == 0L)
        colidx <- which(df_colnames %in% strand.field)
    if (length(colidx) == 0L)
        return(NA_integer_)
    if (length(colidx) >= 2L)
        stop("Cannnot determine strand column unambiguously. ",
             "(You can use\n  'ignore.strand=TRUE' to ignore ",
             "strand information.)")
    colidx
}

.get_data_frame_col_as_numeric <- function(df, col)
{
    ans <- df[[col]]
    if (is(ans, "Rle"))
        ans <- S4Vectors:::decodeRle(ans)
    if (!is.numeric(ans)) {
        if (is.factor(ans))
            ans <- as.character(ans)
        ans <- suppressWarnings(as.numeric(ans))
        if (anyNA(ans))
            stop(wmsg("some values in the ",
                      "\"", names(df)[[col]], "\" ",
                      "column cannot be turned into numeric values"))
    }
    ans
}

.get_strand_from_data_frame <- function(df, corecol_map, ignore.strand)
{
    colidx <- corecol_map[["strand"]]
    if (is.na(colidx) || ignore.strand)
        return(NULL)
    strand <- as.character(df[[colidx]])
    strand[strand %in% "."] <- "*"
    strand
}

.make_dummy_Seqinfo_from_seqnames <- function(seqnames)
{
    ## If 'seqnames' has levels (i.e. is a factor or factor-Rle) then
    ## we use them as-is, that is, in their original order. Otherwise,
    ## we infer them from 'unique(seqnames)' and order them according
    ## to rankSeqlevels().
    seqlevels <- levels(seqnames)
    if (is.null(seqlevels)) {
        seqlevels <- unique(seqnames)
        if (!is.character(seqlevels))
            seqlevels <- as.character(seqlevels)
        seqlevels[rankSeqlevels(seqlevels)] <- seqlevels
    }
    Seqinfo(seqlevels)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .find_core_GRanges_cols()
###

.find_start_end_cols <- function(df_colnames, start.field, end.field,
                                 startend_prefix="")
{
    colidx1 <- which(df_colnames %in% paste0(startend_prefix, start.field))
    colidx2 <- which(df_colnames %in% paste0(startend_prefix, end.field))
    if (length(colidx1) == 1L && length(colidx2) == 1L)
        return(list(c(start=colidx1, end=colidx2), startend_prefix))
    if (length(colidx1) != 0L || length(colidx2) != 0L ||
        startend_prefix != "")
    {
        stop(wmsg("cannnot determine start/end columns"))
    }
    prefixes1 <- .collect_prefixes(df_colnames, start.field)
    prefixes2 <- .collect_prefixes(df_colnames, end.field)
    if (length(prefixes1) != 1L || length(prefixes2) != 1L ||
        prefixes1 != prefixes2)
    {
        stop(wmsg("cannnot determine start/end columns"))
    }
    startend_prefix <- prefixes1
    .find_start_end_cols(df_colnames, start.field, end.field, startend_prefix)
}

.find_width_col <- function(df_colnames, width.field, startend_prefix)
{
    colidx <- which(df_colnames %in% paste0(startend_prefix, width.field))
    if (length(colidx) == 0L)
        colidx <- which(df_colnames %in% width.field)
    if (length(colidx) == 0L)
        return(NA_integer_)
    if (length(colidx) >= 2L) {
        warning("cannnot determine width column unambiguously")
        return(colidx[[1L]])
    }
    colidx
}

### The 5 core GRanges columns are: seqnames, start, end, width, strand.
### Returns a named integer vector where the names are the 5 core GRanges
### columns and the values are valid column indices. The "width"
### and/or "strand" column indices can be NAs.
.find_core_GRanges_cols <-
    function(df_colnames,
             seqnames.field=c("seqnames", "seqname",
                              "chromosome", "chrom",
                              "chr", "chromosome_name",
                              "seqid"),
             start.field="start",
             end.field=c("end", "stop"),
             strand.field="strand",
             ignore.strand=FALSE)
{
    ## The heuristic we use to find the core GRanges columns is
    ## case insensitive.
    df_colnames0 <- tolower(df_colnames)
    seqnames.field0 <- .normarg_field(seqnames.field, "seqnames")
    start.field0 <- .normarg_field(start.field, "start")
    end.field0 <- .normarg_field(end.field, "end")
    strand.field0 <- .normarg_field(strand.field, "strand")

    start_end_cols <- .find_start_end_cols(df_colnames0,
                                           start.field0,
                                           end.field0)
    startend_prefix <- start_end_cols[[2L]]
    ## Name of "width" field is not under user control for now (until we
    ## need that).
    width_col <- .find_width_col(df_colnames0, "width", startend_prefix)
    seqnames_col <- .find_seqnames_col(df_colnames0, seqnames.field0,
                                       startend_prefix)
    strand_col <- .find_strand_col(df_colnames0, strand.field0, ignore.strand,
                                   startend_prefix)

    c(seqnames=seqnames_col, start_end_cols[[1L]], width=width_col,
      strand=strand_col)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### makeGRangesFromDataFrame()
###

.drop_rows_with_na_start_end <- function(df, corecol_map, na.rm)
{
    start_col <- .get_data_frame_col_as_numeric(df, corecol_map[["start"]])
    end_col <- .get_data_frame_col_as_numeric(df, corecol_map[["end"]])
    is_na <- is.na(start_col) | is.na(end_col)
    if (!any(is_na))
        return(df)
    if (na.rm) {
        keep_idx <- which(!is_na)
        return(S4Vectors:::extract_data_frame_rows(df, keep_idx))
    }
    start_colname <- names(df)[[corecol_map[["start"]]]]
    end_colname <- names(df)[[corecol_map[["end"]]]]
    stop(wmsg(
        "The \"", start_colname, "\" and/or \"", end_colname, "\" ",
        "columns contain NAs. Use 'na.rm=TRUE' to ignore the rows ",
        "with NAs."
    ))
}

### 'df' must be a data.frame or DataFrame object.
makeGRangesFromDataFrame <- function(df,
                                     keep.extra.columns=FALSE,
                                     ignore.strand=FALSE,
                                     seqinfo=NULL,
                                     seqnames.field=c("seqnames", "seqname",
                                                      "chromosome", "chrom",
                                                      "chr", "chromosome_name",
                                                      "seqid"),
                                     start.field="start",
                                     end.field=c("end", "stop"),
                                     strand.field="strand",
                                     starts.in.df.are.0based=FALSE,
                                     na.rm=FALSE)
{
    ## Check args.
    if (is.character(df))  # for people that provide the path to a file
        stop("'df' must be a data.frame or DataFrame object")
    if (!(is.data.frame(df) || is(df, "DataFrame")))
        df <- as.data.frame(df)
    if (!isTRUEorFALSE(keep.extra.columns))
        stop("'keep.extra.columns' must be TRUE or FALSE")
    if (!isTRUEorFALSE(ignore.strand))
        stop("'ignore.strand' must be TRUE or FALSE")
    ans_seqinfo <- normarg_seqinfo1(seqinfo)
    if (!isTRUEorFALSE(starts.in.df.are.0based))
        stop("'starts.in.df.are.0based' must be TRUE or FALSE")
    if (!isTRUEorFALSE(na.rm))
        stop(wmsg("'na.rm' must be TRUE or FALSE"))

    corecol_map <- .find_core_GRanges_cols(names(df),
                                           seqnames.field=seqnames.field,
                                           start.field=start.field,
                                           end.field=end.field,
                                           strand.field=strand.field,
                                           ignore.strand=ignore.strand)
    df <- .drop_rows_with_na_start_end(df, corecol_map, na.rm)

    ## Prepare 'ans_ranges'.
    ans_start <- .get_data_frame_col_as_numeric(df, corecol_map[["start"]])
    ans_end <- .get_data_frame_col_as_numeric(df, corecol_map[["end"]])
    if (starts.in.df.are.0based)
        ans_start <- ans_start + 1L
    ans_names <- rownames(df)
    if (identical(ans_names, as.character(seq_len(nrow(df)))))
        ans_names <- NULL
    ans_ranges <- IRanges(ans_start, ans_end, names=ans_names)

    ## Prepare 'ans_seqnames'.
    ans_seqnames <- df[[corecol_map[["seqnames"]]]]

    ## Prepare 'ans_strand'.
    ans_strand <- .get_strand_from_data_frame(df, corecol_map, ignore.strand)

    ## Prepare 'ans_mcols'.
    if (keep.extra.columns) {
        drop_colidx <- c(corecol_map[["seqnames"]],
                         corecol_map[["start"]],
                         corecol_map[["end"]])
        if (!is.na(corecol_map[["width"]]))
            drop_colidx <- c(drop_colidx, corecol_map[["width"]])
        if (!is.na(corecol_map[["strand"]]))
            drop_colidx <- c(drop_colidx, corecol_map[["strand"]])
        ans_mcols <- df[-drop_colidx]
    } else {
        ans_mcols <- NULL
    }

    ## Prepare 'ans_seqinfo'.
    if (is.null(ans_seqinfo))
        ans_seqinfo <- .make_dummy_Seqinfo_from_seqnames(ans_seqnames)

    ## Make and return the GRanges object.
    GRanges(ans_seqnames, ans_ranges, strand=ans_strand, ans_mcols,
            seqinfo=ans_seqinfo)
}

setAs("data.frame", "GRanges",
    function(from) makeGRangesFromDataFrame(from, keep.extra.columns=TRUE)
)

setAs("DataFrame", "GRanges",
    function(from) makeGRangesFromDataFrame(from, keep.extra.columns=TRUE)
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .find_core_GPos_cols()
###

.find_pos_col <- function(df_colnames, pos.field)
{
    colidx <- which(df_colnames %in% pos.field)
    if (length(colidx) == 0L)
        return(NA_integer_)
    if (length(colidx) >= 2L) {
        warning("cannnot determine pos column unambiguously")
        return(colidx[[1L]])
    }
    colidx
}

### The 3 core GPos columns are: seqnames, pos, strand.
### Returns a named integer vector where the names are the 3 core GPos
### columns and the values are valid column indices. The "strand" column
### index can be NA.
.find_core_GPos_cols <-
    function(df_colnames,
             seqnames.field=c("seqnames", "seqname",
                              "chromosome", "chrom",
                              "chr", "chromosome_name",
                              "seqid"),
             pos.field=c("pos", "position", "positions"),
             strand.field="strand",
             ignore.strand=FALSE)
{
    ## The heuristic we use to find the core GPos columns is case insensitive.
    df_colnames0 <- tolower(df_colnames)
    seqnames.field0 <- .normarg_field(seqnames.field, "seqnames")
    pos.field0 <- .normarg_field(pos.field, "pos")
    strand.field0 <- .normarg_field(strand.field, "strand")

    seqnames_col <- .find_seqnames_col(df_colnames0, seqnames.field0)
    pos_col <- .find_pos_col(df_colnames0, pos.field0)
    strand_col <- .find_strand_col(df_colnames0, strand.field0, ignore.strand)

    c(seqnames=seqnames_col, pos=pos_col, strand=strand_col)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### makeGPosFromDataFrame()
###

.drop_rows_with_na_pos <- function(df, corecol_map, na.rm)
{
    pos_col <- .get_data_frame_col_as_numeric(df, corecol_map[["pos"]])
    is_na <- is.na(pos_col)
    if (!any(is_na))
        return(df)
    if (na.rm) {
        keep_idx <- which(!is_na)
        return(S4Vectors:::extract_data_frame_rows(df, keep_idx))
    }
    pos_colname <- names(df)[[corecol_map[["pos"]]]]
    stop(wmsg(
        "The \"", pos_colname, "\" column contains NAs. ",
        "Use 'na.rm=TRUE' to ignore the rows with NAs."
    ))
}

### 'df' must be a data.frame or DataFrame object.
makeGPosFromDataFrame <- function(df,
                                  keep.extra.columns=FALSE,
                                  ignore.strand=FALSE,
                                  seqinfo=NULL,
                                  seqnames.field=c("seqnames", "seqname",
                                                   "chromosome", "chrom",
                                                   "chr", "chromosome_name",
                                                   "seqid"),
                                  pos.field=c("pos", "position", "positions"),
                                  strand.field="strand",
                                  na.rm=FALSE)
{
    ## Check args.
    if (is.character(df))  # for people that provide the path to a file
        stop("'df' must be a data.frame or DataFrame object")
    if (!(is.data.frame(df) || is(df, "DataFrame")))
        df <- as.data.frame(df)
    if (!isTRUEorFALSE(keep.extra.columns))
        stop("'keep.extra.columns' must be TRUE or FALSE")
    if (!isTRUEorFALSE(ignore.strand))
        stop("'ignore.strand' must be TRUE or FALSE")
    ans_seqinfo <- normarg_seqinfo1(seqinfo)
    if (!isTRUEorFALSE(na.rm))
        stop(wmsg("'na.rm' must be TRUE or FALSE"))

    corecol_map <- .find_core_GPos_cols(names(df),
                                        seqnames.field=seqnames.field,
                                        pos.field=pos.field,
                                        strand.field=strand.field,
                                        ignore.strand=ignore.strand)
    df <- .drop_rows_with_na_pos(df, corecol_map, na.rm)

    ## Prepare 'ans_pos'.
    ans_pos <- .get_data_frame_col_as_numeric(df, corecol_map[["pos"]])
    ans_names <- rownames(df)
    if (identical(ans_names, as.character(seq_len(nrow(df)))))
        ans_names <- NULL
    ans_pos <- IPos(ans_pos, names=ans_names)

    ## Prepare 'ans_seqnames'.
    ans_seqnames <- df[[corecol_map[["seqnames"]]]]

    ## Prepare 'ans_strand'.
    ans_strand <- .get_strand_from_data_frame(df, corecol_map, ignore.strand)

    ## Prepare 'ans_mcols'.
    if (keep.extra.columns) {
        drop_colidx <- c(corecol_map[["seqnames"]],
                         corecol_map[["pos"]])
        if (!is.na(corecol_map[["strand"]]))
            drop_colidx <- c(drop_colidx, corecol_map[["strand"]])
        ans_mcols <- df[-drop_colidx]
    } else {
        ans_mcols <- NULL
    }

    ## Prepare 'ans_seqinfo'.
    if (is.null(ans_seqinfo))
        ans_seqinfo <- .make_dummy_Seqinfo_from_seqnames(ans_seqnames)

    ## Make and return the GPos object.
    GPos(ans_seqnames, ans_pos, strand=ans_strand, ans_mcols,
         seqinfo=ans_seqinfo)
}

setAs("data.frame", "GPos",
    function(from) makeGPosFromDataFrame(from, keep.extra.columns=TRUE)
)

setAs("DataFrame", "GPos",
    function(from) makeGPosFromDataFrame(from, keep.extra.columns=TRUE)
)

