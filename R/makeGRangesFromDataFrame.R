### =========================================================================
### makeGRangesFromDataFrame()
### -------------------------------------------------------------------------


.normarg_field <- function(field, what)
{
    if (!is.character(field) || any(is.na(field)))
        stop("'", what, ".field' must be a character vector with no NAs")
    tolower(field)
}

.collect_prefixes <- function(df_colnames, field)
{
    df_colnames_nc <- nchar(df_colnames)
    prefixes <- lapply(field,
        function(suf) {
            pref_nc <- df_colnames_nc - nchar(suf)
            idx <- which(substr(df_colnames, pref_nc + 1L, df_colnames_nc) ==
                         suf)
            substr(df_colnames[idx], 1L, pref_nc[idx])
        })
    unique(unlist(prefixes))
}

.find_start_end_cols <- function(df_colnames, start.field, end.field)
{
    idx1 <- which(df_colnames %in% start.field)
    idx2 <- which(df_colnames %in% end.field)
    if (length(idx1) == 1L && length(idx2) == 1L)
        return(list(c(start=idx1, end=idx2), ""))
    if (length(idx1) == 0L && length(idx2) == 0L) {
        prefixes1 <- .collect_prefixes(df_colnames, start.field)
        prefixes2 <- .collect_prefixes(df_colnames, end.field)
        if (length(prefixes1) == 1L && length(prefixes2) == 1L
         && prefixes1 == prefixes2)
        {
            prefix <- prefixes1
            idx1 <- which(df_colnames %in% paste0(prefix, start.field))
            idx2 <- which(df_colnames %in% paste0(prefix, end.field))
            if (length(idx1) == 1L && length(idx2) == 1L)
                return(list(c(start=idx1, end=idx2), prefix))
        }
    }
    stop("cannnot determine start/end columns")
}

.find_width_col <- function(df_colnames, width.field, prefix)
{
    idx <- which(df_colnames %in% paste0(prefix, width.field))
    if (length(idx) == 0L)
        idx <- which(df_colnames %in% width.field)
    if (length(idx) == 0L)
        return(NA_integer_)
    if (length(idx) >= 2L) {
        warning("cannnot determine width column unambiguously")
        return(idx[[1L]])
    }
    idx
}

.find_seqnames_col <- function(df_colnames, seqnames.field, prefix)
{
    idx <- which(df_colnames %in% paste0(prefix, seqnames.field))
    if (length(idx) == 0L)
        idx <- which(df_colnames %in% seqnames.field)
    if (length(idx) == 0L)
        stop("cannnot find seqnames column")
    if (length(idx) >= 2L)
        stop("cannnot determine seqnames column unambiguously")
    idx
}

.find_strand_col <- function(df_colnames, strand.field, prefix)
{
    idx <- which(df_colnames %in% paste0(prefix, strand.field))
    if (length(idx) == 0L) 
        idx <- which(df_colnames %in% strand.field)
    if (length(idx) == 0L)
        return(NA_integer_)
    if (length(idx) >= 2L)
        stop("Cannnot determine strand column unambiguously. ",
             "(You can use\n  'ignore.strand=TRUE' to ignore ",
             "strand information.)")
    idx
}

### Returns a named integer vector of length 5. Names are: seqnames, start,
### end, width, and strand. The values must be valid column numbers, except
### for the width and strand elements that can also be NAs.
.find_GRanges_cols <- function(df_colnames,
                               seqnames.field=c("seqnames", "seqname",
                                                "chromosome", "chrom",
                                                "chr", "chromosome_name",
                                                "seqid"),
                               start.field="start",
                               end.field=c("end", "stop"),
                               strand.field="strand",
                               ignore.strand=FALSE)
{
    ## Automatic detection of seqnames/start/end/strand columns is case
    ## insensitive.
    df_colnames0 <- tolower(df_colnames)
    seqnames.field0 <- .normarg_field(seqnames.field, "seqnames")
    start.field0 <- .normarg_field(start.field, "start")
    end.field0 <- .normarg_field(end.field, "end")

    start_end_cols <- .find_start_end_cols(df_colnames0,
                                           start.field0,
                                           end.field0)
    prefix <- start_end_cols[[2L]]
    ## Name of "width" field is not under user control for now (until we need
    ## need that).
    width_col <- .find_width_col(df_colnames0, "width", prefix)
    seqnames_col <- .find_seqnames_col(df_colnames0,
                                       seqnames.field0,
                                       prefix)
    if (ignore.strand) {
        strand_col <- NA_integer_
    } else {
        strand.field0 <- .normarg_field(strand.field, "strand")
        strand_col <- .find_strand_col(df_colnames0,
                                       strand.field0,
                                       prefix)
    }
    c(seqnames=seqnames_col, start_end_cols[[1L]], width=width_col,
      strand=strand_col)
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
                                     starts.in.df.are.0based=FALSE)
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

    granges_cols <- .find_GRanges_cols(names(df),
                                       seqnames.field=seqnames.field,
                                       start.field=start.field,
                                       end.field=end.field,
                                       strand.field=strand.field,
                                       ignore.strand=ignore.strand)

    ## Prepare 'ans_seqnames'.
    ans_seqnames <- df[[granges_cols[["seqnames"]]]]

    ## Prepare 'ans_ranges'.
    ans_start <- .get_data_frame_col_as_numeric(df, granges_cols[["start"]])
    ans_end <- .get_data_frame_col_as_numeric(df, granges_cols[["end"]])
    if (starts.in.df.are.0based)
        ans_start <- ans_start + 1L
    ans_names <- rownames(df)
    if (identical(ans_names, as.character(seq_len(nrow(df)))))
        ans_names <- NULL
    ans_ranges <- IRanges(ans_start, ans_end, names=ans_names)

    ## Prepare 'ans_strand'.
    if (is.na(granges_cols[["strand"]]) || ignore.strand) {
        ans_strand <- "*"
    } else {
        ans_strand <- as.character(df[[granges_cols[["strand"]]]])
        ans_strand[ans_strand %in% "."] <- "*"
    }

    ## Prepare 'ans_mcols'.
    if (keep.extra.columns) {
        drop_idx <- c(granges_cols[["seqnames"]],
                      granges_cols[["start"]],
                      granges_cols[["end"]])
        if (!is.na(granges_cols[["width"]]))
            drop_idx <- c(drop_idx, granges_cols[["width"]])
        if (!is.na(granges_cols[["strand"]]))
            drop_idx <- c(drop_idx, granges_cols[["strand"]])
        ans_mcols <- df[-drop_idx]
    } else {
        ans_mcols <- NULL
    }

    ## Prepare 'ans_seqinfo'.
    if (is.null(ans_seqinfo)) {
        ## Only if 'ans_seqnames' is a factor-Rle, we preserve the seqlevels
        ## in the order they are in 'levels(ans_seqnames)'. Otherwise, we
        ## order them according to rankSeqlevels().
        seqlevels <- levels(ans_seqnames)
        if (is.null(seqlevels)) {
            seqlevels <- unique(ans_seqnames)
            if (!is.character(seqlevels))
                seqlevels <- as.character(seqlevels)
        }
        if (!(is(ans_seqnames, "Rle") && is.factor(runValue(ans_seqnames))))
            seqlevels[rankSeqlevels(seqlevels)] <- seqlevels
        ans_seqinfo <- Seqinfo(seqlevels)
    }

    ## Make and return the GRanges object.
    GRanges(ans_seqnames, ans_ranges, strand=ans_strand,
            ans_mcols, seqinfo=ans_seqinfo)
}

setAs("data.frame", "GRanges",
    function(from) makeGRangesFromDataFrame(from, keep.extra.columns=TRUE)
)

setAs("DataFrame", "GRanges",
    function(from) makeGRangesFromDataFrame(from, keep.extra.columns=TRUE)
)

