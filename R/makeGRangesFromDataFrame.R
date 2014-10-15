### =========================================================================
### makeGRangesFromDataFrame()
### -------------------------------------------------------------------------


### Must return NULL or a Seqinfo object.
.normarg_seqinfo <- function(seqinfo)
{
    if (is.null(seqinfo) || is(seqinfo, "Seqinfo"))
        return(seqinfo)
    if (is.character(seqinfo))
        return(Seqinfo(seqinfo))
    if (is.numeric(seqinfo)) {
        seqlevels <- names(seqinfo)
        if (is.null(seqlevels))
            stop("when a numeric vector, 'seqinfo' must have names")
        return(Seqinfo(seqlevels, seqlengths=seqinfo))
    }
    stop("'seqinfo' must be NULL, or a Seqinfo object, or a character vector ",
         "of seqlevels, or a named numeric vector of sequence lengths")
}

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
             "(You can use\n  'ignore.strand=FALSE' to ignore ",
             "strand information.)")
    idx
}

### Returns an integer vector of length 4 with names "seqnames", "start",
### "end", and "strand".
.find_GRanges_cols <- function(df_colnames,
                               seqnames.field=c("seqnames", "seqname",
                                                "chromosome", "chrom",
                                                "chr", "chromosome_name"),
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
    c(seqnames=seqnames_col, start_end_cols[[1L]], strand=strand_col)
}

### 'df' must be a data.frame or DataFrame object.
makeGRangesFromDataFrame <- function(df,
                                     keep.extra.columns=FALSE,
                                     ignore.strand=FALSE,
                                     seqinfo=NULL,
                                     seqnames.field=c("seqnames", "seqname",
                                                      "chromosome", "chrom",
                                                      "chr", "chromosome_name"),
                                     start.field="start",
                                     end.field=c("end", "stop"),
                                     strand.field="strand",
                                     starts.in.df.are.0based=FALSE)
{
    ## Check args.
    if (!is.data.frame(df) && !is(df, "DataFrame"))
        stop("'df' must be a data.frame or DataFrame object")
    if (!isTRUEorFALSE(keep.extra.columns))
        stop("'keep.extra.columns' must be TRUE or FALSE")
    if (!isTRUEorFALSE(ignore.strand))
        stop("'ignore.strand' must be TRUE or FALSE")
    ans_seqinfo <- .normarg_seqinfo(seqinfo)
    if (!isTRUEorFALSE(starts.in.df.are.0based))
        stop("'starts.in.df.are.0based' must be TRUE or FALSE")

    granges_cols <- .find_GRanges_cols(names(df),
                                       seqnames.field=seqnames.field,
                                       start.field=start.field,
                                       end.field=end.field,
                                       strand.field=strand.field,
                                       ignore.strand=ignore.strand)

    ## Prepare the GRanges components.
    ans_seqnames <- df[[granges_cols[["seqnames"]]]]
    ans_start <- df[[granges_cols[["start"]]]]
    ans_end <- df[[granges_cols[["end"]]]]
    if (!is.numeric(ans_start) || !is.numeric(ans_end))
        stop("\"", names(df)[granges_cols[["start"]]], "\" and ",
             "\"", names(df)[granges_cols[["end"]]], "\" columns ",
             "must be numeric")
    if (starts.in.df.are.0based)
        ans_start <- ans_start + 1L
    ans_ranges <- IRanges(ans_start, ans_end)
    if (is.na(granges_cols[["strand"]]) || ignore.strand) {
        ans_strand <- "*"
    } else {
        ans_strand <- df[[granges_cols[["strand"]]]]
    }
    if (keep.extra.columns) {
        drop_idx <- c(granges_cols[["seqnames"]],
                      granges_cols[["start"]],
                      granges_cols[["end"]])
        if (!is.na(granges_cols[["strand"]]))
            drop_idx <- c(drop_idx, granges_cols[["strand"]])
        ans_mcols <- df[-drop_idx]
    } else {
        ans_mcols <- NULL
    }
    ans_names <- rownames(df)
    if (identical(as.character(seq_len(nrow(df))), ans_names))
        ans_names <- NULL

    ## Make the GRanges object and return it.
    ans <- GRanges(ans_seqnames, ans_ranges, strand=ans_strand,
                   ans_mcols, seqinfo=ans_seqinfo)
    if (!is.null(ans_names))
        names(ans) <- ans_names
    ans
}

setAs("data.frame", "GRanges",
    function(from) makeGRangesFromDataFrame(from, keep.extra.columns=TRUE)
)

setAs("DataFrame", "GRanges",
    function(from) makeGRangesFromDataFrame(from, keep.extra.columns=TRUE)
)

