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

.get_field_pos <- function(field, df, what, required=TRUE)
{
    if (!is.character(field) || any(is.na(field)))
        stop("'", what, ".field' must be a character vector with no NAs")
    pos <- match(field, names(df))
    pos <- pos[which(!is.na(pos))[1L]]
    if (required && is.na(pos))
        stop("no field listed in '", what, ".field' is present in 'df'")
    pos
}

### 'df' must be a data.frame or DataFrame object.
makeGRangesFromDataFrame <- function(df,
    keep.extra.columns=FALSE,
    ignore.strand=FALSE,
    seqinfo=NULL,
    seqnames.field=c("seqnames", "chr", "chrom"),
    start.field=c("start", "chromStart"),
    end.field=c("end", "chromEnd", "stop", "chromStop"),
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
    seqnames_fpos <- .get_field_pos(seqnames.field, df, "seqnames")
    start_fpos <- .get_field_pos(start.field, df, "start")
    end_fpos <- .get_field_pos(end.field, df, "end")
    strand_fpos <- .get_field_pos(strand.field, df, "strand", required=FALSE)
    if (!isTRUEorFALSE(starts.in.df.are.0based))
        stop("'starts.in.df.are.0based' must be TRUE or FALSE")

    ## Prepare the GRanges components.
    ans_seqnames <- df[[seqnames_fpos]]
    ans_start <- df[[start_fpos]]
    ans_end <- df[[end_fpos]]
    if (!is.numeric(ans_start) || !is.numeric(ans_end))
        stop("\"", names(df)[start_fpos], "\" and ",
             "\"", names(df)[end_fpos], "\" columns must be numeric")
    if (starts.in.df.are.0based)
        ans_start <- ans_start + 1L
    ans_ranges <- IRanges(ans_start, ans_end)
    if (is.na(strand_fpos) || ignore.strand) {
        ans_strand <- "*"
    } else {
        ans_strand <- df[[strand_fpos]]
    }
    if (keep.extra.columns) {
        drop_idx <- c(seqnames_fpos, start_fpos, end_fpos)
        if (!is.na(strand_fpos))
            drop_idx <- c(drop_idx, strand_fpos)
        ans_mcols <- df[-drop_idx]
    } else {
        ans_mcols <- NULL
    }
    ans_names <- rownames(df)
    if (identical(as.character(seq_len(nrow(df))), ans_names))
        ans_names <- NULL

    ## Make the GRanges object and return it.
    ans <- GRanges(seqnames=ans_seqnames, ranges=ans_ranges, strand=ans_strand,
                   ans_mcols, seqinfo=ans_seqinfo)
    if (!is.null(ans_names))
        names(ans) <- ans_names
    ans
}

