###

test_find_GRanges_cols <- function()
{
    find_GRanges_cols <- GenomicRanges:::.find_GRanges_cols

    df_colnames <- c("chrom", "start", "end")
    target <- c(seqnames=1L, start=2L, end=3L,
                width=NA_integer_, strand=NA_integer_)
    current <- find_GRanges_cols(df_colnames)
    checkIdentical(target, current)

    df_colnames <- c("stop", "start", "Chr")
    target <- c(seqnames=3L, start=2L, end=1L,
                width=NA_integer_, strand=NA_integer_)
    current <- find_GRanges_cols(df_colnames)
    checkIdentical(target, current)

    df_colnames <- c("stop", "width", "start", "Chr")
    target <- c(seqnames=4L, start=3L, end=1L,
                width=2L, strand=NA_integer_)
    current <- find_GRanges_cols(df_colnames)
    checkIdentical(target, current)

    df_colnames <- c("strand", "STOP", "START", "chromosome_name")
    target <- c(seqnames=4L, start=3L, end=2L,
                width=NA_integer_, strand=1L)
    current <- find_GRanges_cols(df_colnames)
    checkIdentical(target, current)

    df_colnames <- c("Seqnames", "strand", "start", "end")
    target <- c(seqnames=1L, start=3L, end=4L,
                width=NA_integer_, strand=2L)
    current <- find_GRanges_cols(df_colnames)
    checkIdentical(target, current)

    df_colnames <- c("Seqnames", "strand", "start", "end")
    target <- c(seqnames=1L, start=3L, end=4L,
                width=NA_integer_, strand=NA_integer_)
    current <- find_GRanges_cols(df_colnames, ignore.strand=TRUE)
    checkIdentical(target, current)

    df_colnames <- c("seqname", "start", "end")
    target <- c(seqnames=1L, start=2L, end=3L,
                width=NA_integer_, strand=NA_integer_)
    current <- find_GRanges_cols(df_colnames)
    checkIdentical(target, current)

    df_colnames <- c("chrom", "strand", "txStart", "txEnd")
    target <- c(seqnames=1L, start=3L, end=4L,
                width=NA_integer_, strand=2L)
    current <- find_GRanges_cols(df_colnames)
    checkIdentical(target, current)

    df_colnames <- c("chrom", "strand", "txStart", "txEnd", "txChrom")
    target <- c(seqnames=5L, start=3L, end=4L,
                width=NA_integer_, strand=2L)
    current <- find_GRanges_cols(df_colnames)
    checkIdentical(target, current)

    df_colnames <- c("strand", "txStart", "txEnd", "txChrom", "txStrand")
    target <- c(seqnames=4L, start=2L, end=3L,
                width=NA_integer_, strand=5L)
    current <- find_GRanges_cols(df_colnames)
    checkIdentical(target, current)

    df_colnames <- c("chrom", "txStrand", "start", "end")
    target <- c(seqnames=1L, start=3L, end=4L,
                width=NA_integer_, strand=NA_integer_)
    current <- find_GRanges_cols(df_colnames)
    checkIdentical(target, current)

    df_colnames <- c("stop", "txEnd", "txStart", "CHR", "start")
    target <- c(seqnames=4L, start=5L, end=1L,
                width=NA_integer_, strand=NA_integer_)
    current <- find_GRanges_cols(df_colnames)
    checkIdentical(target, current)

    df_colnames <- c("txEnd", "txStart", "chromosome_name")
    target <- c(seqnames=3L, start=2L, end=1L,
                width=NA_integer_, strand=NA_integer_)
    current <- find_GRanges_cols(df_colnames)
    checkIdentical(target, current)

    df_colnames <- c("tx_end", "tx_start", "chrom", "tx_chrom")
    target <- c(seqnames=4L, start=2L, end=1L,
                width=NA_integer_, strand=NA_integer_)
    current <- find_GRanges_cols(df_colnames)
    checkIdentical(target, current)

    df_colnames <- c("chrom", "strand", "exon_chrom_start", "exon_chrom_end")
    target <- c(seqnames=1L, start=3L, end=4L,
                width=NA_integer_, strand=2L)
    current <- find_GRanges_cols(df_colnames)
    checkIdentical(target, current)

    df_colnames <- c("chrom", "start", "end", "end")
    checkException(find_GRanges_cols(df_colnames), silent=TRUE)

    df_colnames <- c("chrom", "start", "End", "end")
    checkException(find_GRanges_cols(df_colnames), silent=TRUE)

    df_colnames <- c("chrom", "start", "end", "stop")
    checkException(find_GRanges_cols(df_colnames), silent=TRUE)
    target <- c(seqnames=1L, start=2L, end=4L,
                width=NA_integer_, strand=NA_integer_)
    current <- find_GRanges_cols(df_colnames, end.field="stop")
    checkIdentical(target, current)
    checkException(find_GRanges_cols(df_colnames, end.field=4), silent=TRUE)

    df_colnames <- c("chrom", "start", "start", "stop")
    checkException(find_GRanges_cols(df_colnames), silent=TRUE)

    df_colnames <- c("chrom", "start", "end", "CHR")
    checkException(find_GRanges_cols(df_colnames), silent=TRUE)
    target <- c(seqnames=4L, start=2L, end=3L,
                width=NA_integer_, strand=NA_integer_)
    current <- find_GRanges_cols(df_colnames, seqnames.field="chr")
    checkIdentical(target, current)

    df_colnames <- c("chrom", "start", "end", "chromosome_name")
    checkException(find_GRanges_cols(df_colnames), silent=TRUE)

    df_colnames <- c("chrom", "tx_start", "tx_end", "exon_start")
    checkException(find_GRanges_cols(df_colnames), silent=TRUE)
    target <- c(seqnames=1L, start=2L, end=3L,
                width=NA_integer_, strand=NA_integer_)
    current <- find_GRanges_cols(df_colnames,
                                 start.field="tx_start", end.field="tx_end")
    checkIdentical(target, current)

    df_colnames <- c("chrom", "tx_start", "tx_end", "exon_stop")
    checkException(find_GRanges_cols(df_colnames), silent=TRUE)

    df_colnames <- c("chrom", "tx_start", "tx_end", "tx_end")
    checkException(find_GRanges_cols(df_colnames), silent=TRUE)

    df_colnames <- c("chrom", "tx_start", "tx_end", "tx_stop")
    checkException(find_GRanges_cols(df_colnames), silent=TRUE)

    df_colnames <- c("chrom", "tx_start", "start", "end")
    target <- c(seqnames=1L, start=3L, end=4L,
                width=NA_integer_, strand=NA_integer_)
    current <- find_GRanges_cols(df_colnames)
    checkIdentical(target, current)

    df_colnames <- c("chrom", "strand", "start", "end", "STRAND")
    checkException(find_GRanges_cols(df_colnames), silent=TRUE)
    target <- c(seqnames=1L, start=3L, end=4L,
                width=NA_integer_, strand=NA_integer_)
    current <- find_GRanges_cols(df_colnames, ignore.strand=TRUE)
    checkIdentical(target, current)
}

