### =========================================================================
### makeGPosFromDataFrame()
### -------------------------------------------------------------------------

### 'df' must be a data.frame or DataFrame object.
makeGPosFromDataFrame <- function(df,
                                  keep.extra.columns=FALSE,
                                  ignore.strand=FALSE,
                                  seqinfo=NULL,
                                  seqnames.field=c("seqnames", "seqname",
                                                   "chromosome", "chrom",
                                                   "chr", "chromosome_name",
                                                   "seqid"),
                                  start.field=c("start", "pos"),
                                  end.field=c("end", "stop", "pos"),
                                  strand.field="strand",
                                  starts.in.df.are.0based=FALSE)
{
    .makeXFromDataFrame(df = df, x = "GPos",
        keep.extra.columns=keep.extra.columns,
        ignore.strand=ignore.strand,
        seqinfo=seqinfo,
        seqnames.field=seqnames.field,
        start.field=start.field,
        end.field=end.field,
        strand.field=strand.field,
        starts.in.df.are.0based=starts.in.df.are.0based)
}

setAs("data.frame", "GPos",
    function(from) makeGPosFromDataFrame(from, keep.extra.columns=TRUE)
)

setAs("DataFrame", "GPos",
    function(from) makeGPosFromDataFrame(from, keep.extra.columns=TRUE)
)

