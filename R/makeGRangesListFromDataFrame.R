### =========================================================================
### makeGRangesListFromDataFrame()
### -------------------------------------------------------------------------

### 'df' must be a data.frame or DataFrame object.
makeGRangesListFromDataFrame <-
    function(df,
             split.field = NULL,
             names.field = NULL,
             ...)
{
    splitIdx <- namesIdx <- integer()
    if (!is.null(split.field)) {
        if (!isSingleString(split.field))
            stop("'split.field' must be a single string")
        splitIdx <- which(names(df) %in% split.field)
        if (!length(splitIdx))
            stop("'split.field' is not in 'names(df)'")
        if (length(splitIdx) > 1L)
            stop("'split.field' matched more than one 'names(df)'")
        splitField <- df[[split.field]]
    } else splitField <- seq_len(nrow(df))

    if (!is.null(names.field)) {
        if (!isSingleString(names.field))
            stop("'names.field' must be a single string")
        namesIdx <- which(names(df) %in% names.field)
        if (!length(namesIdx))
            stop("'names.field' is not found in 'names(df)'")
        if (length(namesIdx) > 1L)
            stop("'names.field' matched more than one 'names(df)'")
        namesField <- df[[names.field]]
    } else namesField <- NULL 

    if (length(c(splitIdx, namesIdx)))
        df <- df[, -c(splitIdx, namesIdx)]
    gr <- makeGRangesFromDataFrame(df, ...)
    names(gr) <- namesField 
    S4Vectors::split(gr, splitField)
}
