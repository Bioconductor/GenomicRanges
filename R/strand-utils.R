### =========================================================================
### Strand utilities
### -------------------------------------------------------------------------


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Some "strand" methods
###

setMethod("strand", "missing", function(x) factor(levels=c("+","-","*")))
setMethod("strand", "NULL", function(x) strand())

setMethod("strand", "character",
    function(x)
    {
        lvls <- levels(strand())
        if (!all(x %in% lvls))
            stop("strand values must be in '", paste(lvls, collapse="' '"), "'")
        factor(x, levels=lvls)
    }
)

setMethod("strand", "factor",
    function(x)
    {
        if (any(is.na(x)))
            stop("NA not a valid strand value, use \"*\" instead")
        lvls <- levels(strand())
        x_levels <- levels(x)
        if (identical(x_levels, lvls))
            return(x)
        invalid_levels <- setdiff(x_levels, lvls)
        if (length(invalid_levels) != 0L)
            stop("invalid strand levels in 'x': ",
                 paste(invalid_levels, collapse=", "))
        factor(x, levels=lvls)
    }
)

setMethod("strand", "integer",
    function(x)
    {
        lvls <- c(1L, -1L, NA)
        if (!all(x %in% lvls))
            stop("strand values must be in '", paste(lvls, collapse="' '"), "'")
        ans <- rep.int(strand("*"), length(x))
        ans[x ==  1L] <- "+"
        ans[x == -1L] <- "-"
        ans
    }
)

setMethod("strand", "logical",
    function(x)
    {
        ans <- rep.int(strand("*"), length(x))
        ans[!x] <- "+"
        ans[ x] <- "-"
        ans
    }
)

setMethod("strand", "Rle",
    function(x)
    {
        x_runValue <- runValue(x)
        if (!(is.character(x_runValue) ||
              is.factor(x_runValue) ||
              is.integer(x_runValue) ||
              is.logical(x_runValue)))
            stop("\"strand\" method for Rle objects only works on a ",
                 "character-, factor-, integer-, or logical-Rle object")
        runValue(x) <- strand(x_runValue)
        x
    }
)

setMethod("strand", "RleList",
    function(x) relist(strand(unlist(x, use.names=FALSE)), x)
)

setMethod("strand", "DataFrame",
    function(x)
    {
        ans <- x[["strand"]]
        if (is.null(ans)) {
            ans <- rep.int(strand("*"), nrow(x))
        } else {
            ans <- strand(ans)
        }
        ans
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Some "strand<-" methods
###

normalize_strand_replacement_value <- function(value, x)
{
    if (!is(value, "Rle"))
        value <- Rle(value)
    if (!is.factor(runValue(value))
     || !identical(levels(runValue(value)), levels(strand())))
        runValue(value) <- strand(runValue(value))
    S4Vectors:::V_recycle(value, x, x_what="value", skeleton_what="x")
}

setReplaceMethod("strand", "DataFrame",
    function(x, value)
    {
        x$strand <- normalize_strand_replacement_value(value, seq_len(nrow(x)))
        x
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Some "invertStrand" methods
###

.invert_strand_factor <- function(x)
{
    x_codes <- as.integer(x)
    switch_idx <- which(x_codes <= 2L)
    x[switch_idx] <- structure(3L - x_codes[switch_idx],
                               levels=levels(x), class=class(x))
    x
}

### One method for each of the "strand" methods defined above that return a
### factor (except for the method for "missing").
setMethods("invertStrand", list("NULL",
                                "character",
                                "factor",
                                "integer",
                                "logical"),
    function(x) .invert_strand_factor(strand(x))
)

setMethod("invertStrand", "Rle",
    function(x)
    {
        ans <- strand(x)
        runValue(ans) <- invertStrand(runValue(ans))
        ans
    }
)

setMethod("invertStrand", "RleList",
    function(x) relist(invertStrand(unlist(x, use.names=FALSE)), x)
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### compatibleStrand()
###

compatibleStrand <- function(a, b) {
    a == "*" | b == "*" | a == b
}
