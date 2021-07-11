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
### compatibleStrand() generic and methods
###

setGeneric("compatibleStrand", signature=c("x","y"),  # not exported
    function(x, y) standardGeneric("compatibleStrand")
)

setMethod("compatibleStrand", c("factor", "factor"),  # not exported
    function(x, y)
    {
        lvls <- levels(strand())
        if (length(x) != length(y))
            stop("'x' and 'y' must be of equal length")
        if (!identical(levels(x), lvls) || !identical(levels(y), lvls))
            stop("strand values must be in '", paste(lvls, collapse="' '"), "'")

        levels(x) <- c("1", "-1", "0")
        x <- as.integer(as.character(x))

        levels(y) <- c("1", "-1", "0")
        y <- as.integer(as.character(y))

        ans <- x * y != -1L
        if (S4Vectors:::anyMissing(ans)) {
            fix <- which(is.na(ans))
            ans[fix] <- (x[fix] == 0L) | (y[fix] == 0L)
            if (S4Vectors:::anyMissing(ans))
                ans[is.na(ans)] <- FALSE
        }
        ans
    }
)

setMethod("compatibleStrand", c("Rle", "Rle"),  # not exported
    function(x, y)
    {
        lvls <- levels(strand())
        if (length(x) != length(y))
            stop("'x' and 'y' must be of equal length")
        if (!identical(levels(runValue(x)), lvls) ||
            !identical(levels(runValue(y)), lvls))
            stop("strand values must be in '", paste(lvls, collapse="' '"), "'")

        levels(x) <- c("1", "-1", "0")
        runValue(x) <- as.integer(as.character(runValue(x)))

        levels(y) <- c("1", "-1", "0")
        runValue(y) <- as.integer(as.character(runValue(y)))

        ans <- x * y != -1L
        if (S4Vectors:::anyMissing(runValue(ans))) {
            fix <- which(is.na(runValue(ans)))
            runValue(ans)[fix] <-
              (runValue(x) == 0L)[fix] | (runValue(y)[fix] == 0L)
            if (S4Vectors:::anyMissing(runValue(ans)))
                runValue(ans)[is.na(runValue(ans))] <- FALSE
        }
        ans
    }
)

