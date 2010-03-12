### The strand() generic and methods.

setGeneric("strand", function(x) standardGeneric("strand"))
setGeneric("strand<-", function(x, value) standardGeneric("strand<-"))

setMethod("strand", "missing", function(x) factor(levels=c("+","-","*")))

setMethod("strand", "character",
    function(x) {
        lvls <- levels(strand())
        if (!all(is.na(x) | (x %in% lvls)))
            stop("strand values must be in '", paste(lvls, collapse="' '"), "'")
        factor(x, levels=lvls)
    })

setMethod("strand", "logical",
    function(x) {
        ans <- strand()
        length(ans) <- length(x)
        ans[!x] <- "+"
        ans[x] <- "-"
        ans
    })

setMethod("strand", "DataTable",
    function(x) {
        ans <- x[["strand"]]
        if (is.null(ans))
            ans <- strand(rep(NA_character_, nrow(x)))
        else if (is.character(ans))
            ans <- strand(ans)
        else if (is(ans, "Rle"))
            ans <- as.vector(ans)
        ans
    })


setGeneric("compatableStrand", signature=c("x","y"),
    function(x, y) standardGeneric("compatableStrand")
)

setMethod("compatableStrand", c("factor", "factor"),
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
        if (IRanges:::anyMissing(ans)) {
            fix <- IRanges:::whichAsVector(is.na(ans))
            ans[fix] <- (x[fix] == 0L) | (y[fix] == 0L)
            if (IRanges:::anyMissing(ans))
                ans[is.na(ans)] <- FALSE
        }
        ans
    }
)

setMethod("compatableStrand", c("Rle", "Rle"),
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
        if (IRanges:::anyMissing(runValue(ans))) {
            fix <- IRanges:::whichAsVector(is.na(runValue(ans)))
            runValue(ans)[fix] <-
              (runValue(x) == 0L)[fix] | (runValue(y)[fix] == 0L)
            if (IRanges:::anyMissing(runValue(ans)))
                runValue(ans)[is.na(runValue(ans))] <- FALSE
        }
        ans
    }
)
