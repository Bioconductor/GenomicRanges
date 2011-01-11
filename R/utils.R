### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Some low-level (non exported) utility functions.
###

### Note that, strictly speaking, mergeNamedAtomicVectors() is not
### commutative, i.e., in general 'z1 <- mergeNamedAtomicVectors(x, y)' is
### not identical to 'z2 <- mergeNamedAtomicVectors(y, x)'. However 'z1' and
### 'z2' are both guaranteed to have unique names and to contain the same set
### of name/value pairs (but typically not in the same order).
mergeNamedAtomicVectors <- function(x, y, what=c("key", "values"))
{
    if (!is.atomic(x) || !is.atomic(y) || typeof(x) != typeof(y))
        stop("'x' and 'y' must be atomic vectors of the same type")
    x_names <- names(x)
    y_names <- names(y)
    if (is.null(x_names) || is.null(y_names))
        stop("'x' and 'y' must have names")
    ans_names <- union(x_names, y_names)
    if (any(ans_names %in% c(NA_character_, "")))
        stop("some names in 'x' or 'y' are NA or the empty string")
    ans <- x[ans_names]
    ## Some of 'ans' names can be NA. This is because subsetting a named
    ## vector 'x' by a character vector 'i' returns a vector whose names
    ## are 'i' except for the values in 'i' that are not in 'names(x)'.
    ## So we need to fix this.
    names(ans) <- ans_names
    ans2 <- y[ans_names]
    idx <- which(ans != ans2)
    if (length(idx) != 0L) {
        msg <- c(what[1L], ifelse(length(idx) >= 2, "s", ""), " ",
                 paste(ans_names[idx], collapse=", "), " ",
                 ifelse(length(idx) >= 2, "have", "has"),
                 " incompatible ", what[2L], ":\n  - in 'x': ",
                 paste(ans[idx], collapse=", "), "\n  - in 'y': ",
                 paste(ans2[idx], collapse=", "))
        stop(msg)
    }
    idx <- is.na(ans) & !is.na(ans2)
    ans[idx] <- ans2[idx]
    ans
}

