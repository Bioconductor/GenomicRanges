### =========================================================================
### Some low-level (non exported) utility functions.
### -------------------------------------------------------------------------


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Used by "show" methods.
###

### 'makeNakedMat.FUN' must be a function returning a character matrix.
makePrettyMatrixForCompactPrinting <- function(x, makeNakedMat.FUN)
{
    nhalf <- 5L
    lx <- length(x)
    if (is.null(nhead <- getOption("showHeadLines")))
        nhead <- nhalf 
    if (is.null(ntail <- getOption("showTailLines")))
        ntail <- nhalf 

    if (lx < (nhalf*2 + 1L) | (lx < (nhead + ntail + 1L))) {
        ans <- makeNakedMat.FUN(x)
        ans_rownames <- .rownames(names(x), lx)
    } else {
        top_idx <- 1:nhead
        if (nhead == 0)
            top_idx <- 0 
        bottom_idx=(lx - (ntail - 1L)):lx
        if (ntail == 0)
            bottom_idx <- 0 
        ans_top <- makeNakedMat.FUN(x[top_idx])
        ans_bottom <- makeNakedMat.FUN(x[bottom_idx])
        ans <- rbind(ans_top,
                     matrix(rep.int("...", ncol(ans_top)), nrow=1L),
                     ans_bottom)
        ans_rownames <- .rownames(names(x), lx, top_idx, bottom_idx)
    }
    rownames(ans) <- format(ans_rownames, justify="right")
    ans
}

.rownames <- function(names=NULL, len=NULL, tindex=NULL, bindex=NULL)
{
    if (is.null(tindex) & is.null(bindex)) {
        ## all lines
        if (len == 0L)
            character(0)
        else if (is.null(names))
            paste0("[", seq_len(len), "]")
        else
            names
    } else {
        ## head and tail 
        if (!is.null(names)) {
            c(names[tindex], "...", names[bindex])
        } else {
            s1 <- paste0("[", tindex, "]")
            s2 <- paste0("[", bindex, "]")
            if (all(tindex == 0)) 
                s1 <- character(0) 
            if (all(bindex == 0)) 
                s2 <- character(0) 
            c(s1, "...", s2)
        }
    }
}

makeClassinfoRowForCompactPrinting <- function(x, col2class)
{
    ans_names <- names(col2class)
    no_bracket <- ans_names == ""
    ans_names[no_bracket] <- col2class[no_bracket]
    left_brackets <- right_brackets <- character(length(col2class))
    left_brackets[!no_bracket] <- "<"
    right_brackets[!no_bracket] <- ">"
    ans <- paste0(left_brackets, col2class, right_brackets)
    names(ans) <- ans_names
    if (ncol(mcols(x)) > 0L) {
        tmp <- sapply(mcols(x),
                      function(xx) paste0("<", classNameForDisplay(xx), ">"))
        ans <- c(ans, `|`="|", tmp)
    }
    matrix(ans, nrow=1L, dimnames=list("", names(ans)))
}

showSeqlengths <- function(object, margin="")
{
    seqlens <- seqlengths(object)
    nseq <- length(seqlens)
    halfWidth <- getOption("width") %/% 2L
    first <- max(1L, halfWidth)
    showMatrix <-
      rbind(as.character(head(names(seqlens), first)),
            as.character(head(seqlens, first)))
    if (nseq > first) {
        last <- min(nseq - first, halfWidth)
        showMatrix <-
          cbind(showMatrix,
                rbind(as.character(tail(names(seqlens), last)),
                      as.character(tail(seqlens, last))))
    }
    showMatrix <- format(showMatrix, justify="right")
    cat(margin, "seqlengths:\n", sep="")
    cat(margin, IRanges:::labeledLine("", showMatrix[1L, ], count=FALSE,
                                      labelSep=""), sep="")
    cat(margin, IRanges:::labeledLine("", showMatrix[2L, ], count=FALSE,
                                      labelSep=""), sep="")
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Used by "elementMetadata<-" methods.
###

### Try to turn 'value' into a DataFrame compatible with 'x'.
mk_elementMetadataReplacementValue <- function(x, value)
{
    if (is.null(value))
        return(new("DataFrame", nrows=length(x)))
    if (!is(value, "DataFrame"))
        value <- DataFrame(value)
    if (!is.null(rownames(value)))
        rownames(value) <- NULL
    n <- length(x)
    k <- nrow(value)
    if (k == n)
        return(value)
    if ((k == 0L) || (k > n) || (n %% k != 0L))
        stop(k, " rows in value to replace ", n, " rows")
    value[rep(seq_len(k), length.out=n), , drop=FALSE]
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Other stuff...
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

