### =========================================================================
### The seqinfo() (and releated) generic getters and setters
### -------------------------------------------------------------------------


.reverseNew2old <- function(new2old, new_N, old_N)
{
    if (is.null(new2old))
        return(NULL)
    if (!is.integer(new2old) || length(new2old) != new_N)
        stop("when 'new2old' is not NULL, it must be an integer ",
             "vector of the\n  same length as the supplied 'seqlevels'")
    min_new2old <- suppressWarnings(min(new2old, na.rm=TRUE))
    if (min_new2old != Inf) {
        if (min_new2old < 1L || max(new2old, na.rm=TRUE) > old_N)
            stop("non-NA values in 'new2old' must be >= 1 and <= N, ",
                 "where N is the\n  nb of sequence levels in 'x'")
    }
    if (any(duplicated(new2old) & !is.na(new2old)))
        stop("duplicates are not allowed among non-NA values in 'new2old'")
    IRanges:::reverseIntegerInjection(new2old, old_N)
}

### The dangling seqlevels in 'x' are those seqlevels that the user wants to
### drop but they are in use.
getDanglingSeqlevels <- function(x, new2old=NULL, force=FALSE, new_seqlevels)
{
    if (!is.character(new_seqlevels) || any(is.na(new_seqlevels)))
        stop("the supplied 'seqlevels' must be a character vector with no NAs")
    if (!isTRUEorFALSE(force))
        stop("'force' must be TRUE or FALSE")
    if (is.null(new2old))
        return(character(0))
    new_N <- length(new_seqlevels)
    old_N <- length(seqlevels(x))
    old2new <- .reverseNew2old(new2old, new_N, old_N)
    dangling_seqlevels <- intersect(unique(seqnames(x)),
                                    seqlevels(x)[is.na(old2new)])
    if (!force && length(dangling_seqlevels) != 0L)
        stop("won't drop seqlevels currently in use (",
             paste(dangling_seqlevels, collapse = ", "),
             "), unless you\n",
             "  use 'force=TRUE' to also drop elements in 'x' ",
             "where those seqlevels are used\n",
             "  (e.g. with 'seqlevels(x, force=TRUE) <- new_seqlevels').\n",
             "  Alternatively, you can also subset 'x' first.")
    dangling_seqlevels
}

### Compute the new seqnames resulting from new seqlevels.
### Assumes that 'seqnames(x)' is a 'factor' Rle (which is true if 'x' is a
### GRanges or GappedAlignments object, but not if it's a GRangesList object),
### and returns a 'factor' Rle of the same length (and same runLength vector).
makeNewSeqnames <- function(x, new2old=NULL, new_seqlevels)
{
    if (!is.character(new_seqlevels) || any(is.na(new_seqlevels)))
        stop("the supplied 'seqlevels' must be a character vector with no NAs")
    new_N <- length(new_seqlevels)
    old_N <- length(seqlevels(x))
    x_seqnames <- seqnames(x)
    if (is.null(new2old)) {
        if (new_N < old_N ||
            !identical(new_seqlevels[seq_len(old_N)], seqlevels(x)))
            stop("when 'new2old' is NULL, the first elements in the\n",
                 "  supplied 'seqlevels' must be identical to 'seqlevels(x)'")
        levels(x_seqnames) <- new_seqlevels
        return(x_seqnames)
    }
    old2new <- .reverseNew2old(new2old, new_N, old_N)
    tmp <- runValue(x_seqnames)
    levels(tmp) <- new_seqlevels[old2new]
    runValue(x_seqnames) <- factor(as.character(tmp), levels=new_seqlevels)
    return(x_seqnames)
}

### Returns -2 for "subsetting" mode, -1 for "renaming" mode, or an integer
### vector containing the mapping from the new to the old levels for "general"
### mode (i.e. a combination of renaming and/or subsetting). Note that this
### integer vector is guaranteed to contain no negative values.
getSeqlevelsReplacementMode <- function(new_seqlevels, old_seqlevels)
{
    ## Only check for NAs. Duplicated or zero-length values will be rejected
    ## later by the Seqinfo() constructor.
    if (!is.character(new_seqlevels) || any(is.na(new_seqlevels)))
        stop("the supplied 'seqlevels' must be a character vector with no NAs")
    nsl_names <- names(new_seqlevels)
    if (!is.null(nsl_names)) {
        nonempty_names <- nsl_names[!(nsl_names %in% c(NA, ""))]
        if (any(duplicated(nonempty_names)) ||
            length(setdiff(nonempty_names, old_seqlevels)) != 0L)
            stop("the names of the supplied 'seqlevels' contain duplicates ",
                 "or invalid sequence levels")
        return(match(nsl_names, old_seqlevels))
    }
    if (length(new_seqlevels) != length(old_seqlevels))
        return(-2L)
    is_renamed <- new_seqlevels != old_seqlevels
    tmp <- intersect(new_seqlevels[is_renamed], old_seqlevels[is_renamed])
    if (length(tmp) != 0L)
        return(-2L)
    return(-1L)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### seqinfo() getter and setter.
###

setGeneric("seqinfo", function(x) standardGeneric("seqinfo"))

setGeneric("seqinfo<-", signature="x",
    function(x, new2old=NULL, force=FALSE, value) standardGeneric("seqinfo<-")
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### seqnames() getter and setter.
###

setGeneric("seqnames", function(x) standardGeneric("seqnames"))

setGeneric("seqnames<-", signature="x",
    function(x, value) standardGeneric("seqnames<-")
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### seqlevels() getter and setter.
###

setGeneric("seqlevels", function(x) standardGeneric("seqlevels"))

### Default "seqlevels" method works on any object 'x' with a working
### "seqinfo" method.
setMethod("seqlevels", "ANY", function(x) seqlevels(seqinfo(x)))

setGeneric("seqlevels<-", signature="x",
    function(x, force=FALSE, value) standardGeneric("seqlevels<-")
)

### Default "seqlevels<-" method works on any object 'x' with working
### "seqinfo" and "seqinfo<-" methods.
setReplaceMethod("seqlevels", "ANY",
    function(x, force=FALSE, value)
    {
        ## Make the new Seqinfo object.
        x_seqinfo <- seqinfo(x)
        seqlevels(x_seqinfo) <- value
        ## Map the new sequence levels to the old ones.
        new2old <- getSeqlevelsReplacementMode(value, seqlevels(x))
        if (identical(new2old, -2L)) {
            new2old <- match(value, seqlevels(x))
        } else if (identical(new2old, -1L)) {
            new2old <- seq_len(length(value))
        }
        ## Do the replacement.
        seqinfo(x, new2old=new2old, force=force) <- x_seqinfo
        x
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### seqlengths() getter and setter.
###

setGeneric("seqlengths", function(x) standardGeneric("seqlengths"))

### Default "seqlengths" method works on any object 'x' with a working
### "seqinfo" method.
setMethod("seqlengths", "ANY", function(x) seqlengths(seqinfo(x)))

setGeneric("seqlengths<-", signature="x",
    function(x, value) standardGeneric("seqlengths<-")
)

### Default "seqlengths<-" method works on any object 'x' with working
### "seqinfo" and "seqinfo<-" methods.
setReplaceMethod("seqlengths", "ANY",
    function(x, value)
    {
        seqlengths(seqinfo(x)) <- value
        x
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### isCircular() getter and setter.
###

setGeneric("isCircular", function(x) standardGeneric("isCircular"))

### Default "isCircular" method works on any object 'x' with a working
### "seqinfo" method.
setMethod("isCircular", "ANY", function(x) isCircular(seqinfo(x)))

setGeneric("isCircular<-", signature="x",
    function(x, value) standardGeneric("isCircular<-")
)

### Default "isCircular<-" method works on any object 'x' with working
### "seqinfo" and "seqinfo<-" methods.
setReplaceMethod("isCircular", "ANY",
    function(x, value)
    {
        isCircular(seqinfo(x)) <- value
        x
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### genome() getter and setter.
###

setGeneric("genome", function(x) standardGeneric("genome"))

### Default "genome" method works on any object 'x' with a working
### "seqinfo" method.
setMethod("genome", "ANY", function(x) genome(seqinfo(x)))

setGeneric("genome<-", signature="x",
    function(x, value) standardGeneric("genome<-")
)

### Default "genome<-" method works on any object 'x' with working
### "seqinfo" and "seqinfo<-" methods.
setReplaceMethod("genome", "ANY",
    function(x, value)
    {
        genome(seqinfo(x)) <- value
        x
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### makeSeqnameIds().
###
### Assign a unique ID to each unique sequence name passed in 'seqnames'.
### The returned IDs are guaranteed to span 1:N where N is the number of
### unique sequence names in 'seqnames'.
### Also the function tries hard to assign IDs in a way that is consistent
### with a "good looking" order of the sequence names. This "good looking"
### order is roughly defined by the following complicated and arbitrary set
### of rules (rules apply in the order shown below):
###
###   1. Names starting with CHR < names starting with chr <
###      names starting with CH < names starting with ch.
###
###   2. Names where an arabic number is followed by an L or R will be placed
###      in the group of sequences with an arabic number, right after the
###      sequence with that number.
###
###   3. Every name should fall into exactly 1 of the following 14 groups:
###        (a) roman numbers
###        (b) arabic numbers (possibly followed by L or R)
###        (c) X
###        (d) Y
###        (e) U
###        (f) M
###        (g) MT
###        (h) arabic numbers "followed by something" (not L or R)
###        (i) X "followed by something"
###        (j) Y "followed by something"
###        (k) U "followed by something"
###        (l) M "followed by something"
###        (m) MT "followed by something"
###        (n) anything else
###      Names in first groups are < names in last groups.
###
###   4. In groups (h) to (n), ties are broken by looking at the "followed
###      by something" part (or at the entire name for group (n)): collation
###      defined by LC_COLLATE set to C applies.
###
### 'X.is.seXual' lets the user control whether X refers to the sexual
### chromosome or to chromosome with roman number X.
###
### One of the ugliest functions I ever wrote, sorry...
### NOTE: makeSeqnameIds() returns IDs that are consistent with the order of
### the seqnames in BSgenome.Hsapiens.UCSC.hg19, BSgenome.Celegans.UCSC.ce2,
### BSgenome.Dmelanogaster.UCSC.dm3 and BSgenome.Scerevisiae.UCSC.sacCer2,
### i.e. 'makeSeqnameIds(seqnames(Hsapiens))' is identical to
### 'seq_len(length(seqnames(Hsapiens))))', etc...

makeSeqnameIds <- function(seqnames, X.is.seXual=NA)
{
    if (!is.character(seqnames))
        stop("'seqnames' must be a character vector")
    if (!is.logical(X.is.seXual) || length(X.is.seXual) != 1L)
        stop("'X.is.seXual' must be a single logical")

    seqnames_sm <- match(seqnames, seqnames)
    uidx <- which(seqnames_sm == seq_len(length(seqnames_sm)))
    useqnames <- seqnames[uidx]

    ## Set LC_COLLATE to C so the order of the factor levels will be the same
    ## for everybody (i.e. for every user on any machine in any country).
    prev_locale <- Sys.getlocale("LC_COLLATE")
    Sys.setlocale("LC_COLLATE", "C")
    on.exit(Sys.setlocale("LC_COLLATE", prev_locale))

    ## Provisional ids.
    prov_ids <- rep.int(NA_integer_, length(useqnames))
    last_prov_id <- 0L

    ## 'i' indices of elements in 'prov_ids' to set.
    ## 'ints' integer vector of length > 0 with no NAs. Is recycled to the
    ## length of 'i'.
    assignNewProvIds <- function(i, ints=0L)
    {
        new_prov_ids <- last_prov_id + ints - min(ints) + 1L
        prov_ids2 <- prov_ids
        prov_ids2[i] <- new_prov_ids
        last_prov_id2 <- max(new_prov_ids)
        assign("prov_ids", prov_ids2, inherits=TRUE)
        assign("last_prov_id", last_prov_id2, inherits=TRUE)
    }

    hasPrefix <- function(seqnames, prefix)
    {
        substr(seqnames, start=1L, stop=nchar(prefix)) == prefix
    }

    REGEXP0 <- "[1-9][0-9]*"
    REGEXP1 <- "[0-9]*"
    isNb <- function(seqnames, abc="")
    {
        regexp <- paste0("^", REGEXP0, abc, "$")
        grepl(regexp, seqnames)
    }
    getIntPrefix <- function(seqnames)
    {
        regexp <- paste0("^(", REGEXP1, ")(.*)$")
        sub(regexp, "\\1", seqnames)
    }
    getIntPrefixTail <- function(seqnames)
    {
        regexp <- paste0("^(", REGEXP1, ")(.*)$")
        sub(regexp, "\\2", seqnames)
    }

    assignProvIdsToUniqueSeqnames <- function(useqnames, prefix)
    {
        idx <- which(is.na(prov_ids) & hasPrefix(useqnames, prefix))
        if (length(idx) == 0L)
            return()
        useqnames_suffix <- useqnames[idx]
        useqnames_suffix <- substr(useqnames_suffix,
                                   start=nchar(prefix)+1L,
                                   stop=nchar(useqnames_suffix))
        suppressWarnings(suff_as_roman <- as.roman(useqnames_suffix))
        suff_is_nb <- isNb(useqnames_suffix)
        suff_is_nbL <- isNb(useqnames_suffix, abc="L")
        suff_is_nbR <- isNb(useqnames_suffix, abc="R")
        suff_is_nb_or_nbL_or_nbR <- suff_is_nb | suff_is_nbL | suff_is_nbR
        suff_is_nbother <- isNb(useqnames_suffix, abc=".*") &
                           !suff_is_nb_or_nbL_or_nbR
        suff_is_X <- useqnames_suffix == "X"
        suff_is_Y <- useqnames_suffix == "Y"
        if (is.na(X.is.seXual)) {
            if (any(suff_is_X)) {
                X.is.seXual <- any(suff_is_Y) ||
                               any(suff_is_nb_or_nbL_or_nbR) ||
                               any(suff_is_nbother)
            } else {
                X.is.seXual <- TRUE  # or FALSE, won't make any difference
            }
        }
        suff_is_seXual <- suff_is_X & X.is.seXual
        suff_is_roman <- !is.na(suff_as_roman) & !suff_is_seXual
        suff_is_U <- useqnames_suffix == "U"
        suff_is_MT <- useqnames_suffix == "MT"
        suff_is_M <- useqnames_suffix == "M"
        suff_is_Xother <- hasPrefix(useqnames_suffix, "X") & !suff_is_X &
                          !suff_is_roman
        suff_is_Yother <- hasPrefix(useqnames_suffix, "Y") & !suff_is_Y &
                          !suff_is_roman
        suff_is_Uother <- hasPrefix(useqnames_suffix, "U") & !suff_is_U &
                          !suff_is_roman
        suff_is_MTother <- hasPrefix(useqnames_suffix, "MT") & !suff_is_MT &
                           !suff_is_roman
        suff_is_Mother <- hasPrefix(useqnames_suffix, "M") & !suff_is_M &
                          !suff_is_MT & !suff_is_MTother & !suff_is_roman
        ## The groups below must define a partitioning of 'useqnames_suffix'
        ## i.e. any element in 'useqnames_suffix' must fall in exactly 1 group.
        suff_is_other <- !suff_is_roman & !suff_is_nb_or_nbL_or_nbR &
                         !suff_is_seXual & !suff_is_Y &
                         !suff_is_U & !suff_is_M & !suff_is_MT &
                         !suff_is_nbother &
                         !suff_is_Xother & !suff_is_Yother &
                         !suff_is_Uother & !suff_is_Mother & !suff_is_MTother
        if (any(suff_is_roman)) {
            ints <- as.integer(suff_as_roman[suff_is_roman])
            assignNewProvIds(idx[suff_is_roman], ints=ints)
        }
        if (any(suff_is_nb_or_nbL_or_nbR)) {
            suff0 <- useqnames_suffix[suff_is_nb_or_nbL_or_nbR]
            is_nb <- suff_is_nb[suff_is_nb_or_nbL_or_nbR]
            is_nbL <- suff_is_nbL[suff_is_nb_or_nbL_or_nbR]
            is_nbR <- suff_is_nbR[suff_is_nb_or_nbL_or_nbR]
            nb_ints <- as.integer(suff0[is_nb])
            suff0L <- suff0[is_nbL]
            nbL_ints <- as.integer(substr(suff0L,
                                          start=1L,
                                          stop=nchar(suff0L)-1L))
            suff0R <- suff0[is_nbR]
            nbR_ints <- as.integer(substr(suff0R,
                                          start=1L,
                                          stop=nchar(suff0R)-1L))
            ints <- integer(length(suff0))
            ints[is_nb] <- 3L * nb_ints
            ints[is_nbL] <- 3L * nbL_ints + 1L
            ints[is_nbR] <- 3L * nbR_ints + 2L
            assignNewProvIds(idx[suff_is_nb_or_nbL_or_nbR], ints=ints)
        }
        if (any(suff_is_seXual))
            assignNewProvIds(idx[suff_is_seXual])
        if (any(suff_is_Y))
            assignNewProvIds(idx[suff_is_Y])
        if (any(suff_is_U))
            assignNewProvIds(idx[suff_is_U])
        if (any(suff_is_M))
            assignNewProvIds(idx[suff_is_M])
        if (any(suff_is_MT))
            assignNewProvIds(idx[suff_is_MT])
        if (any(suff_is_nbother)) {
            suff0 <- useqnames_suffix[suff_is_nbother]
            ints1 <- as.integer(getIntPrefix(suff0))
            ints2 <- as.integer(factor(getIntPrefixTail(suff0)))
            ints <- (max(ints2) + 1L) * ints1 + ints2
            assignNewProvIds(idx[suff_is_nbother], ints=ints)
        }
        if (any(suff_is_Xother)) {
            suff0 <- useqnames_suffix[suff_is_Xother]
            ints <- as.integer(factor(substr(suff0,
                                             start=2L,
                                             stop=nchar(suff0))))
            assignNewProvIds(idx[suff_is_Xother], ints=ints)
        }
        if (any(suff_is_Yother)) {
            suff0 <- useqnames_suffix[suff_is_Yother]
            ints <- as.integer(factor(substr(suff0,
                                             start=2L,
                                             stop=nchar(suff0))))
            assignNewProvIds(idx[suff_is_Yother], ints=ints)
        }
        if (any(suff_is_Uother)) {
            suff0 <- useqnames_suffix[suff_is_Uother]
            ints <- as.integer(factor(substr(suff0,
                                             start=2L,
                                             stop=nchar(suff0))))
            assignNewProvIds(idx[suff_is_Uother], ints=ints)
        }
        if (any(suff_is_Mother)) {
            suff0 <- useqnames_suffix[suff_is_Mother]
            ints <- as.integer(factor(substr(suff0,
                                             start=3L,
                                             stop=nchar(suff0))))
            assignNewProvIds(idx[suff_is_Mother], ints=ints)
        }
        if (any(suff_is_MTother)) {
            suff0 <- useqnames_suffix[suff_is_MTother]
            ints <- as.integer(factor(substr(suff0,
                                             start=2L,
                                             stop=nchar(suff0))))
            assignNewProvIds(idx[suff_is_MTother], ints=ints)
        }
        if (any(suff_is_other)) {
            suff0 <- useqnames_suffix[suff_is_other]
            ints <- as.integer(factor(suff0))
            assignNewProvIds(idx[suff_is_other], ints=ints)
        }
    }

    ## Longest prefixes first.
    assignProvIdsToUniqueSeqnames(useqnames, "CHR")
    assignProvIdsToUniqueSeqnames(useqnames, "chr")
    assignProvIdsToUniqueSeqnames(useqnames, "CH")
    assignProvIdsToUniqueSeqnames(useqnames, "ch")
    assignProvIdsToUniqueSeqnames(useqnames, "")

    ## Unique ids.
    uids <- integer(length(prov_ids))
    oo <- order(prov_ids)
    uids[oo] <- seq_len(length(uids))

    ## Assign ids to all seqnames.
    uids[seqnames_sm]
}

