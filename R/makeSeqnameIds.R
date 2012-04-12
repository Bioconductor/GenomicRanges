### =========================================================================
### makeSeqnameIds()
### -------------------------------------------------------------------------
###
### Assign a unique ID to each unique sequence name passed in 'seqnames'.
### The returned IDs are guaranteed to span 1:N where N is the number of
### unique sequence names in 'seqnames'.
### Also the function tries hard to assign IDs in a way that is consistent
### with a "good looking" order of the sequence names. This "good looking"
### order is roughly defined by the following (complicated and arbitrary)
### set of rules (rules apply in the order shown below):
###
###   1. Every name should fall into exactly 1 of the 5 following "super
###      groups":
###        (A) starts with CHR
###        (B) starts with chr
###        (C) starts with CH
###        (D) starts with ch
###        (E) anything else
###      Names in early super groups are ranked before names in late super
###      groups.
###
###   2. Within each super group, and after the prefix corresponding to the
###      super group has been dropped (nothing is dropped for super group (E)),
###      every name should fall into exactly 1 of the 14 following groups:
###        (a) roman number
###        (b) arabic number (possibly followed by L or R)
###        (c) X
###        (d) Y
###        (e) U
###        (f) M
###        (g) MT
###        (h) arabic number "followed by something" (not L or R)
###        (i) X "followed by something"
###        (j) Y "followed by something"
###        (k) U "followed by something"
###        (l) M "followed by something"
###        (m) MT "followed by something"
###        (n) anything else
###      Names in early groups are ranked before names in late groups.
###
###   3. A name in group (b) that ends with L or R is ranked right after the
###      name obtained by dropping the L or R.
###
###   4. In groups (h) to (n), ties are broken by looking at the "followed
###      by something" part (or at the entire name for group (n)): collation
###      defined by LC_COLLATE set to C applies.
###
### 'X.is.seXual' lets the user control whether X refers to the sexual
### chromosome or to chromosome with roman number X.
###
### One of the ugliest functions I ever wrote, sorry...
###
### NOTE: makeSeqnameIds() was successfully tested on the BSgenome data
### packages for hg19, mm10, ce2, dm3, sacCer1, sacCer2, sacCer3 and rheMac2
### i.e. the IDs returned on the seqnames defined in those packages match the
### ranks of the seqnames.
### For example, for hg19, 'makeSeqnameIds(seqnames(Hsapiens))' is identical
### to 'seq_len(length(seqnames(Hsapiens))))'.
### TODO: Add unit test for makeSeqnameIds().

makeSeqnameIds <- function(seqnames, X.is.seXual=NA)
{
    if (!is.character(seqnames))
        stop("'seqnames' must be a character vector")
    if (!is.logical(X.is.seXual) || length(X.is.seXual) != 1L)
        stop("'X.is.seXual' must be a single logical")

    seqnames_sm <- match(seqnames, seqnames)  # self-match
    uidx <- which(seqnames_sm == seq_len(length(seqnames_sm)))
    seqlevels <- seqnames[uidx]  # unique seqnames

    ## Set LC_COLLATE to C so the order of the factor levels will be the same
    ## for everybody (i.e. for every user on any machine in any country).
    prev_locale <- Sys.getlocale("LC_COLLATE")
    Sys.setlocale("LC_COLLATE", "C")
    on.exit(Sys.setlocale("LC_COLLATE", prev_locale))

    ## Provisional ids.
    prov_ids <- rep.int(NA_integer_, length(seqlevels))
    last_prov_id <- 0L

    ## 'i' indices of elements in 'prov_ids' to set.
    ## 'ints' integer vector of length > 0 with no NAs. Is recycled to the
    ## length of 'i'.
    makeAndAssignProvIds <- function(i, ints=0L)
    {
        new_prov_ids <- last_prov_id + ints - min(ints) + 1L
        prov_ids2 <- prov_ids
        prov_ids2[i] <- new_prov_ids
        last_prov_id2 <- max(new_prov_ids)
        assign("prov_ids", prov_ids2, inherits=TRUE)
        assign("last_prov_id", last_prov_id2, inherits=TRUE)
    }

    ## Some simple helpers for low-level string manipulation.
    hasPrefix <- function(seqnames, prefix)
        substr(seqnames, start=1L, stop=nchar(prefix)) == prefix
    dropPrefix <- function(seqnames, nchar)
        substr(seqnames, start=nchar+1L, stop=nchar(seqnames))
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

    assignProvIdsForSuperGroup <- function(seqlevels, prefix)
    {
        sgidx <- which(is.na(prov_ids) & hasPrefix(seqlevels, prefix))
        if (length(sgidx) == 0L)
            return()
        sgsuffix <- seqlevels[sgidx]
        sgsuffix <- dropPrefix(sgsuffix, nchar(prefix))
        suppressWarnings(suff_as_roman <- as.roman(sgsuffix))
        is_nb <- isNb(sgsuffix)
        is_nbL <- isNb(sgsuffix, abc="L")
        is_nbR <- isNb(sgsuffix, abc="R")
        is_nb_or_nbL_or_nbR <- is_nb | is_nbL | is_nbR
        is_nbxxx <- isNb(sgsuffix, abc=".*") & !is_nb_or_nbL_or_nbR
        is_X <- sgsuffix == "X"
        is_Y <- sgsuffix == "Y"
        if (is.na(X.is.seXual)) {
            if (any(is_X)) {
                X.is.seXual <- any(is_Y) ||
                               any(is_nb_or_nbL_or_nbR) ||
                               any(is_nbxxx)
            } else {
                X.is.seXual <- TRUE  # or FALSE, won't make any difference
            }
        }
        is_seXual <- is_X & X.is.seXual
        is_roman <- !is.na(suff_as_roman) & !is_seXual
        is_U <- sgsuffix == "U"
        is_MT <- sgsuffix == "MT"
        is_M <- sgsuffix == "M"
        is_Xxxx <- hasPrefix(sgsuffix, "X") & !is_X & !is_roman
        is_Yxxx <- hasPrefix(sgsuffix, "Y") & !is_Y & !is_roman
        is_Uxxx <- hasPrefix(sgsuffix, "U") & !is_U & !is_roman
        is_MTxxx <- hasPrefix(sgsuffix, "MT") & !is_MT & !is_roman
        is_Mxxx <- hasPrefix(sgsuffix, "M") & !is_M &
                   !is_MT & !is_MTxxx & !is_roman
        ## The groups below must define a partitioning of the current super
        ## group i.e. any element in 'sgsuffix' must fall in exactly 1 group.
        is_xxx <- !is_roman & !is_nb_or_nbL_or_nbR &
                  !is_seXual & !is_Y & !is_U & !is_M & !is_MT &
                  !is_nbxxx &
                  !is_Xxxx & !is_Yxxx & !is_Uxxx & !is_Mxxx & !is_MTxxx
        ## Group (a).
        if (any(is_roman)) {
            ints <- as.integer(suff_as_roman[is_roman])
            makeAndAssignProvIds(sgidx[is_roman], ints=ints)
        }
        ## Group (b).
        if (any(is_nb_or_nbL_or_nbR)) {
            gsuffix <- sgsuffix[is_nb_or_nbL_or_nbR]
            isnb_idx <- which(is_nb[is_nb_or_nbL_or_nbR])
            isnbL_idx <- which(is_nbL[is_nb_or_nbL_or_nbR])
            isnbR_idx <- which(is_nbR[is_nb_or_nbL_or_nbR])
            nb_ints <- as.integer(gsuffix[isnb_idx])
            gsuffixL <- gsuffix[isnbL_idx]
            nbL_ints <- as.integer(substr(gsuffixL,
                                          start=1L,
                                          stop=nchar(gsuffixL)-1L))
            gsuffixR <- gsuffix[isnbR_idx]
            nbR_ints <- as.integer(substr(gsuffixR,
                                          start=1L,
                                          stop=nchar(gsuffixR)-1L))
            ints <- integer(length(gsuffix))
            ints[isnb_idx] <- 3L * nb_ints
            ints[isnbL_idx] <- 3L * nbL_ints + 1L
            ints[isnbR_idx] <- 3L * nbR_ints + 2L
            makeAndAssignProvIds(sgidx[is_nb_or_nbL_or_nbR], ints=ints)
        }
        ## Group (c).
        if (any(is_seXual))
            makeAndAssignProvIds(sgidx[is_seXual])
        ## Group (d).
        if (any(is_Y))
            makeAndAssignProvIds(sgidx[is_Y])
        ## Group (e).
        if (any(is_U))
            makeAndAssignProvIds(sgidx[is_U])
        ## Group (f).
        if (any(is_M))
            makeAndAssignProvIds(sgidx[is_M])
        ## Group (g).
        if (any(is_MT))
            makeAndAssignProvIds(sgidx[is_MT])
        ## Group (h).
        if (any(is_nbxxx)) {
            gsuffix <- sgsuffix[is_nbxxx]
            ints1 <- as.integer(getIntPrefix(gsuffix))
            ints2 <- as.integer(factor(getIntPrefixTail(gsuffix)))
            ints <- (max(ints2) + 1L) * ints1 + ints2
            makeAndAssignProvIds(sgidx[is_nbxxx], ints=ints)
        }
        ## Group (i).
        if (any(is_Xxxx)) {
            gsuffix <- sgsuffix[is_Xxxx]
            ints <- as.integer(factor(dropPrefix(gsuffix, 1L)))
            makeAndAssignProvIds(sgidx[is_Xxxx], ints=ints)
        }
        ## Group (j).
        if (any(is_Yxxx)) {
            gsuffix <- sgsuffix[is_Yxxx]
            ints <- as.integer(factor(dropPrefix(gsuffix, 1L)))
            makeAndAssignProvIds(sgidx[is_Yxxx], ints=ints)
        }
        ## Group (k).
        if (any(is_Uxxx)) {
            gsuffix <- sgsuffix[is_Uxxx]
            ints <- as.integer(factor(dropPrefix(gsuffix, 1L)))
            makeAndAssignProvIds(sgidx[is_Uxxx], ints=ints)
        }
        ## Group (l).
        if (any(is_Mxxx)) {
            gsuffix <- sgsuffix[is_Mxxx]
            ints <- as.integer(factor(dropPrefix(gsuffix, 1L)))
            makeAndAssignProvIds(sgidx[is_Mxxx], ints=ints)
        }
        ## Group (m).
        if (any(is_MTxxx)) {
            gsuffix <- sgsuffix[is_MTxxx]
            ints <- as.integer(factor(dropPrefix(gsuffix, 2L)))
            makeAndAssignProvIds(sgidx[is_MTxxx], ints=ints)
        }
        ## Group (n).
        if (any(is_xxx)) {
            gsuffix <- sgsuffix[is_xxx]
            ints <- as.integer(factor(gsuffix))
            makeAndAssignProvIds(sgidx[is_xxx], ints=ints)
        }
    }

    ## Longest prefixes first.
    assignProvIdsForSuperGroup(seqlevels, "CHR")
    assignProvIdsForSuperGroup(seqlevels, "chr")
    assignProvIdsForSuperGroup(seqlevels, "CH")
    assignProvIdsForSuperGroup(seqlevels, "ch")
    assignProvIdsForSuperGroup(seqlevels, "")

    ## Seqlevel ids.
    uids <- integer(length(prov_ids))
    oo <- order(prov_ids)
    uids[oo] <- seq_len(length(uids))

    ## Seqname ids.
    uids[seqnames_sm]
}

