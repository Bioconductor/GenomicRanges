##

.TARGET_names <- paste0(letters[1:10], 1L:10L)
.TARGET_seqlevels <- c("chr1", "chr2", "chr3", "chrX")
.TARGET_seqnames <- Rle(factor(c("chr1", "chr3", "chr1", "chrX"),
                               levels = .TARGET_seqlevels),
                        c(2, 4, 2, 2))
.TARGET_start <- 11:20
.TARGET_end <- 12:21
.TARGET_ranges <- IRanges(.TARGET_start, .TARGET_end, names = .TARGET_names)
.TARGET_strand <- Rle(strand(c("-", "+", "*", "+", "-")), c(1, 2, 3, 2, 2))
.TARGET_mcols <- DataFrame(score = 10:19, GC = seq(1, 0, length = 10))
.TARGET_classifier <- c("a", "a", "b", "b", "c", "c", "d", "d", "a", "b")
.TARGET_classifierF <-
    factor(c("a", "a", "b", "b", "c", "c", "d", "d", "a", "b"))

dfr <- data.frame(seqnames = .TARGET_seqnames, .TARGET_ranges,
                 strand = .TARGET_strand, .TARGET_mcols, .TARGET_classifier)
dfac <- data.frame(seqnames = .TARGET_seqnames, .TARGET_ranges,
                 strand = .TARGET_strand, .TARGET_mcols, .TARGET_classifierF)

DF <- DataFrame(seqnames = .TARGET_seqnames, start = .TARGET_start,
                end = .TARGET_end, width = width(.TARGET_ranges),
                names = .TARGET_names, strand = .TARGET_strand, .TARGET_mcols,
                .TARGET_classifier)
DFac <- DataFrame(seqnames = .TARGET_seqnames, start = .TARGET_start,
                end = .TARGET_end, width = width(.TARGET_ranges),
                names = .TARGET_names, strand = .TARGET_strand, .TARGET_mcols,
                .TARGET_classifierF)

.make_TARGET_GRangesList_from_dataframe <- function()
{
        makeGRangesListFromDataFrame(dfr,
                                     split.field = ".TARGET_classifier",
                                     names.field = "names",
                                     keep.extra.columns = TRUE)
}

.make_TARGET_GRangesList_from_DataFrame <- function()
{
        makeGRangesListFromDataFrame(DF,
                                     split.field = ".TARGET_classifier",
                                     names.field = "names",
                                     keep.extra.columns = TRUE)
}

test_makeGRangesListFromDataFrame <- function()
{
    checkException(
        makeGRangesListFromDataFrame(dfr,
                                     split.field = ".TARGET_classifier",
                                     names.field = "Bad_Field"),
        silent = TRUE)
    checkException(
        makeGRangesListFromDataFrame(dfr,
                                     split.field = "Bad_Field"),
        silent = TRUE)
    checkTrue(validObject(
        makeGRangesListFromDataFrame(dfr,
                                     split.field =
                                         ".TARGET_classifier")
    ))
    checkTrue(validObject(
        makeGRangesListFromDataFrame(DF,
                                     split.field =
                                         ".TARGET_classifier")
    ))
    checkTrue(validObject(
        .make_TARGET_GRangesList_from_dataframe()
    ))
    checkTrue(validObject(
        .make_TARGET_GRangesList_from_DataFrame()
    ))
    # Test with factor in data.frame
    checkTrue(validObject(
        makeGRangesListFromDataFrame(
            dfac,
            split.field = ".TARGET_classifierF",
            names.field = "names",
            keep.extra.columns = TRUE)
    ))
    # Test with factor in DataFrame
    checkTrue(validObject(
        makeGRangesListFromDataFrame(
            DFac,
            split.field = ".TARGET_classifierF",
            names.field = "names",
            keep.extra.columns = TRUE)
    ))
    # Test GRangesList length to split.field unique/level length
    checkTrue(identical(length(unique(.TARGET_classifier)),
                        length(.make_TARGET_GRangesList_from_dataframe())))

    checkTrue(identical(length(levels(.TARGET_classifierF)),
                        length(.make_TARGET_GRangesList_from_dataframe())))

    checkTrue(identical(length(unique(.TARGET_classifier)),
                        length(.make_TARGET_GRangesList_from_DataFrame())))

    checkTrue(identical(length(levels(.TARGET_classifierF)),
                        length(.make_TARGET_GRangesList_from_DataFrame())))

    # Names check with and without keep.extra.columns
    checkTrue(identical(
        Reduce(intersect,
               lapply(.make_TARGET_GRangesList_from_dataframe(),
                      function(x)
                          names(mcols(x))
               )),
        c("score", "GC")
    ))
    checkTrue(identical(
        Reduce(intersect,
               lapply(.make_TARGET_GRangesList_from_DataFrame(),
                      function(x)
                          names(mcols(x))
               )),
        c("score", "GC")
    ))
    checkTrue(identical(
        Reduce(intersect,
               lapply(
                   makeGRangesListFromDataFrame(dfr,
                                                split.field =
                                                    ".TARGET_classifier",
                                                names.field = "names"),
                   function(x)
                       names(mcols(x))
               )),
        character(0L)
    ))
    checkTrue(identical(
        Reduce(intersect,
               lapply(
                   makeGRangesListFromDataFrame(DF,
                                                split.field =
                                                    ".TARGET_classifier",
                                                names.field = "names"),
                   function(x)
                       names(mcols(x))
               )),
        character(0L)
    ))
}
