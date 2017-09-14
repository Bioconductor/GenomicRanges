###

findOverlaps_GNCList <- GenomicRanges:::findOverlaps_GNCList

### We need some of the helper functions defined for the NCList unit tests
### in IRanges.
source(system.file("unitTests", "test_NCList-class.R", package="IRanges"))

.get_query_overlaps2 <- function(query, subject,
            maxgap=-1L, minoverlap=0L,
            type=c("any", "start", "end", "within", "extend", "equal"),
            ignore.strand=FALSE)
{
    ok <- .get_query_overlaps(query, subject,
                              maxgap=maxgap, minoverlap=minoverlap,
                              type=type)
    ok <- ok & seqnames(query) == seqnames(subject)
    if (ignore.strand || as.logical(strand(query) == "*"))
        return(ok)
    ok & (strand(subject) == "*" | strand(query) == strand(subject))
}

### Redefine the .findOverlaps_naive() function we got from sourcing
### test_NCList-class.R above.
.findOverlaps_naive <- function(query, subject,
                                maxgap=-1L, minoverlap=0L,
                                type=c("any", "start", "end",
                                       "within", "extend", "equal"),
                                select=c("all", "first", "last", "arbitrary",
                                         "count"),
                                ignore.strand=FALSE)
{
    type <- match.arg(type)
    select <- match.arg(select)
    hits_per_query <- lapply(seq_along(query),
        function(i)
            which(.get_query_overlaps2(query[i], subject,
                                       maxgap=maxgap, minoverlap=minoverlap,
                                       type=type, ignore.strand=ignore.strand)))
    hits <- .make_Hits_from_q2s(hits_per_query, length(subject))
    selectHits(hits, select=select)
}

test_GNCList <- function()
{
    x <- GRanges(Rle(c("chrM", "chr1", "chrM", "chr1"), 4:1),
                 IRanges(1:10, width=5, names=LETTERS[1:10]),
                 strand=rep(c("+", "-"), 5),
                 score=seq(0.7, by=0.045, length.out=10))
    gnclist <- GNCList(x)

    checkTrue(is(gnclist, "GNCList"))
    checkTrue(validObject(gnclist, complete=TRUE))
    checkIdentical(granges(x), granges(gnclist))
    checkIdentical(x, granges(gnclist, use.mcols=TRUE))
    checkIdentical(length(x), length(gnclist))
    checkIdentical(names(x), names(gnclist))
    checkIdentical(seqnames(x), seqnames(gnclist))
    checkIdentical(start(x), start(gnclist))
    checkIdentical(end(x), end(gnclist))
    checkIdentical(width(x), width(gnclist))
    checkIdentical(ranges(x), ranges(gnclist))
    checkIdentical(ranges(x, use.names=FALSE), ranges(gnclist, use.names=FALSE))
    checkIdentical(ranges(x, use.mcols=TRUE), ranges(gnclist, use.mcols=TRUE))
    checkIdentical(strand(x), strand(gnclist))
    checkIdentical(seqinfo(x), seqinfo(gnclist))
    checkIdentical(x, as(gnclist, "GRanges"))
    checkIdentical(x[-6], as(gnclist[-6], "GRanges"))
}

test_findOverlaps_GNCList <- function()
{
    q_ranges <- IRanges(-3:7, width=3)
    s_ranges <- IRanges(rep.int(1:6, 6:1), c(0:5, 1:5, 2:5, 3:5, 4:5, 5))

    query <- GRanges(
        Rle(c("chr1", "chr2", "chrM"), rep(length(q_ranges), 3)),
        rep(q_ranges, 3),
        strand=Rle(c("+", "+", "-"), rep(length(q_ranges), 3)))
    subject <- GRanges(
        Rle(c("chr1", "chr2", "chrM"), rep(length(s_ranges), 3)),
        rep(s_ranges, 3),
        strand=Rle(c("+", "-", "*"), rep(length(s_ranges), 3)))

    for (ignore.strand in c(FALSE, TRUE)) {
        target0 <- .findOverlaps_naive(query, subject,
                                       ignore.strand=ignore.strand)
        current <- findOverlaps_GNCList(query, GNCList(subject),
                                        ignore.strand=ignore.strand)
        checkTrue(.compare_hits(target0, current))
        current <- findOverlaps_GNCList(GNCList(query), subject,
                                        ignore.strand=ignore.strand)
        checkTrue(.compare_hits(target0, current))
        current <- findOverlaps_GNCList(query, subject,
                                        ignore.strand=ignore.strand)
        checkTrue(.compare_hits(target0, current))

        ## Shuffle query and/or subject elements.
        permute_input <- function(q_perm, s_perm) {
            q_revperm <- integer(length(q_perm))
            q_revperm[q_perm] <- seq_along(q_perm)
            s_revperm <- integer(length(s_perm))
            s_revperm[s_perm] <- seq_along(s_perm)
            target <- remapHits(target0, Lnodes.remapping=q_revperm,
                                         new.nLnode=length(q_perm),
                                         Rnodes.remapping=s_revperm,
                                         new.nRnode=length(s_perm))
            current <- findOverlaps_GNCList(query[q_perm],
                                            GNCList(subject[s_perm]),
                                            ignore.strand=ignore.strand)
            checkTrue(.compare_hits(target, current))
            current <- findOverlaps_GNCList(GNCList(query[q_perm]),
                                            subject[s_perm],
                                            ignore.strand=ignore.strand)
            checkTrue(.compare_hits(target, current))
            current <- findOverlaps_GNCList(query[q_perm], subject[s_perm],
                                            ignore.strand=ignore.strand)
            checkTrue(.compare_hits(target, current))
        }

        q_perm <- rev(seq_along(query))
        s_perm <- rev(seq_along(subject))
        permute_input(q_perm, seq_along(subject))  # reverse query
        permute_input(seq_along(query), s_perm)    # reverse subject
        permute_input(q_perm, s_perm)              # reverse both

        set.seed(97)
        for (i in 1:17) {
            ## random permutations
            q_perm <- sample(length(query))
            s_perm <- sample(length(subject))
            permute_input(q_perm, seq_along(subject))
            permute_input(seq_along(query), s_perm)
            permute_input(q_perm, s_perm)
        }
    }
}

