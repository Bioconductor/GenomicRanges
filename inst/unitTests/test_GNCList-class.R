###

findOverlaps_GNCList <- GenomicRanges:::findOverlaps_GNCList

### We need some of the helper functions defined for the NCList unit tests
### in IRanges.
source(system.file("unitTests", "test_NCList-class.R", package="IRanges"))

### Redefine the .get_query_overlaps() and .findOverlaps_naive() functions
### we got from sourcing test_NCList-class.R above.
.get_query_overlaps <- function(query, subject, min.score, type_codes)
{
    if (S4Vectors:::decodeRle(strand(query) == "*")) {
        strand(subject) <- "*"
    } else {
        strand(subject)[strand(subject) == "*"] <- strand(query)
    }
    ok1 <- .overlap_score(query, subject) >= min.score
    ok2 <- rangeComparisonCodeToLetter(compare(query, subject)) %in% type_codes
    which(ok1 & ok2)
}

.findOverlaps_naive <- function(query, subject, min.score=1L,
                                type=c("any", "start", "end",
                                       "within", "extend", "equal"),
                                select=c("all", "first", "last", "arbitrary",
                                         "count"),
                                ignore.strand=FALSE)
{
    if (ignore.strand)
        strand(query) <- "*"
    type <- match.arg(type)
    select <- match.arg(select)
    type_codes <- switch(type,
        "any"    = letters[1:13],
        "start"  = c("f", "g", "h"),
        "end"    = c("d", "g", "j"),
        "within" = c("f", "g", "i", "j"),
        "extend" = c("d", "e", "g", "h"),
        "equal"  = "g"
    )
    hits_per_query <- lapply(seq_along(query),
        function(i) .get_query_overlaps(query[i], subject,
                                        min.score, type_codes))
    hits <- .make_Hits_from_q2s(hits_per_query, length(subject))
    selectHits(hits, select=select)
}

test_findOverlaps_NCList <- function()
{
    q_ranges <- IRanges(-3:7, width=3)
    s_ranges <- IRanges(rep.int(1:5, 5:1), c(1:5, 2:5, 3:5, 4:5, 5))

    query <- GRanges(
        Rle(c("chr1", "chr2", "chrM"), rep(length(q_ranges), 3)),
        rep(q_ranges, 3),
        Rle(c("+", "+", "-"), rep(length(q_ranges), 3)))
    subject <- GRanges(
        Rle(c("chr1", "chr2", "chrM"), rep(length(s_ranges), 3)),
        rep(s_ranges, 3),
        Rle(c("+", "-", "*"), rep(length(s_ranges), 3)))

    target0 <- .findOverlaps_naive(query, subject)
    current <- findOverlaps_GNCList(query, GNCList(subject))
    checkTrue(.compare_hits(target0, current))
    current <- findOverlaps_GNCList(GNCList(query), subject)
    checkTrue(.compare_hits(target0, current))
    current <- findOverlaps_GNCList(query, subject)
    checkTrue(.compare_hits(target0, current))

    ## Shuffle query and/or subject elements.
    permute_input <- function(q_perm, s_perm) {
        q_revperm <- integer(length(q_perm))
        q_revperm[q_perm] <- seq_along(q_perm)
        s_revperm <- integer(length(s_perm))
        s_revperm[s_perm] <- seq_along(s_perm)
        target <- remapHits(target0, query.map=q_revperm,
                                     new.queryLength=length(q_perm),
                                     subject.map=s_revperm,
                                     new.subjectLength=length(s_perm))
        current <- findOverlaps_GNCList(query[q_perm], GNCList(subject[s_perm]))
        checkTrue(.compare_hits(target, current))
        current <- findOverlaps_GNCList(GNCList(query[q_perm]), subject[s_perm])
        checkTrue(.compare_hits(target, current))
        current <- findOverlaps_GNCList(query[q_perm], subject[s_perm])
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

