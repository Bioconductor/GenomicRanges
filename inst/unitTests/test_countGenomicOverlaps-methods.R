rng1 <- function(s, w)
GRanges(seq="chr1", IRanges(s, width=w), strand="+")
rng2 <- function(s, w)
GRanges(seq="chr2", IRanges(s, width=w), strand="+")
rng3 <- function(s, w)
GRanges(seq="chr2", IRanges(s, width=w), strand="-")

make_query <- function() 
{
    GRangesList(
        A=rng1(1000, 500),
        B=rng2(2000, 900),
        C=rng1(c(3000, 3600), c(500, 300)),
        D=rng2(c(7000, 7500), c(600, 300)),
        E1=rng1(4000, 500), E2=rng1(c(4300, 4500), c(400, 400)),
        F=rng2(3000, 500),
        G=rng1(c(5000, 5600), c(500, 300)),
        H1=rng1(6000, 500), H2=rng1(6600, 400),
        I=rng3(8000, 200))
}

make_GRsubject <- function() 
{
    GRangesList(
        a=rng1(1400, 500),
        b=rng2(2700, 100),
        c=rng1(3400, 300),
        d=rng2(7100, 600),
        e=rng1(4200, 500),
        f=rng2(c(3100, 3300), 50),
        g=rng1(c(5400, 5600), 50),
        h=rng1(c(6400, 6600), 50),
        i=rng2(8000, 100))
}

make_GAsubject <- function() 
{
    GappedAlignments(
        rname=(c("chr1", "chr2", "chr1", "chr2", "chr1", "chr2", 
               "chr1", "chr1", "chr2")),
        pos=as.integer(c(1400, 2700, 3400, 7100, 4200, 3100, 5400, 6400, 8000)),
        cigar=c("500M", "100M", "300M", "600M", "500M", "50M150N50M", 
                "50M150N50M", "50M150N50M", "100M"),
        strand=strand(rep.int("+", 9L))) 
}

.getCounts <- function(res)
{
    as.vector(assays(res)$counts)
}

test_input_type <- function()
{
    query <- make_query()
    GRsubject <- make_GRsubject()
    GAsubject <- make_GAsubject() 
    GRLq_GRs <- countGenomicOverlaps(query, GRsubject)
    GRLq_GAs <- countGenomicOverlaps(query, GAsubject)
    GRq_GRs <- countGenomicOverlaps(unlist(query), GRsubject)
    GRq_GAs <- countGenomicOverlaps(unlist(query), GAsubject)

    checkIdentical(GRLq_GRs, GRLq_GAs) 
    checkIdentical(GRq_GRs, GRq_GAs) 
}

test_typeAny <- function(type="any", ...)
{
    query <- make_query()
    subject <- make_GRsubject()

    .ignoreStrandFalse <- function(ignore.strand=FALSE, ...)
    {
        anyNone <- countGenomicOverlaps(query, subject, type, 
                                        resolution="none", 
                                        ignore.strand=ignore.strand)
        anyDivide <- countGenomicOverlaps(query, subject, type, 
                                          resolution="divide",
                                          ignore.strand=ignore.strand)
        anyUD <- countGenomicOverlaps(query, subject, type, 
                                      resolution="uniqueDisjoint",
                                      ignore.strand=ignore.strand)
 
        checkIdentical(round(.getCounts(anyNone), 2), c(1, 1, integer(7), 1, 
            rep.int(0.5, 4), 0))
        checkIdentical(round(.getCounts(anyDivide), 2), c(1, 1, rep.int(0.5, 4), 
                       rep.int(0.33, 3), 1, rep.int(0.5, 4), 0))
        checkIdentical(round(.getCounts(anyUD), 2), c(1, 1, integer(4), 1, 
                       0, 0, 1, rep.int(0.5, 4), 0))
    }

    .ignoreStrandTrue <- function(ignore.strand=TRUE, ...)
    {
        anyNone <- countGenomicOverlaps(query, subject, type, 
                                        resolution="none",
                                        ignore.strand=ignore.strand)
        anyDivide <- countGenomicOverlaps(query, subject, type, 
                                          resolution="divide", 
                                          ignore.strand=ignore.strand)
        anyUD <- countGenomicOverlaps(query, subject, type, 
                                      resolution="uniqueDisjoint", 
                                      ignore.strand=ignore.strand)
 
        checkIdentical(round(.getCounts(anyNone), 2), c(1, 1, integer(7), 1, 
            rep(0.5, 4), 1))
        checkIdentical(round(.getCounts(anyDivide), 2), c(1, 1, rep.int(0.5, 4), 
                       rep.int(0.33, 3), 1, rep.int(0.5, 4), 1))
        checkIdentical(round(.getCounts(anyUD), 2), c(1, 1, integer(4), 1, 
                       0, 0, 1, rep.int(0.5, 4), 1))
    }
    .ignoreStrandFalse()
    .ignoreStrandTrue()
} 


test_typeWithin <- function(type="within", ignore.strand=FALSE, ...)
{
    query <- make_query()
    subject <- make_GRsubject()
 
    withinNone <- countGenomicOverlaps(query, subject, type=type, 
        resolution="none", ignore.strand=ignore.strand)
    withinDivide <- countGenomicOverlaps(query, subject, type=type, 
        resolution="divide", ignore.strand=ignore.strand)
 
    checkIdentical(round(.getCounts(withinNone), 2), c(0, 1, integer(7),
        1, rep.int(0.5, 4), 0))
    checkIdentical(round(.getCounts(withinDivide), 2), 
                   c(0, 1, integer(7), 1, rep.int(0.5, 4), 0))
} 




