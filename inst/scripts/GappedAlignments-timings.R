###

bam_dir <- "~biocdev/data_store/1000genomes/extdata"

f1 <- "NA19239.SLX.maq.SRP000033.2009_09.subset.bam"
f2 <- "NA19240.chrom1.454.ssaha2.SRP000032.2009_10.bam"
f3 <- "NA19240.chrom6.SLX.maq.SRP000032.2009_07.subset.bam"

timingsGappedAlignments <- function(class, bam_filename)
{
    cat("\n")
    cat("================================================================\n")
    cat("Timings (in seconds) for class: ", class, "\n",
        "and file: ", bam_filename, "\n", sep="")
    cat("----------------------------------------------------------------\n")
    f <- file.path(bam_dir, bam_filename)

    ## readGappedAlignments:
    T0 <- system.time(x <-
            readGappedAlignments(f, ans.subtype=class))[["elapsed"]]
    cat("x <- readGappedAlignments(file): ", T0, "\n", sep="")
    cat("object.size(x): ", object.size(x), "\n", sep="")

    ## Accessors:
    TA1 <- system.time(x_rname <- rname(x))[["elapsed"]]
    cat("x_rname <- rname(x): ", TA1, "\n", sep="")
    TA2 <- system.time(x_strand <- strand(x))[["elapsed"]]
    cat("x_strand <- strand(x): ", TA2, "\n", sep="")
    TA3 <- system.time(x_cigar <- cigar(x))[["elapsed"]]
    cat("x_cigar <- cigar(x): ", TA3, "\n", sep="")
    TA4 <- system.time(x_qwidth <- qwidth(x))[["elapsed"]]
    cat("x_qwidth <- qwidth(x): ", TA4, "\n", sep="")
    TA5 <- system.time(x_grglist <- grglist(x))[["elapsed"]]
    cat("x_grglist <- grglist(x): ", TA5, "\n", sep="")
    TA6 <- system.time(x_rglist <- rglist(x))[["elapsed"]]
    cat("x_rglist <- rglist(x): ", TA6, "\n", sep="")
    TA7 <- system.time(x_start <- start(x))[["elapsed"]]
    cat("x_start <- start(x): ", TA7, "\n", sep="")
    TA8 <- system.time(x_end <- end(x))[["elapsed"]]
    cat("x_end <- end(x): ", TA8, "\n", sep="")
    TA9 <- system.time(x_width <- width(x))[["elapsed"]]
    cat("x_width <- width(x): ", TA9, "\n", sep="")
    TA10 <- system.time(x_ngap <- ngap(x))[["elapsed"]]
    cat("x_ngap <- ngap(x): ", TA10, "\n", sep="")

    ## show():
    Tshow <- system.time({sink("/dev/null"); show(x); sink(NULL)})[["elapsed"]]
    cat("show(x): ", Tshow, "\n", sep="")

    ## Subsetting:
    TS1 <- system.time(xs1 <- x[1000:1])[["elapsed"]]
    cat("xs1 <- x[1000:1]: ", TS1, "\n", sep="")

    ## Mixing basic operations:
    TM1 <- system.time(rname(x) <- sub("seq", "chr", rname(x)))[["elapsed"]]
    cat("rname(x) <- sub(\"seq\", \"chr\", rname(x)): ", TM1, "\n", sep="")
    TM2 <- system.time(x[strand(x) == "-"])[["elapsed"]]
    cat("x[strand(x) == \"-\"]: ", TM2, "\n", sep="")

    ## qnarrow():
    TQN1 <- system.time(qnx1 <- qnarrow(x, start=4, end=-6))[["elapsed"]]
    cat("qnx1 <- qnarrow(x, start=4, end=-6): ", TQN1, "\n", sep="")
    TQN2 <- system.time(qnx2 <- qnarrow(x, start=14, end=-16))[["elapsed"]]
    cat("qnx2 <- qnarrow(x, start=14, end=-16): ", TQN2, "\n", sep="")

    ## narrow():
    TN1 <- system.time(nx1 <- narrow(x, start=4, end=-6))[["elapsed"]]
    cat("nx1 <- narrow(x, start=4, end=-6): ", TN1, "\n", sep="")
    TN2 <- system.time(nx2 <- narrow(x, start=14, end=-13))[["elapsed"]]
    cat("nx2 <- narrow(x, start=14, end=-13): ", TN2, "\n", sep="")

    ## coverage():
    TCVG <- system.time(cvg <- coverage(x))[["elapsed"]]
    cat("cvg <- coverage(x): ", TCVG, "\n", sep="")

    ## findOverlaps()/countOverlaps()/subsetByOverlaps():
    TFO1 <- system.time(fo1 <- findOverlaps(x, grg(x)[1]))[["elapsed"]]
    cat("fo1 <- findOverlaps(x, grg(x)[1]): ", TFO1, "\n", sep="")
    TCO1 <- system.time(co1 <- countOverlaps(x, grg(x)[1]))[["elapsed"]]
    cat("co1 <- countOverlaps(x, grg(x)[1]): ", TCO1, "\n", sep="")
    TSO1 <- system.time(so1 <- subsetByOverlaps(x, grg(x)[1]))[["elapsed"]]
    cat("so1 <- subsetByOverlaps(x, grg(x)[1]): ", TSO1, "\n", sep="")
}

for (file in c(f1, f2, f3)) {
    for (class in c("Alignments0", "Alignments1", "Alignments2"))
        timingsGappedAlignments(class, file)
}

