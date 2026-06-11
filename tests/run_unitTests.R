require("GenomicRanges") || stop("unable to load GenomicRanges package")
require("RUnit") || stop("unable to load RUnit package")
GenomicRanges:::.test()
