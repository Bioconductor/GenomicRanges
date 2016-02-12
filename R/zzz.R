.onUnload <- function(libpath)
{
    library.dynam.unload("GenomicRanges", libpath)
}

.test <- function() BiocGenerics:::testPackage("GenomicRanges")

