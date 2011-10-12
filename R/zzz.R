###

.onUnload <- function(libpath)
{
    library.dynam.unload("GenomicRanges", libpath)
}
