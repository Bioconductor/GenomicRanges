### =========================================================================
### countGenomicOverlaps
### -------------------------------------------------------------------------

countGenomicOverlaps <- function(query, subject,
                     type = c("any", "start", "end", "within", "equal"),
                     resolution = c("none", "divide", "uniqueDisjoint"),
                     ignore.strand = FALSE, splitreads = TRUE, ...)
{
    .Defunct("summarizeOverlaps")
}

