setMethod("findOverlaps", c("GenomicRanges", "GIntervalTree"),
  function(query, subject, maxgap=0L, minoverlap=1L,
         type=c("any","start","end","within","equal"),
         select=c("all","first","last","arbitrary"),
         ignore.strand=FALSE) {

    if (any(isCircular(query), na.rm=TRUE)) 
     stop("'GIntervalTree' not supported for circular sequences")

    if (!isSingleNumber(maxgap) || maxgap < 0L)
      stop("'maxgap' must be a non-negative integer")
    type <- match.arg(type)
    select <- match.arg(select)
    
    ## merge() also checks that 'query' and 'subject' are based on the
    ## same reference genome.
    seqinfo <- merge(seqinfo(query), seqinfo(subject))
    
    hits <- findOverlaps(ranges(query), subject@ranges,
                            maxgap=maxgap,minoverlap=minoverlap,
                            type=type,select="all",
                            partition=seqnames(query))
     
    if (!ignore.strand) {
      q_strand <- .strandAsSignedNumber(strand(query))
      s_strand <- .strandAsSignedNumber(strand(subject))
      
      compatible_strand <- q_strand[queryHits(hits)] *
        s_strand[subjectHits(hits)] != -1L
      
      hits <- hits[compatible_strand]
    }
    
    q_hits <- queryHits(hits)
    s_hits <- subjectHits(hits)
    q_len <- length(query)
    s_len <- length(subject)
    
    if (select == "arbitrary") {
      ans <- rep.int(NA_integer_, q_len)
      ans[q_hits] <- s_hits
      return(ans)
    }
    if (select == "first") {
      ans <- rep.int(NA_integer_, q_len)
      oo <- IRanges:::orderIntegerPairs(q_hits, s_hits, decreasing=TRUE)
      ans[q_hits[oo]] <- s_hits[oo]
      return(ans)
    }
    oo <- IRanges:::orderIntegerPairs(q_hits, s_hits)
    q_hits <- q_hits[oo]
    s_hits <- s_hits[oo]
    if (select == "last") {
      ans <- rep.int(NA_integer_, q_len)
      ans[q_hits] <- s_hits
      return(ans)
    }
       
    new2("Hits", queryHits=q_hits, subjectHits=s_hits,
         queryLength=q_len, subjectLength=s_len,
         check=FALSE)
})

