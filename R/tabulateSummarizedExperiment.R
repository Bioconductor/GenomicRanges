## a directory oriented way to compute overlaps from a GRanges style
## annotation and a directory of Bam files


##  The following must be written as a method that can take either a
##  BamFileList or a BamViews and do the right thing.  It should also support
##  a dir (character) of course.  Internally, I might want to just do it by
##  actually calling the countBam method of BamFileList, which raises the
##  question of how much of this I actually need to be doing...

## countBam() is really slow (367 seconds for 2 files vs 140 seconds). so for
## now, I am going to go with the loop below .  My guess is that it has to do
## with indexing of bam files somehow.
## bfs <- open(bfs)
## param <- ScanBamParam(which=annot)
## foo = countBam(bfs, param=param)



## For now inernal function will take a character vector of files and process
## those into a grid of numbers (then put that into a SummarizedExperiment
## obj.)
.tabulateBamFiles <- function(x, rowData, colData,
                              annotType=c("exon_id", "gene_id", "tx_id"),
                             ignore.strand=TRUE, ...){
  ## TODO: better argument checking 
  ## check args
  annotType <- match.arg(annotType)
  if(is(rowData,"GRanges")){
    if(length(unique(values(rowData)[[annotType]]))!=length(rowData))
    stop("You must select an 'annotType' value that matches the 'rowData' used.")
  }

  len <- seqlengths(rowData)
  ## TODO: also make sure that x is a list of open BamFiles?
  #len <- len[names(len) %in% seqnames(x[[1]])]  ## requires that x be BamFileList???
  which <- GRanges(names(len), IRanges(1, len))
  ## count
  cnt <- lapply(x, function(fl, rowData) {
    message(fl)        
    ga <- readGappedAlignments(fl, which=which)
    if(ignore.strand == TRUE){
      strand(ga) <- "*"
    }
    ## list(region = unlist(countGenomicOverlaps(rowData, ga, ...)) )
    list(region = unlist(countOverlaps(rowData, ga)) )
  }, rowData)
  ## If there are sample names, lets use them
  if(!is.null(names(x))){
    names(cnt) <- names(x)
  }else{ ## name the files after the files I guess
    names(cnt) <- x
  }
  Cnt <- sapply(cnt, "[[", "region")

  
  ## if there is an annotType and it's a GRanges object then use that
  if(!is.null(annotType) &&  is(rowData,"GRanges")){ 
    rownames(Cnt) <- values(rowData)[[annotType]] 
  }
  ## if it is a GRangesList, then use the names for the row names
  if(is(rowData,"GRangesList")){ 
    rownames(Cnt) <- names(rowData)
  }  
  ## Then assemble the object
  colData <- as(colData, "DataFrame")
  SummarizedExperiment(assays=SimpleList(counts = Cnt),
                       rowData=rowData, colData=colData)
}




## TODO: move to allGenerics.R
## change default arguments
## change workhorse function
## add ... to the generic etc.
## add internal check for the BamFileList being open
## rename???
## move to Rsamtools?
setGeneric("tabulateBamFiles", signature="x",
    function(x, rowData, colData,annotType,ignore.strand){standardGeneric("tabulateBamFiles")})






## Define methods:
## for a named character vector (names/paths to files)
setMethod("tabulateBamFiles", "character",
    function(x, rowData, colData,annotType,ignore.strand){
      .tabulateBamFiles(x, rowData, colData, annotType=NULL,
                                       ignore.strand=TRUE)}
)

## x <- BamFileList(fls)
## for a BamFileList
setMethod("tabulateBamFiles", "BamFileList",
    function(x, rowData, colData,annotType,ignore.strand){
      ## convert to character list
      x <- unlist(lapply(x,path))
      .tabulateBamFiles(x, rowData, colData, annotType=NULL,
                        ignore.strand=TRUE)
    }
)


## x <- BamViews(fls, bamRanges=rngs)
## for a BamViews
setMethod("tabulateBamFiles", "BamViews",
    function(x, rowData, colData,annotType,ignore.strand){
      ## convert to character list
      x <- bamPaths(x)
      .tabulateBamFiles(x, rowData, colData, annotType=NULL,
                        ignore.strand=TRUE)
    }
)




## TODO: get better smaller example for the manual page and describe all of the arguments properly.


