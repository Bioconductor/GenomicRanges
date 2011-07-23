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
.tabulateBamFiles <- function(x, annot, columnData, annotType=NULL,
                             strandAgnostic=TRUE){
  ## TODO: better argument checking 
  ## check args
  annotType <- match.arg(annotType, c("exon_id","gene_id", "tx_id", NULL))
  if(is(annot,"GRanges")){
    if(length(unique(values(annot)[[annotType]]))!=length(annot))
    stop("You must select an 'annotType' value that matches the 'annot' used.")
  }
  ## count
  cnt <- lapply(x, function(fl, annot) {
    print(fl)        
    ga <- readGappedAlignments(fl)
    if(strandAgnostic == TRUE){
      strand(ga) <- "*"
    }
    idx <- rname(ga) %in% names(seqlengths(annot))
    ga <- ga[idx]
    list(region = unlist(countOverlaps(annot, ga)) )
  }, annot)
  ## If there are sample names, lets use them
  if(!is.null(names(x))){
    names(cnt) <- names(x)
  }else{ ## name the files after the files I guess
    names(cnt) <- x
  }
  Cnt <- sapply(cnt, "[[", "region")

  
  ## if there is an annotType and it's a GRanges object then use that
  if(!is.null(annotType) &&  is(annot,"GRanges")){ 
    rownames(Cnt) <- values(annot)[[annotType]] 
  }
  ## if it is a GRangesList, then use the names for the row names
  if(is(annot,"GRangesList")){ 
    rownames(Cnt) <- names(annot)
  }  
  ## Then assemble the object
  columnData <- as(columnData, "DataFrame")
  SummarizedExperiment(assays=SimpleList(counts = Cnt),
                       rowData=annot, colData=columnData)
}




## TODO: move to allGenerics.R
setGeneric("tabulateBamFiles", signature="x",
    function(x, annot, columnData,annotType,strandAgnostic){standardGeneric("tabulateBamFiles")})

## Define methods:
## for a named character vector (names/paths to files)
setMethod("tabulateBamFiles", "character",
    function(x, annot, columnData,annotType,strandAgnostic){
      .tabulateBamFiles(x, annot, columnData, annotType=NULL,
                                       strandAgnostic=TRUE)}
)

## x <- BamFileList(fls)
## for a BamFileList
setMethod("tabulateBamFiles", "BamFileList",
    function(x, annot, columnData,annotType,strandAgnostic){
      ## convert to character list
      x <- unlist(lapply(x,path))
      .tabulateBamFiles(x, annot, columnData, annotType=NULL,
                        strandAgnostic=TRUE)
    }
)


## x <- BamViews(fls, bamRanges=rngs)
## for a BamViews
setMethod("tabulateBamFiles", "BamViews",
    function(x, annot, columnData,annotType,strandAgnostic){
      ## convert to character list
      x <- bamPaths(x)
      .tabulateBamFiles(x, annot, columnData, annotType=NULL,
                        strandAgnostic=TRUE)
    }
)




## TODO: get better smaller example for the manual page and describe all of the arguments properly.


