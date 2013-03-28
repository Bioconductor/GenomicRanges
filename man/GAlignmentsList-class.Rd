\name{GAlignmentsList-class}
\docType{class}

% Class
\alias{class:GAlignmentsList}
\alias{GAlignmentsList-class}
\alias{GAlignmentsList}

% Constructors:
\alias{GAlignmentsList}
\alias{makeGAlignmentsListFromFeatureFragments}
\alias{updateObject,GAlignmentsList-method}
\alias{readGAlignmentsList}

% Accessors:
\alias{names,GAlignmentsList-method}
\alias{names<-,GAlignmentsList-method}
\alias{seqnames,GAlignmentsList-method}
\alias{seqnames<-,GAlignmentsList-method}
\alias{rname,GAlignmentsList-method}
\alias{rname<-,GAlignmentsList-method}
\alias{strand,GAlignmentsList-method}
\alias{strand<-,GAlignmentsList-method}
\alias{cigar,GAlignmentsList-method}
\alias{qwidth,GAlignmentsList-method}
\alias{ngap,GAlignmentsList-method}
\alias{elementMetadata,GAlignmentsList-method}
\alias{elementMetadata<-,GAlignmentsList-method}
\alias{seqinfo,GAlignmentsList-method}
\alias{seqinfo<-,GAlignmentsList-method}
\alias{start,GAlignmentsList-method}
\alias{end,GAlignmentsList-method}
\alias{width,GAlignmentsList-method}

% Coercion:
\alias{as.data.frame,GAlignmentsList-method}
\alias{coerce,GAlignmentsList,GRangesList-method}
\alias{coerce,GAlignmentsList,GRanges-method}
\alias{coerce,GAlignmentsList,RangesList-method}
\alias{coerce,GAlignmentsList,Ranges-method}
\alias{grglist,GAlignmentsList-method}
\alias{granges,GAlignmentsList-method}
\alias{rglist,GAlignmentsList-method}
\alias{ranges,GAlignmentsList-method}
\alias{introns,GAlignmentsList-method}

% Subsetting:
\alias{[,GAlignmentsList-method}
\alias{[<-,GAlignmentsList,ANY,ANY,ANY-method}
\alias{[[<-,GAlignmentsList-method}
\alias{[[<-,GAlignmentsList,ANY,ANY-method}

% Combining:
\alias{c,GAlignmentsList-method}

% show:
\alias{show,GAlignmentsList-method}

% Other methods:
\alias{narrow,GAlignmentsList-method}
\alias{qnarrow,GAlignmentsList-method}

\title{GAlignmentsList objects}

\description{
  The GAlignmentsList class is a container for storing a collection of
  GappedAlignments objects.
}

\details{
  A GAlignmentsList object contains a list of GappedAlignments. The
  majority of operations on this page are described in more detail
  on the GappedAlignments man page, see ?\code{GappedAlignments}.
} 

\section{Constructors}{
  \describe{
    \item{}{
      \code{GAlignmentsList(...)}:
      Creates a GAlignmentsList from a list of GappedAlignments objects.
    }
    \item{}{
      \code{readGAlignmentsList(file, format="BAM", use.names=FALSE, ...)}:
      Read a file containing aligned reads as a GAlignmentsList object.
      Note that this function is just a front-end that delegates to the
      \code{\link[Rsamtools]{readBamGAlignmentsList}} function defined
      in the Rsamtools package. 
      See \code{?\link[Rsamtools]{readBamGAlignmentsList}} for
      more information.
 
      The \code{param} and \code{asProperPairs} arguments described on the
      \code{\link[Rsamtools]{readBamGAlignmentsList}} man page fine tune the 
      records are returned. \code{param} is specifed by the standard
      \code{ScanBamParam()} options. When \code{asProperPairs=TRUE} (default), 
      only the proper pairs are returned. This is the same QC filtering that 
      is done by \code{readBamGappedAlignmentPairs} but the records are
      returned in a \code{GAlignmentsList} instead of
      \code{GappedAlignmentsPairs} object. If \code{asProperPairs=FALSE}, all 
      records in the file are retained (i.e., no QC) and the list structure is 
      split by read id (QNAME in SAM/BAM). This format allows the flexibility 
      to retain singletons, pairs, and multi-fragment reads in a single class. 
    }

    \item{}{
      \code{makeGAlignmentsListFromFeatureFragments(seqnames=Rle(factor()), 
                                                    fragmentPos=list(),
                                                    fragmentCigar=list(),
                                                    strand=character(0),
                                                    sep=",")}:
      Constructs a GAlignmentsList from a list of fragmented features.
    }
  }
}

\section{Accessors}{
  In the code snippets below, \code{x} is a GAlignmentsList object.

  \describe{
    \item{}{
      \code{length(x)}:
      Return the number of elements in \code{x}.
    }
    \item{}{
      \code{names(x)}, \code{names(x) <- value}:
      Get or set the names of the elements of \code{x}.
    }
    \item{}{
      \code{seqnames(x)}, \code{seqnames(x) <- value}:
      Get or set the name of the reference sequences of the
      alignments in each element of \code{x}.
    }
    \item{}{
      \code{rname(x)}, \code{rname(x) <- value}:
      Same as \code{seqnames(x)} and \code{seqnames(x) <- value}.
    }
    \item{}{
      \code{strand(x)}, \code{strand(x) <- value}:
      Get or set the strand of the alignments in each element 
      of \code{x}.
    }
    \item{}{
      \code{cigar(x)}:
      Returns a character list of length \code{length(x)}
      containing the CIGAR string for the alignments in
      each element of \code{x}.
    }
    \item{}{
      \code{qwidth(x)}:
      Returns an integer list of length \code{length(x)}
      containing the length of the alignments in each element of
      \code{x} *after* hard clipping (i.e. the length of the 
      query sequence that is stored in the corresponding SAM/BAM record).
    }
    \item{}{
      \code{start(x)}, \code{end(x)}:
      Returns an integer list of length \code{length(x)}
      containing the "start" and "end" (respectively) of the 
      alignments in each element of \code{x}. 
    }
    \item{}{
      \code{width(x)}:
      Returns an integer list of length \code{length(x)} containing
      the "width" of the alignments in each element of \code{x}.
    }
    \item{}{
      \code{ngap(x)}:
      Returns an integer list of length \code{x} containing the number
      of gaps (i.e. N operations in the CIGAR) for the alignments
      in each element of \code{x}.
    }
    \item{}{
      \code{seqinfo(x)}, \code{seqinfo(x) <- value}:
      Get or set the information about the underlying sequences in each
      element of \code{x}. \code{value} must be a list of \link{Seqinfo} 
      objects.
    }
    \item{}{
      \code{seqlevels(x)}, \code{seqlevels(x) <- value}:
      Get or set the sequence levels of the alignments in each element
      of \code{x}.
    }
    \item{}{
      \code{seqlengths(x)}, \code{seqlengths(x) <- value}:
      Get or set the sequence lengths for each element of \code{x}.
      \code{seqlengths(x)} is equivalent to \code{seqlengths(seqinfo(x))}.
      \code{value} can be a named non-negative integer or numeric vector
      eventually with NAs.
    }
    \item{}{
      \code{isCircular(x)}, \code{isCircular(x) <- value}:
      Get or set the circularity flags for the alignments in each
      element in \code{x}. \code{value} must be a named logical list 
      eventually with NAs.
    }
    \item{}{
      \code{genome(x)}, \code{genome(x) <- value}:
      Get or set the genome identifier or assembly name for the alignments 
      in each element of \code{x}. \code{value} must be a named character 
      list eventually with NAs.
    }
    \item{}{
      \code{seqnameStyle(x)}:
      Get or set the seqname style for alignments in each element of \code{x}.
    }
  }
}

\section{Coercion}{
  In the code snippets below, \code{x} is a GAlignmentsList object.

  \describe{
    \item{}{
      \code{granges(x, ignore.strand=FALSE)}, 
      \code{ranges(x)}:
      Return either a \link{GRanges} or a \link[IRanges]{IRanges}
      object of length \code{length(x)}. Note this coercion IGNORES 
      the cigar information. The resulting ranges span the entire
      range, including any gaps or spaces between paired-end reads.

      The \code{granges} coercion supports \code{ignore.strand} to
      allow ranges of opposite strand to be combined (see examples).
    }
    \item{}{
      \code{grglist(x)}, 
      \code{rglist(x)}:
      Return either a \link{GRangesList} or a \link[IRanges]{IRangesList}
      object of length \code{length(x)}. This coercion RESPECTS the cigar 
      information. The resulting ranges are fragments of the original ranges 
      that do not include gaps or spaces between paired-end reads. 
    }
    \item{}{
      \code{as(x, "GRangesList")}, \code{as(x, "GRanges")},
      \code{as(x, "RangesList")}, \code{as(x, "Ranges")}:
      Alternate ways of doing \code{grglist(x)}, \code{granges(x)},
      \code{rglist(x)}, \code{ranges(x)}, respectively.
    }
  }
}

\section{Subsetting and related operations}{
  In the code snippets below, \code{x} is a GAlignmentsList object.

  \describe{
    \item{}{
      \code{x[i]}, \code{x[i] <- value}:
      Get or set list elements \code{i}. \code{i} can be a numeric 
      or logical vector. \code{value} must be a GappedAlignments.
    }
    \item{}{
      \code{x[[i]]}, \code{x[[i]] <- value}:
      Same as \code{x[i]}, \code{x[i] <- value}. 
    }
    \item{}{
      \code{x[i, j]}, \code{x[i, j] <- value}:
      Get or set list elements \code{i} with optional metadata columns
      \code{j}. \code{i} can be a numeric, logical or missing. 
      \code{value} must be a GappedAlignments.
    }
  }
}

\section{Combining}{
  \describe{
    \item{}{
      \code{c(...)}:
      Concatenates the GAlignmentsList objects in \code{...}.
    }
  }
}

\references{
  \url{http://samtools.sourceforge.net/}
}

\author{Valerie Obenchain <vobencha@fhcrc.org}

\seealso{
  \itemize{
    \item \link{GappedAlignments-class}.
    \item \link{GappedAlignmentPairs-class}.
    \item \code{\link[Rsamtools]{readBamGappedAlignments}}.
    \item \code{\link[Rsamtools]{readBamGappedAlignmentPairs}}.
    \item \link{findOverlaps-methods}.
  }
}

\examples{
gal1 <- GappedAlignments(
    seqnames=Rle(factor(c("chr1", "chr2", "chr1", "chr3")),
        c(1, 3, 2, 4)),
    pos=1:10, cigar=paste0(10:1, "M"),
    strand=Rle(strand(c("-", "+", "*", "+", "-")), c(1, 2, 2, 3, 2)),
    names=head(letters, 10), score=1:10)

gal2 <- GappedAlignments(
    seqnames=Rle(factor(c("chr2", "chr4")), c(3, 4)), pos=1:7,
    cigar=c("5M", "3M2N3M2N3M", "5M", "10M", "5M1N4M", "8M2N1M", "5M"),
    strand=Rle(strand(c("-", "+")), c(4, 3)),
    names=tail(letters, 7), score=1:7)

galist <- GAlignmentsList(noGaps=gal1, Gaps=gal2)


## ---------------------------------------------------------------------
## A. BASIC MANIPULATION
## ---------------------------------------------------------------------

length(galist)
names(galist)
seqnames(galist)
strand(galist)
head(cigar(galist))
head(qwidth(galist))
head(start(galist))
head(end(galist))
head(width(galist))
head(ngap(galist))
seqlevels(galist)

## Rename the reference sequences:
seqlevels(galist) <- sub("chr", "seq", seqlevels(galist))
seqlevels(galist)

grglist(galist)  # a GRangesList object
rglist(galist)   # an IRangesList object

## ---------------------------------------------------------------------
## B. SUBSETTING
## ---------------------------------------------------------------------

galist[strand(galist) == "-"]
gaps <- sapply(galist, function(x) any(grepl("N", cigar(x), fixed=TRUE)))
galist[gaps]

## Different ways to subset:
galist[2]             # a GappedAlignments object of length 1
galist[[2]]           # a GappedAlignments object of length 1
grglist(galist[2])  # a GRangesList object of length 1
rglist(galist[2])   # a NormalIRangesList object of length 1

## ---------------------------------------------------------------------
## C. mcols()/elementMetadata()
## ---------------------------------------------------------------------

## Metadata can be defined on the individual GappedAlignment elements
## and the overall GAlignmentsList object. By default, 'level=between' 
## extracts the GALignmentsList metadata. Using 'level=within' 
## will extract the metadata on the individual GappedAliginments.

mcols(galist) ## no metadata on the GAlignmentsList object
mcols(galist, level="within")


## ---------------------------------------------------------------------
## D. readBamGAlignmentsList() 
## ---------------------------------------------------------------------

library(pasillaBamSubset)
library(Rsamtools)
fl <- sortBam(untreated3_chr4(), tempfile(), byQname=TRUE)
bf <- BamFile(fl, index=character(0), yieldSize=50, obeyQname=TRUE)
galist <- readGAlignmentsList(bf, asProperPairs=FALSE)
galist

## ---------------------------------------------------------------------
## E. COERCION 
## ---------------------------------------------------------------------

## The granges() coercion supports 'ignore.strand' to allow ranges
## from different strand to be combined. In this example paired-end
## reads aligned to opposite strands were read into a GAlignmentsList. 
## If the desired operation is to combine these ranges, reguardless of 
## gaps or the space between pairs, 'ignore.strand' must be TRUE.
granges(galist[1], ignore.strand=TRUE)
granges(galist[1], ignore.strand=FALSE)

## grglist() splits ranges by gap and the space bein each list element
galist <- GAlignmentsList(noGaps=gal1, Gaps=gal2)
grglist(galist)
}

\keyword{methods}
\keyword{classes}

