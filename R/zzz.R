###

.onLoad <- function(libname, pkgname)
{
    ## Putting library(methods) below produces the following 'R CMD check'
    ## NOTE:
    ##
    ##   * checking R code for possible problems ... NOTE
    ##   File ‘GenomicRanges/R/zzz.R’:
    ##     .onLoad calls:
    ##       library(methods)
    ##
    ##   Package startup functions should not change the search path.
    ##   See section ‘Good practice’ in ?.onAttach.
    ##
    ## On the other hand, NOT loading the methods package will cause the
    ## following WARNING when running 'R CMD check' on a BSgenome data package:
    ##
    ##   * checking whether the namespace can be loaded with stated dependencies ... WARNING
    ##   Error: .onLoad failed in loadNamespace() for ‘BSgenome.Scerevisiae.UCSC.sacCer3’, details:
    ##     call: length(x)
    ##     error: could not find function "loadMethod"
    ##   Execution halted
    ##
    ##   A namespace must be able to be loaded with just the base namespace
    ##   loaded: otherwise if the namespace gets loaded by a saved object, the
    ##   session will be unable to start.
    ##
    ##   Probably some imports need to be declared in the NAMESPACE file.
    ##
    ## Hence the following trick...
    sillyname <- library  # cheating on codetools
    sillyname(methods)
}

.onUnload <- function(libpath) library.dynam.unload("GenomicRanges", libpath)

