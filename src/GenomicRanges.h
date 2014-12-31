#include <Rdefines.h>


/* transcript_utils.c */

SEXP transcript_widths(
	SEXP exonStarts,
	SEXP exonEnds
);

SEXP tlocs2rlocs(
	SEXP tlocs,
	SEXP exonStarts,
	SEXP exonEnds,
	SEXP strand,
	SEXP decreasing_rank_on_minus_strand,
	SEXP error_if_out_of_bounds 
);

