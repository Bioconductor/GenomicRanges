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
	SEXP decreasing_rank_on_minus_strand
);

SEXP extract_transcripts(
	SEXP classname,
	SEXP x,
	SEXP exonStarts,
	SEXP exonEnds,
	SEXP strand,
	SEXP decreasing_rank_on_minus_strand,
	SEXP lkup
);

