#include <Rdefines.h>


/* cigar_utils.c */

SEXP valid_cigar(SEXP cigar, SEXP ans_type);

SEXP split_cigar(SEXP cigar);

SEXP cigar_op_table(SEXP cigar);

SEXP cigar_to_qwidth(SEXP cigar, SEXP before_hard_clipping);

SEXP cigar_to_width(SEXP cigar);

SEXP cigar_qnarrow(SEXP cigar, SEXP left_qwidth, SEXP right_qwidth);

SEXP cigar_narrow(SEXP cigar, SEXP left_width, SEXP right_width);

SEXP cigar_to_IRanges(SEXP cigar, SEXP drop_D_ranges, SEXP merge_ranges);

SEXP cigar_to_list_of_IRanges_by_alignment(
	SEXP cigar,
	SEXP pos,
	SEXP flag,
	SEXP drop_D_ranges
);

SEXP cigar_to_list_of_IRanges_by_rname(
	SEXP cigar,
	SEXP rname,
	SEXP pos,
	SEXP flag,
	SEXP drop_D_ranges,
	SEXP merge_ranges
);

