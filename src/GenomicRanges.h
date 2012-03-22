#include <Rdefines.h>


/* cigar_utils.c */

SEXP valid_cigar(SEXP cigar, SEXP ans_type);

SEXP split_cigar(SEXP cigar);

SEXP cigar_op_table(SEXP cigar);

SEXP cigar_to_qwidth(SEXP cigar, SEXP before_hard_clipping);

SEXP cigar_to_width(SEXP cigar);

SEXP cigar_qnarrow(SEXP cigar, SEXP left_qwidth, SEXP right_qwidth);

SEXP cigar_narrow(SEXP cigar, SEXP left_width, SEXP right_width);

SEXP cigar_to_IRanges(
	SEXP cigar,
	SEXP drop_D_ranges,
	SEXP drop_empty_ranges,
	SEXP reduce_ranges
);

SEXP cigar_to_list_of_IRanges_by_alignment(
	SEXP cigar,
	SEXP pos,
	SEXP flag,
	SEXP drop_D_ranges,
	SEXP drop_empty_ranges,
	SEXP reduce_ranges
);

SEXP cigar_to_list_of_IRanges_by_rname(
	SEXP cigar,
	SEXP rname,
	SEXP pos,
	SEXP flag,
	SEXP drop_D_ranges,
	SEXP drop_empty_ranges,
	SEXP reduce_ranges
);

SEXP ref_locs_to_query_locs(
        SEXP ref_locs,
        SEXP cigar,
        SEXP pos
);

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

