#include "GenomicRanges.h"
#include "IRanges_interface.h"

static char errmsg_buf[200];

static int get_nexons(SEXP starts, SEXP ends)
{
	int nstarts, nends;

	if (starts == R_NilValue) {
		nstarts = 0;
	} else if (IS_INTEGER(starts)) {
		nstarts = LENGTH(starts);
	} else {
		snprintf(errmsg_buf, sizeof(errmsg_buf),
			 "'exonStarts' has invalid elements");
		return -1;
	}
	if (ends == R_NilValue) {
		nends = 0;
	} else if (IS_INTEGER(ends)) {
		nends = LENGTH(ends);
	} else {
		snprintf(errmsg_buf, sizeof(errmsg_buf),
			 "'exonEnds' has invalid elements");
		return -1;
	}
	if (nstarts != nends) {
		snprintf(errmsg_buf, sizeof(errmsg_buf),
			 "'exonStarts' and 'exonEnds' have different shapes");
		return -1;
	}
	return nstarts;
}

static int get_transcript_width(SEXP starts, SEXP ends, int ref_length)
{
	int nexons, transcript_width, j, start, end, width;

	if ((nexons = get_nexons(starts, ends)) == -1)
		return -1;
	transcript_width = 0;
	for (j = 0; j < nexons; j++) {
		start = INTEGER(starts)[j];
		end = INTEGER(ends)[j];
		if (start == NA_INTEGER || end == NA_INTEGER) {
			snprintf(errmsg_buf, sizeof(errmsg_buf),
			  "'exonStarts' and/or 'exonEnds' contain NAs'");
			return -1;
		}
		if (ref_length != -1) {
			if (start < 1 || start > ref_length + 1) {
				snprintf(errmsg_buf, sizeof(errmsg_buf),
				  "'exonStarts' contains \"out of limits\" "
				  "values");
				return -1;
			}
			if (end < 0 || end > ref_length) {
				snprintf(errmsg_buf, sizeof(errmsg_buf),
				  "'exonEnds' contains \"out of limits\" "
				  "values");
				return -1;
			}
		}
		width = end - start + 1;
		if (width < 0) {
			snprintf(errmsg_buf, sizeof(errmsg_buf),
			  "'exonStarts/exonEnds' define exons "
			  "of negative length");
			return -1;
		}
		transcript_width += width;
	}
	return transcript_width;
}

static SEXP mk_transcript_widths(SEXP exonStarts, SEXP exonEnds, int ref_length)
{
	int ntranscripts, i, transcript_width;
	SEXP ans, starts, ends;

	ntranscripts = LENGTH(exonStarts);
	PROTECT(ans = NEW_INTEGER(ntranscripts));
	for (i = 0; i < ntranscripts; i++) {
		starts = VECTOR_ELT(exonStarts, i);
		ends = VECTOR_ELT(exonEnds, i);
		transcript_width = get_transcript_width(starts, ends,
					ref_length);
		if (transcript_width == -1) {
			UNPROTECT(1);
			error("%s", errmsg_buf);
		}
		INTEGER(ans)[i] = transcript_width;
	}
	UNPROTECT(1);
	return ans;
}

static int strand_is_minus(SEXP strand, int i)
{
	SEXP strand_elt;

	strand_elt = STRING_ELT(strand, i);
	if (strand_elt != NA_STRING && LENGTH(strand_elt) == 1) {
		switch (CHAR(strand_elt)[0]) {
		    case '+': return 0;
		    case '-': return 1;
		}
	}
	snprintf(errmsg_buf, sizeof(errmsg_buf),
		 "'strand' elements must be \"+\" or \"-\"");
	return -1;
}

static int tloc2rloc(int tloc,
		SEXP starts, SEXP ends,
		int on_minus_strand, int decreasing_rank_on_minus_strand)
{
	int nexons, j, start, end, width;

	nexons = LENGTH(starts);
	if (on_minus_strand && decreasing_rank_on_minus_strand) {
		for (j = nexons - 1; j >= 0; j--) {
			start = INTEGER(starts)[j];
			end = INTEGER(ends)[j];
			width = end - start + 1;
			if (tloc <= width)
				break;
			tloc -= width;
		}
	} else {
		for (j = 0; j < nexons; j++) {
			start = INTEGER(starts)[j];
			end = INTEGER(ends)[j];
			width = end - start + 1;
			if (tloc <= width)
				break;
			tloc -= width;
		}
	}
	tloc--;
	return on_minus_strand ? end - tloc : start + tloc;
}


static int copy_exon(char * out, const cachedCharSeq *in,
		int start, int end, int on_minus_strand, SEXP lkup)
{
	int width;

	start--; end--;
	width = end - start + 1;
	if (on_minus_strand) {
		Orevcopy_bytes_from_i1i2_with_lkup(start, end,
			out, width,
			in->seq, in->length,
			INTEGER(lkup), LENGTH(lkup));
	} else {
		Ocopy_bytes_from_i1i2_with_lkup(start, end,
			out, width,
			in->seq, in->length,
			NULL, 0);
	}
	return width;
}

static void copy_transcript(char * out, const cachedCharSeq *in,
		SEXP starts, SEXP ends,
		int on_minus_strand, int decreasing_rank_on_minus_strand,
		SEXP lkup)
{
	int nexons, j, start, end;

	nexons = LENGTH(starts);
	if (on_minus_strand && decreasing_rank_on_minus_strand) {
		for (j = nexons - 1; j >= 0; j--) {
			start = INTEGER(starts)[j];
			end = INTEGER(ends)[j];
			out += copy_exon(out, in,
				start, end, on_minus_strand, lkup);
		}
	} else {
		for (j = 0; j < nexons; j++) {
			start = INTEGER(starts)[j];
			end = INTEGER(ends)[j];
			out += copy_exon(out, in,
				start, end, on_minus_strand, lkup);
		}
	}
	return;
}


/****************************************************************************
 *                        --- .Call ENTRY POINTS ---                        *
 ****************************************************************************/

SEXP transcript_widths(SEXP exonStarts, SEXP exonEnds)
{
	return mk_transcript_widths(exonStarts, exonEnds, -1);
}

SEXP tlocs2rlocs(SEXP tlocs, SEXP exonStarts, SEXP exonEnds,
		SEXP strand, SEXP decreasing_rank_on_minus_strand)
{
	SEXP ans, starts, ends, ans_elt;
	int decreasing_rank_on_minus_strand0, ans_length,
	    i, transcript_width, on_minus_strand, nlocs, j, tloc;

	decreasing_rank_on_minus_strand0 =
		LOGICAL(decreasing_rank_on_minus_strand)[0];
	ans_length = LENGTH(tlocs);
	PROTECT(ans = duplicate(tlocs));
	for (i = 0; i < ans_length; i++) {
		starts = VECTOR_ELT(exonStarts, i);
		ends = VECTOR_ELT(exonEnds, i);
		transcript_width = get_transcript_width(starts, ends, -1);
		if (transcript_width == -1) {
			UNPROTECT(1);
			error("%s", errmsg_buf);
		}
		on_minus_strand = strand_is_minus(strand, i);
		if (on_minus_strand == -1) {
			UNPROTECT(1);
			error("%s", errmsg_buf);
		}
		ans_elt = VECTOR_ELT(ans, i);
		if (ans_elt == R_NilValue) {
			nlocs = 0;
		} else if (IS_INTEGER(ans_elt)) {
			nlocs = LENGTH(ans_elt);
		} else {
			UNPROTECT(1);
			error("'tlocs' has invalid elements");
		}
		for (j = 0; j < nlocs; j++) {
			tloc = INTEGER(ans_elt)[j];
			if (tloc == NA_INTEGER)
				continue;
			if (tloc < 1 || tloc > transcript_width) {
				UNPROTECT(1);
				error("'tlocs[[%d]]' contains \"out of limits\" "
				      "transcript locations (length of "
				      "transcript is %d)", j + 1, transcript_width);
			}
			INTEGER(ans_elt)[j] = tloc2rloc(tloc,
				starts, ends,
				on_minus_strand,
				decreasing_rank_on_minus_strand0);
		}
	}
	UNPROTECT(1);
	return ans;
}

SEXP extract_transcripts(SEXP classname, SEXP x,
		SEXP exonStarts, SEXP exonEnds, SEXP strand,
		SEXP decreasing_rank_on_minus_strand, SEXP lkup)
{
	cachedCharSeq X, Y;
	SEXP ans_width, ans, starts, ends;
	cachedXVectorList cached_ans;
	int decreasing_rank_on_minus_strand0, ans_length,
	    i, on_minus_strand;

	X = cache_XRaw(x);
	decreasing_rank_on_minus_strand0 =
		LOGICAL(decreasing_rank_on_minus_strand)[0];
	PROTECT(ans_width = mk_transcript_widths(exonStarts,
					exonEnds, X.length));
	PROTECT(ans = alloc_XRawList(CHAR(STRING_ELT(classname, 0)),
					get_classname(x), ans_width));
	cached_ans = cache_XVectorList(ans);
	ans_length = get_cachedXVectorList_length(&cached_ans);
	for (i = 0; i < ans_length; i++) {
		starts = VECTOR_ELT(exonStarts, i);
		if (starts == R_NilValue || LENGTH(starts) == 0)
			continue;
		ends = VECTOR_ELT(exonEnds, i);
		on_minus_strand = strand_is_minus(strand, i);
		if (on_minus_strand == -1) {
			UNPROTECT(2);
			error("%s", errmsg_buf);
		}
		Y = get_cachedXRawList_elt(&cached_ans, i);
		/* Y.seq is a const char * so we need to cast it to
		   char * before we can write to it */
		copy_transcript((char *) Y.seq, &X,
			starts, ends,
			on_minus_strand, decreasing_rank_on_minus_strand0,
			lkup);
	}
	UNPROTECT(2);
	return ans;
}

