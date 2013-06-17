#include "GenomicRanges.h"
#include "IRanges_interface.h"

#include <ctype.h> /* for isdigit() */

static char errmsg_buf[200];

/* Return the number of chars that was read, or 0 if there is no more char
   to read (i.e. cig0[offset] is '\0'), or -1 in case of a parse error.
   Zero-length operations are ignored. */
static int next_cigar_OP(const char *cig0, int offset, char *OP, int *OPL)
{
	char c;
	int offset0, opl;

	if (!cig0[offset])
		return 0;
	offset0 = offset;
	do {
		/* Extract *OPL */
		opl = 0;
		while (isdigit(c = cig0[offset])) {
			offset++;
			opl *= 10;
			opl += c - '0';
		}
		/* Extract *OP */
		if (!(*OP = cig0[offset])) {
			snprintf(errmsg_buf, sizeof(errmsg_buf),
				 "unexpected CIGAR end after char %d",
				 offset);
			return -1;
		}
		offset++;
	} while (opl == 0);
	*OPL = opl;
	return offset - offset0;
}

/* Return the number of chars that was read, or 0 if there is no more char
   to read (i.e. offset is 0), or -1 in case of a parse error.
   Zero-length operations are ignored. */
static int prev_cigar_OP(const char *cig0, int offset, char *OP, int *OPL)
{
	char c;
	int offset0, opl, powof10;

	if (offset == 0)
		return 0;
	offset0 = offset;
	do {
		/* Extract *OP */
		offset--;
		*OP = cig0[offset];
		/* Extract *OPL */
		if (offset == 0) {
			snprintf(errmsg_buf, sizeof(errmsg_buf),
				 "no CIGAR operation length before char %d",
				 offset + 1);
			return -1;
		}
		offset--;
		opl = 0;
		powof10 = 1;
		while (offset >= 0 && isdigit(c = cig0[offset])) {
			opl += (c - '0') * powof10;
			powof10 *= 10;
			offset--;
		}
		offset++;
	} while (opl == 0);
	*OPL = opl;
	return offset0 - offset;
}

static const char *split_cigar_string(SEXP cigar_string,
		CharAE *OPbuf, IntAE *OPLbuf)
{
	const char *cig0;
	int offset, n, OPL /* Operation Length */;
	char OP /* Operation */;

	cig0 = CHAR(cigar_string);
	offset = 0;
	while ((n = next_cigar_OP(cig0, offset, &OP, &OPL))) {
		if (n == -1)
			return errmsg_buf;
		if (OPbuf != NULL)
			CharAE_insert_at(OPbuf, CharAE_get_nelt(OPbuf), OP);
		if (OPLbuf != NULL)
			IntAE_insert_at(OPLbuf, IntAE_get_nelt(OPLbuf), OPL);
		offset += n;
	}
	return NULL;
}

#define QUERY_BEFORE_HARD_CLIPPING	0
#define QUERY				1
#define QUERY_AFTER_SOFT_CLIPPING	2
#define PAIRWISE			3
#define REFERENCE			4

static int get_cigar_OP_width(char OP, int OPL, int space)
{
	if (OP == 'M')
		return OPL;
	switch (space) {
	case QUERY_BEFORE_HARD_CLIPPING:
		if (OP == 'H')
			return OPL;
		/* fall through */
	case QUERY:
		if (OP == 'S')
			return OPL;
		/* fall through */
	case QUERY_AFTER_SOFT_CLIPPING:
		if (OP == 'I')
			return OPL;
		break;
	case PAIRWISE:
		if (OP == 'I')
			return OPL;
		/* fall through */
	case REFERENCE:
		if (OP == 'N' || OP == 'D')
			return OPL;
		break;
	}
	if (OP == '=' || OP == 'X')
		return OPL;
	return 0;
}

static int is_in_ops(char OP, SEXP ops)
{
	int ops_len, i;

	if (ops == R_NilValue)
		return 1;
	ops_len = LENGTH(ops);
	for (i = 0; i < ops_len; i++) {
		if (OP == CHAR(STRING_ELT(ops, i))[0])
			return 1;
	}
	return 0;
}

static void drop_or_append_or_merge_range(int start, int width,
		int drop_empty_range, int merge_range, int nelt0,
		RangeAE *range_buf, const char *OP, CharAEAE *OP_buf)
{
	int buf_nelt, buf_nelt_minus_1, prev_end_plus_1;
	CharAE OP_buf_new_elt, *OP_buf_prev_elt;

	if (drop_empty_range && width == 0)  /* Drop. */
		return;
	buf_nelt = RangeAE_get_nelt(range_buf);
	if (merge_range && buf_nelt > nelt0) {
		/* The incoming range should never overlap with the previous
		   incoming range i.e. 'start' should always be > the end of
		   the previous incoming range. */
		buf_nelt_minus_1 = buf_nelt - 1;
		prev_end_plus_1 = range_buf->start.elts[buf_nelt_minus_1] +
				  range_buf->width.elts[buf_nelt_minus_1];
		if (start == prev_end_plus_1) {
			/* Merge. */
			range_buf->width.elts[buf_nelt_minus_1] += width;
			if (OP_buf != NULL) {
				OP_buf_prev_elt = OP_buf->elts +
						  buf_nelt_minus_1;
				CharAE_insert_at(OP_buf_prev_elt,
					CharAE_get_nelt(OP_buf_prev_elt), *OP);
			}
			return;
		}
	}
	/* Append. */
	RangeAE_insert_at(range_buf, buf_nelt, start, width);
	if (OP_buf != NULL) {
		OP_buf_new_elt = new_CharAE(1);
		CharAE_insert_at(&OP_buf_new_elt, 0, *OP);
		CharAEAE_insert_at(OP_buf, buf_nelt, &OP_buf_new_elt);
	}
	return;
}

static const char *parse_cigar_ranges(const char *cigar_string,
		SEXP ops, int space, int pos,
		int drop_empty_ranges, int reduce_ranges,
		RangeAE *range_buf, CharAEAE *OP_buf)
{
	int buf_nelt0, cigar_offset, n, OPL /* Operation Length */,
	    start, width;
	char OP /* Operation */;

	buf_nelt0 = RangeAE_get_nelt(range_buf);
	cigar_offset = 0;
	start = pos;
	while ((n = next_cigar_OP(cigar_string, cigar_offset, &OP, &OPL))) {
		if (n == -1)
			return errmsg_buf;
		width = get_cigar_OP_width(OP, OPL, space);
		if (is_in_ops(OP, ops))
			drop_or_append_or_merge_range(start, width,
						      drop_empty_ranges,
						      reduce_ranges, buf_nelt0,
						      range_buf, &OP, OP_buf);
		start += width;
		cigar_offset += n;
	}
	return NULL;
}

/*
 * Only the M, =, X, I, and D CIGAR operations produce ranges (1 range per
 * operation). The I operation always produces an empty range. The D operation
 * only produces a range if Ds_as_Ns is FALSE.
 */
static const char *cigar_string_to_ranges(SEXP cigar_string, int pos_elt,
		int Ds_as_Ns, int drop_empty_ranges, int reduce_ranges,
		RangeAE *out)
{
	const char *cig0;
	int out_nelt0, offset, n, OPL /* Operation Length */, start, width;
	char OP /* Operation */;

	cig0 = CHAR(cigar_string);
	out_nelt0 = RangeAE_get_nelt(out);
	offset = 0;
	start = pos_elt;
	while ((n = next_cigar_OP(cig0, offset, &OP, &OPL))) {
		if (n == -1)
			return errmsg_buf;
		width = -1;
		switch (OP) {
		/* Alignment match (can be a sequence match or mismatch) */
		    case 'M': case '=': case 'X': width = OPL; break;
		/* Insertion to the reference */
		    case 'I': width = 0; break;
		/* Deletion from the reference */
		    case 'D':
			if (Ds_as_Ns)
				start += OPL;
			else
				width = OPL;
		    break;
		/* Skipped region from the reference */
		    case 'N': start += OPL; break;
		/* Soft/hard clip on the read */
		    case 'S': case 'H': break;
		/* Silent deletion from the padded reference */
		    case 'P': break;
		    default:
			snprintf(errmsg_buf, sizeof(errmsg_buf),
				 "unknown CIGAR operation '%c' at char %d",
				 OP, offset + 1);
			return errmsg_buf;
		}
		if (width != -1) {
			drop_or_append_or_merge_range(start, width,
						      drop_empty_ranges,
						      reduce_ranges, out_nelt0,
						      out, NULL, NULL);
			start += width;
		}
		offset += n;
	}
	return NULL;
}

static SEXP make_CompressedIRangesList(const RangeAE *range_buf,
		const CharAEAE *OP_buf, SEXP partitioning_end)
{
	SEXP ans, ans_unlistData, ans_unlistData_names, ans_partitioning;

	PROTECT(ans_unlistData =
			new_IRanges_from_RangeAE("IRanges", range_buf));
	if (OP_buf != NULL) {
		PROTECT(ans_unlistData_names =
				new_CHARACTER_from_CharAEAE(OP_buf));
		set_IRanges_names(ans_unlistData, ans_unlistData_names);
		UNPROTECT(1);
	}
	PROTECT(ans_partitioning =
			new_PartitioningByEnd("PartitioningByEnd",
					      partitioning_end,
					      NULL));
	PROTECT(ans = new_CompressedList(
				"CompressedIRangesList",
				ans_unlistData, ans_partitioning));
	UNPROTECT(3);
	return ans;
}

static const char *parse_cigar_width(const char *cigar_string, int space,
		int *width)
{
	int cigar_offset, n, OPL /* Operation Length */;
	char OP /* Operation */;

	*width = cigar_offset = 0;
	while ((n = next_cigar_OP(cigar_string, cigar_offset, &OP, &OPL))) {
		if (n == -1)
			return errmsg_buf;
		*width += get_cigar_OP_width(OP, OPL, space);
		cigar_offset += n;
	}
	return NULL;
}

static const char *Lqnarrow_cigar_string(SEXP cigar_string,
		int *Lqwidth, int *Loffset, int *rshift)
{
	const char *cig0;
	int offset, n, OPL /* Operation Length */;
	char OP /* Operation */;

	if (cigar_string == NA_STRING)
		return "CIGAR string is NA";
	if (LENGTH(cigar_string) == 0)
		return "CIGAR string is empty";
	cig0 = CHAR(cigar_string);
	*rshift = offset = 0;
	while ((n = next_cigar_OP(cig0, offset, &OP, &OPL))) {
		if (n == -1)
			return errmsg_buf;
		switch (OP) {
		/* Alignment match (can be a sequence match or mismatch) */
		    case 'M': case '=': case 'X':
			if (*Lqwidth < OPL) {
				*Loffset = offset;
				*rshift += *Lqwidth;
				return NULL;
			}
			*Lqwidth -= OPL;
			*rshift += OPL;
		    break;
		/* Insertion to the reference or soft/hard clip on the read */
		    case 'I': case 'S': case 'H':
			if (*Lqwidth < OPL) {
				*Loffset = offset;
				return NULL;
			}
			*Lqwidth -= OPL;
		    break;
		/* Deletion (or skipped region) from the reference */
		    case 'D': case 'N':
			*rshift += OPL;
		    break;
		/* Silent deletion from the padded reference */
		    case 'P': break;
		    default:
			snprintf(errmsg_buf, sizeof(errmsg_buf),
				 "unknown CIGAR operation '%c' at char %d",
				 OP, offset + 1);
			return errmsg_buf;
		}
		offset += n;
	}
	snprintf(errmsg_buf, sizeof(errmsg_buf),
		 "CIGAR is empty after qnarrowing");
	return errmsg_buf;
}

static const char *Rqnarrow_cigar_string(SEXP cigar_string,
		int *Rqwidth, int *Roffset)
{
	const char *cig0;
	int offset, n, OPL /* Operation Length */;
	char OP /* Operation */;

	if (cigar_string == NA_STRING)
		return "CIGAR string is NA";
	if (LENGTH(cigar_string) == 0)
		return "CIGAR string is empty";
	cig0 = CHAR(cigar_string);
	offset = LENGTH(cigar_string);
	while ((n = prev_cigar_OP(cig0, offset, &OP, &OPL))) {
		if (n == -1)
			return errmsg_buf;
		offset -= n;
		switch (OP) {
		/* M, =, X, I, S, H */
		    case 'M': case '=': case 'X': case 'I': case 'S': case 'H':
			if (*Rqwidth < OPL) {
				*Roffset = offset;
				return NULL;
			}
			*Rqwidth -= OPL;
		    break;
		/* Deletion (or skipped region) from the reference,
		   or silent deletion from the padded reference */
		    case 'D': case 'N': case 'P':
		    break;
		    default:
			snprintf(errmsg_buf, sizeof(errmsg_buf),
				 "unknown CIGAR operation '%c' at char %d",
				 OP, offset + 1);
			return errmsg_buf;
		}
	}
	snprintf(errmsg_buf, sizeof(errmsg_buf),
		 "CIGAR is empty after qnarrowing");
	return errmsg_buf;
}

/* FIXME: 'cigar_buf' is under the risk of a buffer overflow! */
static const char *qnarrow_cigar_string(SEXP cigar_string,
		int Lqwidth, int Rqwidth, char *cigar_buf, int *rshift)
{
	int Loffset, Roffset, buf_offset;
	const char *cig0;
	int offset, n, OPL /* Operation Length */;
	char OP /* Operation */;
	const char *errmsg;

	//Rprintf("qnarrow_cigar_string():\n");
	errmsg = Lqnarrow_cigar_string(cigar_string, &Lqwidth, &Loffset,
				       rshift);
	if (errmsg != NULL)
		return errmsg;
	//Rprintf("  Lqwidth=%d Loffset=%d *rshift=%d\n",
	//	Lqwidth, Loffset, *rshift);
	errmsg = Rqnarrow_cigar_string(cigar_string, &Rqwidth, &Roffset);
	if (errmsg != NULL)
		return errmsg;
	//Rprintf("  Rqwidth=%d Roffset=%d\n", Rqwidth, Roffset);
	if (Roffset < Loffset) {
		snprintf(errmsg_buf, sizeof(errmsg_buf),
			 "CIGAR is empty after qnarrowing");
		return errmsg_buf;
	}
	buf_offset = 0;
	cig0 = CHAR(cigar_string);
	for (offset = Loffset; offset <= Roffset; offset += n) {
		n = next_cigar_OP(cig0, offset, &OP, &OPL);
		if (offset == Loffset)
			OPL -= Lqwidth;
		if (offset == Roffset)
			OPL -= Rqwidth;
		if (OPL <= 0) {
			snprintf(errmsg_buf, sizeof(errmsg_buf),
				 "CIGAR is empty after qnarrowing");
			return errmsg_buf;
		}
		buf_offset += sprintf(cigar_buf + buf_offset,
				      "%d%c", OPL, OP);
	}
	return NULL;
}

static const char *Lnarrow_cigar_string(SEXP cigar_string,
		int *Lwidth, int *Loffset, int *rshift)
{
	const char *cig0;
	int offset, n, OPL /* Operation Length */;
	char OP /* Operation */;

	if (cigar_string == NA_STRING)
		return "CIGAR string is NA";
	if (LENGTH(cigar_string) == 0)
		return "CIGAR string is empty";
	cig0 = CHAR(cigar_string);
	*rshift = offset = 0;
	while ((n = next_cigar_OP(cig0, offset, &OP, &OPL))) {
		if (n == -1)
			return errmsg_buf;
		switch (OP) {
		/* Alignment match (can be a sequence match or mismatch) */
		    case 'M': case '=': case 'X':
			if (*Lwidth < OPL) {
				*Loffset = offset;
				*rshift += *Lwidth;
				return NULL;
			}
			*Lwidth -= OPL;
			*rshift += OPL;
		    break;
		/* Insertion to the reference or soft/hard clip on the read */
		    case 'I': case 'S': case 'H':
		    break;
		/* Deletion (or skipped region) from the reference */
		    case 'D': case 'N':
			if (*Lwidth < OPL)
				*Lwidth = 0;
			else
				*Lwidth -= OPL;
			*rshift += OPL;
		    break;
		/* Silent deletion from the padded reference */
		    case 'P': break;
		    default:
			snprintf(errmsg_buf, sizeof(errmsg_buf),
				 "unknown CIGAR operation '%c' at char %d",
				 OP, offset + 1);
			return errmsg_buf;
		}
		offset += n;
	}
	snprintf(errmsg_buf, sizeof(errmsg_buf),
		 "CIGAR is empty after narrowing");
	return errmsg_buf;
}

static const char *Rnarrow_cigar_string(SEXP cigar_string,
		int *Rwidth, int *Roffset)
{
	const char *cig0;
	int offset, n, OPL /* Operation Length */;
	char OP /* Operation */;

	if (cigar_string == NA_STRING)
		return "CIGAR string is NA";
	if (LENGTH(cigar_string) == 0)
		return "CIGAR string is empty";
	cig0 = CHAR(cigar_string);
	offset = LENGTH(cigar_string);
	while ((n = prev_cigar_OP(cig0, offset, &OP, &OPL))) {
		if (n == -1)
			return errmsg_buf;
		offset -= n;
		switch (OP) {
		/* Alignment match (can be a sequence match or mismatch) */
		    case 'M': case '=': case 'X':
			if (*Rwidth < OPL) {
				*Roffset = offset;
				return NULL;
			}
			*Rwidth -= OPL;
		    break;
		/* Insertion to the reference or soft/hard clip on the read */
		    case 'I': case 'S': case 'H':
		    break;
		/* Deletion (or skipped region) from the reference */
		    case 'D': case 'N':
			if (*Rwidth < OPL)
				*Rwidth = 0;
			else
				*Rwidth -= OPL;
		    break;
		/* Silent deletion from the padded reference */
		    case 'P': break;
		    default:
			snprintf(errmsg_buf, sizeof(errmsg_buf),
				 "unknown CIGAR operation '%c' at char %d",
				 OP, offset + 1);
			return errmsg_buf;
		}
	}
	snprintf(errmsg_buf, sizeof(errmsg_buf),
		 "CIGAR is empty after narrowing");
	return errmsg_buf;
}

/* FIXME: 'cigar_buf' is under the risk of a buffer overflow! */
static const char *narrow_cigar_string(SEXP cigar_string,
		int Lwidth, int Rwidth, char *cigar_buf, int *rshift)
{
	int Loffset, Roffset, buf_offset;
	const char *cig0;
	int offset, n, OPL /* Operation Length */;
	char OP /* Operation */;
	const char *errmsg;

	//Rprintf("narrow_cigar_string():\n");
	errmsg = Lnarrow_cigar_string(cigar_string, &Lwidth, &Loffset,
				      rshift);
	if (errmsg != NULL)
		return errmsg;
	//Rprintf("  Lwidth=%d Loffset=%d *rshift=%d\n",
	//	Lwidth, Loffset, *rshift);
	errmsg = Rnarrow_cigar_string(cigar_string, &Rwidth, &Roffset);
	if (errmsg != NULL)
		return errmsg;
	//Rprintf("  Rwidth=%d Roffset=%d\n", Rwidth, Roffset);
	if (Roffset < Loffset) {
		snprintf(errmsg_buf, sizeof(errmsg_buf),
			 "CIGAR is empty after narrowing");
		return errmsg_buf;
	}
	buf_offset = 0;
	cig0 = CHAR(cigar_string);
	for (offset = Loffset; offset <= Roffset; offset += n) {
		n = next_cigar_OP(cig0, offset, &OP, &OPL);
		if (offset == Loffset)
			OPL -= Lwidth;
		if (offset == Roffset)
			OPL -= Rwidth;
		if (OPL <= 0) {
			snprintf(errmsg_buf, sizeof(errmsg_buf),
				 "CIGAR is empty after narrowing");
			return errmsg_buf;
		}
		buf_offset += sprintf(cigar_buf + buf_offset,
				      "%d%c", OPL, OP);
	}
	return NULL;
}

static const char *cigar_string_op_table(SEXP cigar_string, const char *allOPs,
		int *table_row, int table_nrow)
{
	const char *cig0, *tmp;
	int offset, n, OPL /* Operation Length */;
	char OP /* Operation */;

	if (cigar_string == NA_STRING)
		return "CIGAR string is NA";
	if (LENGTH(cigar_string) == 0)
		return "CIGAR string is empty";
	cig0 = CHAR(cigar_string);
	offset = 0;
	while ((n = next_cigar_OP(cig0, offset, &OP, &OPL))) {
		if (n == -1)
			return errmsg_buf;
		tmp = strchr(allOPs, (int) OP);
		if (tmp == NULL) {
			snprintf(errmsg_buf, sizeof(errmsg_buf),
				 "unknown CIGAR operation '%c' at char %d",
				 OP, offset + 1);
			return errmsg_buf;
		}
		*(table_row + (tmp - allOPs) * table_nrow) += OPL;
		offset += n;
	}
	return NULL;
}


/****************************************************************************
 * --- .Call ENTRY POINT ---
 * Args:
 *   cigar: character vector containing the extended CIGAR string for each
 *          read;
 *   ans_type: a single integer specifying the type of answer to return:
 *     0: 'ans' is a string describing the first validity failure or NULL;
 *     1: 'ans' is logical vector with TRUE values for valid elements
 *        in 'cigar'.
 */
SEXP valid_cigar(SEXP cigar, SEXP ans_type)
{
	SEXP ans, cigar_elt;
	int cigar_len, ans_type0, i, width;
	const char *cigar_string, *errmsg;
	char string_buf[200];

	cigar_len = LENGTH(cigar);
	ans_type0 = INTEGER(ans_type)[0];
	if (ans_type0 == 1)
		PROTECT(ans = NEW_LOGICAL(cigar_len));
	else
		ans = R_NilValue;
	for (i = 0; i < cigar_len; i++) {
		cigar_elt = STRING_ELT(cigar, i);
		if (cigar_elt == NA_STRING) {
			if (ans_type0 == 1)
				LOGICAL(ans)[i] = 1;
			continue;
		}
		cigar_string = CHAR(cigar_elt);
		if (strcmp(cigar_string, "*") == 0) {
			if (ans_type0 == 1)
				LOGICAL(ans)[i] = 1;
			continue;
		}
		/* We use parse_cigar_width() here just for its ability
                   to parse and detect ill-formed CIGAR strings */
		errmsg = parse_cigar_width(cigar_string, 0L, &width);
		if (ans_type0 == 1) {
			LOGICAL(ans)[i] = errmsg == NULL;
			continue;
		}
		if (errmsg != NULL) {
			snprintf(string_buf, sizeof(string_buf),
				 "element %d is invalid (%s)", i + 1, errmsg);
			return mkString(string_buf);
		}
	}
	if (ans_type0 == 1)
		UNPROTECT(1);
	return ans;
}


/****************************************************************************
 * --- .Call ENTRY POINTS ---
 *    - explode_cigar_ops()
 *    - explode_cigar_op_lengths()
 * Arg:
 *   cigar: character vector containing the extended CIGAR strings to
 *          explode.
 * Both functions return a list of the same length as 'cigar' where each
 * list element is a character vector (for explode_cigar_ops()) or an integer
 * vector (for explode_cigar_op_lengths()). The 2 lists have the same shape,
 * that is, same length() and same elementLengths(). The i-th character vector
 * in the list returned by explode_cigar_ops() contains one single-letter
 * string per CIGAR operation in 'cigar[i]'. The i-th integer vector in the
 * list returned by explode_cigar_op_lengths() contains the corresponding
 * CIGAR operation lengths. Zero-length operations are ignored.
 */
SEXP explode_cigar_ops(SEXP cigar)
{
	SEXP ans, cigar_string, ans_elt, ans_elt_elt;
	int cigar_length, ans_elt_len, i, j;
	CharAE OPbuf;
	const char *errmsg;

	cigar_length = LENGTH(cigar);
	PROTECT(ans = NEW_LIST(cigar_length));
	OPbuf = new_CharAE(0);
	for (i = 0; i < cigar_length; i++) {
		cigar_string = STRING_ELT(cigar, i);
		if (cigar_string == NA_STRING) {
			UNPROTECT(1);
			error("'cigar' contains NAs");
		}
		CharAE_set_nelt(&OPbuf, 0);
		errmsg = split_cigar_string(cigar_string, &OPbuf, NULL);
		if (errmsg != NULL) {
			UNPROTECT(1);
			error("in 'cigar' element %d: %s", i + 1, errmsg);
		}
		ans_elt_len = CharAE_get_nelt(&OPbuf);
		PROTECT(ans_elt = NEW_CHARACTER(ans_elt_len));
		for (j = 0; j < ans_elt_len; j++) {
			PROTECT(ans_elt_elt = mkCharLen(OPbuf.elts + j, 1));
			SET_STRING_ELT(ans_elt, j, ans_elt_elt);
			UNPROTECT(1);
		}
		SET_VECTOR_ELT(ans, i, ans_elt);
		UNPROTECT(1);
	}
	UNPROTECT(1);
	return ans;
}

SEXP explode_cigar_op_lengths(SEXP cigar)
{
	SEXP ans, cigar_string, ans_elt;
	int cigar_length, i;
	IntAE OPLbuf;
	const char *errmsg;

	cigar_length = LENGTH(cigar);
	PROTECT(ans = NEW_LIST(cigar_length));
	OPLbuf = new_IntAE(0, 0, 0);
	for (i = 0; i < cigar_length; i++) {
		cigar_string = STRING_ELT(cigar, i);
		if (cigar_string == NA_STRING) {
			UNPROTECT(1);
			error("'cigar' contains NAs");
		}
		IntAE_set_nelt(&OPLbuf, 0);
		errmsg = split_cigar_string(cigar_string, NULL, &OPLbuf);
		if (errmsg != NULL) {
			UNPROTECT(1);
			error("in 'cigar' element %d: %s", i + 1, errmsg);
		}
		PROTECT(ans_elt = new_INTEGER_from_IntAE(&OPLbuf));
		SET_VECTOR_ELT(ans, i, ans_elt);
		UNPROTECT(1);
	}
	UNPROTECT(1);
	return ans;
}


/****************************************************************************
 * --- .Call ENTRY POINT ---
 * Arg:
 *   cigar: character vector containing the extended CIGAR strings to split.
 * Return a list of the same length as 'cigar' where each element is itself
 * a list with 2 elements of the same lengths, the 1st one being a raw
 * vector containing the CIGAR operations and the 2nd one being an integer
 * vector containing the lengths of the CIGAR operations.
 */
SEXP split_cigar(SEXP cigar)
{
	SEXP ans, cigar_string, ans_elt, ans_elt_elt0, ans_elt_elt1;
	int cigar_length, i;
	CharAE OPbuf;
	IntAE OPLbuf;
	const char *errmsg;

	cigar_length = LENGTH(cigar);
	PROTECT(ans = NEW_LIST(cigar_length));
	OPbuf = new_CharAE(0);
	OPLbuf = new_IntAE(0, 0, 0);
	for (i = 0; i < cigar_length; i++) {
		cigar_string = STRING_ELT(cigar, i);
		if (cigar_string == NA_STRING) {
			UNPROTECT(1);
			error("'cigar' contains NAs");
		}
		CharAE_set_nelt(&OPbuf, 0);
		IntAE_set_nelt(&OPLbuf, 0);
		errmsg = split_cigar_string(cigar_string, &OPbuf, &OPLbuf);
		if (errmsg != NULL) {
			UNPROTECT(1);
			error("in 'cigar' element %d: %s", i + 1, errmsg);
		}
		PROTECT(ans_elt = NEW_LIST(2));
		PROTECT(ans_elt_elt0 = new_RAW_from_CharAE(&OPbuf));
		PROTECT(ans_elt_elt1 = new_INTEGER_from_IntAE(&OPLbuf));
		SET_VECTOR_ELT(ans_elt, 0, ans_elt_elt0);
		SET_VECTOR_ELT(ans_elt, 1, ans_elt_elt1);
		SET_VECTOR_ELT(ans, i, ans_elt);
		UNPROTECT(3);
	}
	UNPROTECT(1);
	return ans;
}


/****************************************************************************
 * --- .Call ENTRY POINT ---
 * Args:
 *   cigar: character vector containing extended CIGAR strings.
 *   flag:  NULL or an integer vector containing the SAM flag for each read.
 *          Serves only as a way to indicate whether a read is mapped or not.
 *          According to the SAM Spec v1.4, flag bit 0x4 is the only reliable
 *          place to tell whether a segment (or read) is mapped (bit is 0)
 *          or not (bit is 1).
 *   ops:   NULL or character vector containing the CIGAR operations to
 *          translate to ranges. If NULL, then all CIGAR operations are
 *          translated.
 *   space: single integer (0: reference, 1:query, 2:pairwise)
 *   pos:   integer vector of length 1 or same length as 'cigar' containing
 *          the 1-based leftmost position/coordinate of the clipped read
 *          sequences.
 *   drop_empty_ranges: TRUE or FALSE.
 *   reduce_ranges: TRUE or FALSE.
 *   with_ops: TRUE or FALSE indicating whether the returned ranges should be
 *          named with their corresponding CIGAR operation.
 * Both functions return a CompressedIRangesList object of the same length as
 * 'cigar'.
 */
SEXP cigar_ranges(SEXP cigar, SEXP flag, SEXP ops, SEXP space, SEXP pos,
		SEXP drop_empty_ranges, SEXP reduce_ranges, SEXP with_ops)
{
	int cigar_len, space0, pos_len,
	    drop_empty_ranges0, reduce_ranges0, with_ops0,
	    i, *end_elt, flag_elt;
	RangeAE range_buf;
	CharAEAE OP_buf, *OP_buf_p;
	SEXP ans, ans_partitioning_end, cigar_elt;
	const int *pos_elt;
	const char *cigar_string, *errmsg;

	cigar_len = LENGTH(cigar);
	space0 = INTEGER(space)[0];
	pos_len = LENGTH(pos);
	pos_elt = INTEGER(pos);
	drop_empty_ranges0 = LOGICAL(drop_empty_ranges)[0];
	reduce_ranges0 = LOGICAL(reduce_ranges)[0];
	with_ops0 = LOGICAL(with_ops)[0];
	/* We will generate at least 'cigar_len' ranges. */
	range_buf = new_RangeAE(cigar_len, 0);
	if (with_ops0) {
		OP_buf = new_CharAEAE(cigar_len, 0);
		OP_buf_p = &OP_buf;
	} else {
		OP_buf_p = NULL;
	}
	PROTECT(ans_partitioning_end = NEW_INTEGER(cigar_len));
	for (i = 0, end_elt = INTEGER(ans_partitioning_end);
	     i < cigar_len;
	     i++, *(end_elt++) = RangeAE_get_nelt(&range_buf))
	{
		if (flag != R_NilValue) {
			flag_elt = INTEGER(flag)[i];
			if (flag_elt == NA_INTEGER) {
				UNPROTECT(1);
				error("'flag' contains NAs");
			}
			if (flag_elt & 0x004) {
				/* An unmapped read translates to an empty
				   list element (i.e. IRanges of length 0)
				   in the returned IRangesList object. */
				continue;
			}
		}
		cigar_elt = STRING_ELT(cigar, i);
		if (cigar_elt == NA_STRING) {
			UNPROTECT(1);
			error("'cigar[%d]' is NA", i + 1);
		}
		cigar_string = CHAR(cigar_elt);
		if (strcmp(cigar_string, "*") == 0) {
			UNPROTECT(1);
			error("'cigar[%d]' is \"*\"", i + 1);
		}
		if (*pos_elt == NA_INTEGER || *pos_elt == 0) {
			UNPROTECT(1);
			error("'pos[%d]' is NA or 0", i + 1);
		}
		errmsg = parse_cigar_ranges(cigar_string, ops, space0,
				*pos_elt,
				drop_empty_ranges0, reduce_ranges0,
				&range_buf, OP_buf_p);
		if (errmsg != NULL) {
			UNPROTECT(1);
			error("in 'cigar[%d]': %s", i + 1, errmsg);
		}
		if (pos_len != 1)
			pos_elt++;
	}
	PROTECT(ans = make_CompressedIRangesList(&range_buf, OP_buf_p,
						 ans_partitioning_end));
	UNPROTECT(2);
	return ans;
}


/****************************************************************************
 * --- .Call ENTRY POINT ---
 * Args:
 *   cigar: character vector containing the extended CIGAR string for each
 *          read;
 *   pos:   integer vector containing the 1-based leftmost position/coordinate
 *          of the clipped read sequence;
 *   flag:  NULL or an integer vector containing the SAM flag for each
 *          read;
 *   drop_D_ranges: TRUE or FALSE indicating whether Ds should be treated
 *          like Ns or not;
 * 'cigar', 'pos' and 'flag' (when not NULL) are assumed to have the same
 * length (which is the number of aligned reads).
 *
 * Returns a CompressedIRangesList object of the same length as the input.
 * NOTE: See note for cigar_to_list_of_IRanges_by_rname() below about the
 * strand.
 * TODO: Support character factor 'cigar' in addition to current character
 *       vector format.
 */
SEXP cigar_to_list_of_IRanges_by_alignment(SEXP cigar, SEXP pos, SEXP flag,
		SEXP drop_D_ranges, SEXP drop_empty_ranges, SEXP reduce_ranges)
{
	int cigar_len, pos_len, Ds_as_Ns, drop_empty_ranges0, reduce_ranges0,
	    i, *end_elt, flag_elt;
	RangeAE range_buf;
	SEXP ans, ans_partitioning_end, cigar_elt;
	const int *pos_elt;
	const char *cigar_string, *errmsg;

	cigar_len = LENGTH(cigar);
	pos_len = LENGTH(pos);
	pos_elt = INTEGER(pos);
	Ds_as_Ns = LOGICAL(drop_D_ranges)[0];
	drop_empty_ranges0 = LOGICAL(drop_empty_ranges)[0];
	reduce_ranges0 = LOGICAL(reduce_ranges)[0];
	/* We will generate at least 'cigar_len' ranges. */
	range_buf = new_RangeAE(cigar_len, 0);
	PROTECT(ans_partitioning_end = NEW_INTEGER(cigar_len));
	for (i = 0, end_elt = INTEGER(ans_partitioning_end);
	     i < cigar_len;
	     i++, *(end_elt++) = RangeAE_get_nelt(&range_buf))
	{
		if (flag != R_NilValue) {
			flag_elt = INTEGER(flag)[i];
			if (flag_elt == NA_INTEGER) {
				UNPROTECT(1);
				error("'flag' contains NAs");
			}
			if (flag_elt & 0x004) {
				/* An unmapped read translates to an empty
				   list element (i.e. IRanges of length 0)
				   in the returned IRangesList object. */
				continue;
			}
		}
		cigar_elt = STRING_ELT(cigar, i);
		if (cigar_elt == NA_STRING) {
			UNPROTECT(1);
			error("'cigar[%d]' is NA", i + 1);
		}
		cigar_string = CHAR(cigar_elt);
		if (strcmp(cigar_string, "*") == 0) {
			UNPROTECT(1);
			error("'cigar[%d]' is \"*\"", i + 1);
		}
		if (*pos_elt == NA_INTEGER || *pos_elt == 0) {
			UNPROTECT(1);
			error("'pos[%d]' is NA or 0", i + 1);
		}
		errmsg = cigar_string_to_ranges(cigar_elt, *pos_elt,
				Ds_as_Ns, drop_empty_ranges0, reduce_ranges0,
				&range_buf);
		if (errmsg != NULL) {
			UNPROTECT(1);
			error("in 'cigar[%d]': %s", i + 1, errmsg);
		}
		if (pos_len != 1)
			pos_elt++;
	}
	PROTECT(ans = make_CompressedIRangesList(&range_buf, NULL,
						 ans_partitioning_end));
	UNPROTECT(2);
	return ans;
}


/****************************************************************************
 * --- .Call ENTRY POINT ---
 * Args:
 *   cigar: character vector containing the extended CIGAR string for each
 *          read;
 *   rname: character factor containing the name of the reference sequence
 *          associated with each read (i.e. the name of the sequence the
 *          read has been aligned to);
 *   pos:   integer vector containing the 1-based leftmost position/coordinate
 *          of the clipped read sequence;
 *   flag:  NULL or an integer vector containing the SAM flag for each
 *          read;
 *   drop_D_ranges: TRUE or FALSE indicating whether Ds should be treated
 *          like Ns or not;
 *   reduce_ranges: TRUE or FALSE indicating whether adjacent ranges coming
 *          from the same cigar should be merged or not.
 * 'cigar', 'rname', 'pos' and 'flag' (when not NULL) are assumed to have
 * the same length (which is the number of aligned reads).
 *
 * Return a list of IRanges objects named with the factor levels in 'rname'.
 *
 * NOTE: According to the SAM Format Specification (0.1.2-draft 20090820),
 *   the CIGAR (and the read sequence) stored in the SAM file are represented
 *   on the + strand of the reference sequence. This means that, for a read
 *   that aligns to the - strand, the bases have been reverse complemented
 *   from the unmapped read sequence, and that the corresponding CIGAR has
 *   been reversed. So it seems that, for now, we don't need to deal with the
 *   strand information at all (as long as we are only interested in
 *   returning a list of IRanges objects that is suitable for coverage
 *   extraction).
 *
 * TODO:
 * - Support 'rname' of length 1.
 * - Support character factor 'cigar' in addition to current character vector
 *   format.
 */
SEXP cigar_to_list_of_IRanges_by_rname(SEXP cigar, SEXP rname, SEXP pos,
		SEXP flag,
		SEXP drop_D_ranges, SEXP drop_empty_ranges, SEXP reduce_ranges)
{
	SEXP rname_levels, cigar_string, ans, ans_names;
	int ans_length, nreads, Ds_as_Ns, drop_empty_ranges0, reduce_ranges0,
	    i, level, flag_elt;
	RangeAEAE range_aeae;
	const int *pos_elt;
	const char *errmsg;

	rname_levels = GET_LEVELS(rname);
	ans_length = LENGTH(rname_levels);
	range_aeae = new_RangeAEAE(ans_length, ans_length);
	nreads = LENGTH(pos);
	pos_elt = INTEGER(pos);
	Ds_as_Ns = LOGICAL(drop_D_ranges)[0];
	drop_empty_ranges0 = LOGICAL(drop_empty_ranges)[0];
	reduce_ranges0 = LOGICAL(reduce_ranges)[0];
	for (i = 0; i < nreads; i++) {
		if (flag != R_NilValue) {
			flag_elt = INTEGER(flag)[i];
			if (flag_elt == NA_INTEGER)
				error("'flag' contains NAs");
			if (flag_elt & 0x004)
				continue;
		}
		cigar_string = STRING_ELT(cigar, i);
		if (cigar_string == NA_STRING)
			error("'cigar' contains %sNAs",
			      flag != R_NilValue ? "unexpected " : "");
		level = INTEGER(rname)[i];
		if (level == NA_INTEGER)
			error("'rname' contains %sNAs",
			      flag != R_NilValue ? "unexpected " : "");
		if (*pos_elt == NA_INTEGER)
			error("'pos' contains %sNAs",
			      flag != R_NilValue ? "unexpected " : "");
		errmsg = cigar_string_to_ranges(cigar_string, *pos_elt,
				Ds_as_Ns, drop_empty_ranges0, reduce_ranges0,
				range_aeae.elts + level - 1);
		if (errmsg != NULL)
			error("in 'cigar' element %d: %s", i + 1, errmsg);
		if (nreads != 1)
			pos_elt++;
	}
	PROTECT(ans = new_list_of_IRanges_from_RangeAEAE(
				"IRanges", &range_aeae));
	PROTECT(ans_names = duplicate(rname_levels));
	SET_NAMES(ans, ans_names);
	UNPROTECT(2);
	return ans;
}


/****************************************************************************
 * --- .Call ENTRY POINT ---
 * Args:
 *   cigar: character vector containing the extended CIGAR string for each
 *          read.
 *   space: single integer (0: reference, 1:query, 2:pairwise).
 * Return an integer vector of the same length as 'cigar' containing the
 * widths of the alignments as inferred from the cigar information.
 */
SEXP cigar_width(SEXP cigar, SEXP space)
{
	SEXP ans, cigar_elt;
	int cigar_len, space0, i, width;
	const char *cigar_string, *errmsg;

	cigar_len = LENGTH(cigar);
	space0 = INTEGER(space)[0];
	PROTECT(ans = NEW_INTEGER(cigar_len));
	for (i = 0; i < cigar_len; i++) {
		cigar_elt = STRING_ELT(cigar, i);
		if (cigar_elt == NA_STRING) {
			INTEGER(ans)[i] = NA_INTEGER;
			continue;
		}
		cigar_string = CHAR(cigar_elt);
		if (strcmp(cigar_string, "*") == 0) {
			INTEGER(ans)[i] = NA_INTEGER;
			continue;
		}
		errmsg = parse_cigar_width(cigar_string, space0, &width);
		if (errmsg != NULL) {
			UNPROTECT(1);
			error("in 'cigar[%d]': %s", i + 1, errmsg);
		}
		INTEGER(ans)[i] = width;
	}
	UNPROTECT(1);
	return ans;
}


/****************************************************************************
 * --- .Call ENTRY POINT ---
 */
SEXP cigar_narrow(SEXP cigar, SEXP left_width, SEXP right_width)
{
	SEXP ans, ans_cigar, ans_cigar_string, ans_rshift, cigar_string;
	int cigar_length, i;
	static char cigar_buf[1024];
	const char *errmsg;

	cigar_length = LENGTH(cigar);
	PROTECT(ans_cigar = NEW_CHARACTER(cigar_length));
	PROTECT(ans_rshift = NEW_INTEGER(cigar_length));
	for (i = 0; i < cigar_length; i++) {
		cigar_string = STRING_ELT(cigar, i);
		if (cigar_string == NA_STRING) {
			SET_STRING_ELT(ans_cigar, i, NA_STRING);
			INTEGER(ans_rshift)[i] = NA_INTEGER;
			continue;
		}
		errmsg = narrow_cigar_string(cigar_string,
				INTEGER(left_width)[i],
				INTEGER(right_width)[i],
				cigar_buf, INTEGER(ans_rshift) + i);
		if (errmsg != NULL) {
			UNPROTECT(2);
			error("in 'cigar' element %d: %s", i + 1, errmsg);
		}
		PROTECT(ans_cigar_string = mkChar(cigar_buf));
		SET_STRING_ELT(ans_cigar, i, ans_cigar_string);
		UNPROTECT(1);
	}

	PROTECT(ans = NEW_LIST(2));
	SET_VECTOR_ELT(ans, 0, ans_cigar);
	SET_VECTOR_ELT(ans, 1, ans_rshift);
	UNPROTECT(3);
	return ans;
}


/****************************************************************************
 * --- .Call ENTRY POINT ---
 * Return a list of 2 elements: 1st elt is the narrowed cigar vector, 2nd elt
 * is the 'rshift' vector i.e. the integer vector of the same length as 'cigar'
 * that would need to be added to the 'pos' field of a SAM/BAM file as a
 * consequence of this narrowing.
 */
SEXP cigar_qnarrow(SEXP cigar, SEXP left_qwidth, SEXP right_qwidth)
{
	SEXP ans, ans_cigar, ans_cigar_string, ans_rshift, cigar_string;
	int cigar_length, i;
	static char cigar_buf[1024];
	const char *errmsg;

	cigar_length = LENGTH(cigar);
	PROTECT(ans_cigar = NEW_CHARACTER(cigar_length));
	PROTECT(ans_rshift = NEW_INTEGER(cigar_length));
	for (i = 0; i < cigar_length; i++) {
		cigar_string = STRING_ELT(cigar, i);
		if (cigar_string == NA_STRING) {
			SET_STRING_ELT(ans_cigar, i, NA_STRING);
			INTEGER(ans_rshift)[i] = NA_INTEGER;
			continue;
		}
		errmsg = qnarrow_cigar_string(cigar_string,
				INTEGER(left_qwidth)[i],
				INTEGER(right_qwidth)[i],
				cigar_buf, INTEGER(ans_rshift) + i);
		if (errmsg != NULL) {
			UNPROTECT(2);
			error("in 'cigar' element %d: %s", i + 1, errmsg);
		}
		PROTECT(ans_cigar_string = mkChar(cigar_buf));
		SET_STRING_ELT(ans_cigar, i, ans_cigar_string);
		UNPROTECT(1);
	}

	PROTECT(ans = NEW_LIST(2));
	SET_VECTOR_ELT(ans, 0, ans_cigar);
	SET_VECTOR_ELT(ans, 1, ans_rshift);
	UNPROTECT(3);
	return ans;
}


/****************************************************************************
 * --- .Call ENTRY POINT ---
 * Args:
 *   cigar: character vector containing the extended CIGAR string for each
 *          read;
 * Return an integer matrix with the number of rows equal to the length of
 * 'cigar' and 7 columns, one for each extended CIGAR operation containing
 * a frequency count for the operations for each element of 'cigar'.
 */
SEXP cigar_op_table(SEXP cigar)
{
	SEXP cigar_string, ans, ans_dimnames, ans_colnames;
	int cigar_length, allOPs_length, i, j, *ans_row;
	const char *allOPs = "MIDNSHP", *errmsg;
	char OPstrbuf[2];

	cigar_length = LENGTH(cigar);
	allOPs_length = strlen(allOPs);
	PROTECT(ans = allocMatrix(INTSXP, cigar_length, allOPs_length));
	memset(INTEGER(ans), 0, LENGTH(ans) * sizeof(int));
	ans_row = INTEGER(ans);
	for (i = 0, ans_row = INTEGER(ans); i < cigar_length; i++, ans_row++) {
		cigar_string = STRING_ELT(cigar, i);
		if (cigar_string == NA_STRING) {
			INTEGER(ans)[i] = NA_INTEGER;
			continue;
		}
		errmsg = cigar_string_op_table(cigar_string, allOPs,
				ans_row, cigar_length);
		if (errmsg != NULL) {
			UNPROTECT(1);
			error("in 'cigar' element %d: %s", i + 1, errmsg);
		}
	}

	PROTECT(ans_colnames = NEW_CHARACTER(7));
	OPstrbuf[1] = '\0';
	for (j = 0; j < allOPs_length; j++) {
		OPstrbuf[0] = allOPs[j];
		SET_STRING_ELT(ans_colnames, j, mkChar(OPstrbuf));
	}
	PROTECT(ans_dimnames = NEW_LIST(2));
	SET_ELEMENT(ans_dimnames, 0, R_NilValue);
	SET_ELEMENT(ans_dimnames, 1, ans_colnames);
	SET_DIMNAMES(ans, ans_dimnames);
	UNPROTECT(3);
	return ans;
}


/****************************************************************************
 * --- .Call ENTRY POINT ---
 * Args:
 *   ref_locs: global positions in the reference that we will map
 *   cigar: character string containing the extended CIGAR;
 *   pos: reference position at which the query alignment begins
 *        (after clip)
 *   narrow_left: whether to narrow to the left (or right) side of a gap
 * Returns the local query positions. This assumes that the reference
 * positions actually occur in the read alignment region, outside of
 * any deletions or insertions. 
 */
SEXP ref_locs_to_query_locs(SEXP ref_locs, SEXP cigar, SEXP pos,
                            SEXP narrow_left)
{
  int nlocs, i;
  SEXP query_locs;
  Rboolean _narrow_left = asLogical(narrow_left);
  
  nlocs = LENGTH(ref_locs);
  PROTECT(query_locs = allocVector(INTSXP, nlocs));
  
  for (i = 0; i < nlocs; i++) {
    int query_loc = INTEGER(ref_locs)[i] - INTEGER(pos)[i] + 1;
    int n, offset = 0, OPL, query_consumed = 0;
    char OP;
    const char *cig0 = CHAR(STRING_ELT(cigar, i));
    while (query_consumed < query_loc &&
           (n = next_cigar_OP(cig0, offset, &OP, &OPL)))
    {
      switch (OP) {
        /* Alignment match (can be a sequence match or mismatch) */
      case 'M': case '=': case 'X':
        query_consumed += OPL;
        break;
        /* Insertion to the reference */
      case 'I':
        /* Soft clip on the read */
      case 'S':
        query_loc += OPL;
        query_consumed += OPL;
        break;
        /* Deletion from the reference */
      case 'D':
      case 'N': /* Skipped region from reference; narrow to query */
        {
          Rboolean query_loc_past_gap = query_loc - query_consumed > OPL;
          if (query_loc_past_gap) {
            query_loc -= OPL;
          } else {
            if (_narrow_left) {
              query_loc = query_consumed;
            } else {
              query_loc = query_consumed + 1;
            }
          }
        }
        break;
       /* Hard clip on the read */
      case 'H':
        break;
        /* Silent deletion from the padded reference */
      case 'P':
        break;
      default:
        break;
      }
      offset += n;
    }
    if (n == 0)
      error("hit end of cigar string %d: %s", i + 1, cig0);
    INTEGER(query_locs)[i] = query_loc;
  }

  UNPROTECT(1);
  return query_locs;
}


/****************************************************************************
 * --- .Call ENTRY POINT ---
 * Args:
 *   query_locs: local positions in the read that we will map
 *   cigar: character string containing the extended CIGAR;
 *   pos: reference position at which the query alignment begins
 *        (after clip)
 *   narrow_left: whether to narrow to the left (or right) side of a gap
 * Returns the local query positions. This assumes that the reference
 * positions actually occur in the read alignment region, outside of
 * any deletions or insertions. 
 */
SEXP query_locs_to_ref_locs(SEXP query_locs, SEXP cigar, SEXP pos,
                            SEXP narrow_left)
{
  int nlocs, i;
  SEXP ref_locs;
  Rboolean _narrow_left = asLogical(narrow_left);
  
  nlocs = LENGTH(query_locs);
  PROTECT(ref_locs = allocVector(INTSXP, nlocs));
  
  for (i = 0; i < nlocs; i++) {
    int query_loc_i = INTEGER(query_locs)[i];
    int ref_loc = query_loc_i + INTEGER(pos)[i] - 1;
    int n, offset = 0, OPL, query_consumed = 0;
    char OP;
    const char *cig0 = CHAR(STRING_ELT(cigar, i));
    while (query_consumed < query_loc_i &&
           (n = next_cigar_OP(cig0, offset, &OP, &OPL)))
      {
        switch (OP) {
          /* Alignment match (can be a sequence match or mismatch) */
        case 'M': case '=': case 'X':
          query_consumed += OPL;
          break;
          /* Insertion to the reference */
        case 'I': {
          /* Soft clip on the read */
          int width_from_insertion_start = query_loc_i - query_consumed;
          Rboolean query_loc_past_insertion = width_from_insertion_start > OPL;
          if (query_loc_past_insertion) {
            ref_loc -= OPL;
          } else {
            ref_loc -= width_from_insertion_start;
            if (!_narrow_left) {
              ref_loc += 1;
            }
          }
          query_consumed += OPL;
          break;
        }
        case 'S':
          query_consumed += OPL;
          ref_loc -= OPL;
          break;
          /* Deletion from the reference */
        case 'D':
        case 'N': /* Skipped region from reference; narrow to query */
          ref_loc += OPL;
          break;
          /* Hard clip on the read */
        case 'H':
          break;
          /* Silent deletion from the padded reference */
        case 'P':
          break;
        default:
          break;
        }
        offset += n;
      }
    if (n == 0)
      error("hit end of cigar string %d: %s", i + 1, cig0);
    INTEGER(ref_locs)[i] = ref_loc;
  }

  UNPROTECT(1);
  return ref_locs;
}

