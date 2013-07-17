#include "GenomicRanges.h"
#include "IRanges_interface.h"

#include <ctype.h> /* for isdigit() */

/* The 8 supported spaces. */
#define REFERENCE			1
#define REFERENCE_N_REGIONS_REMOVED	2
#define QUERY				3
#define QUERY_BEFORE_HARD_CLIPPING	4
#define QUERY_AFTER_SOFT_CLIPPING	5
#define PAIRWISE			6
#define PAIRWISE_N_REGIONS_REMOVED	7
#define PAIRWISE_DENSE			8

static char errmsg_buf[200];

/* Return the number of chars that was read, or 0 if there is no more char
   to read (i.e. cigar_string[offset] is '\0'), or -1 in case of a parse error.
   Zero-length operations are ignored. */
static int next_cigar_OP(const char *cigar_string, int offset,
		char *OP, int *OPL)
{
	char c;
	int offset0, opl;

	if (!cigar_string[offset])
		return 0;
	offset0 = offset;
	do {
		/* Extract *OPL */
		opl = 0;
		while (isdigit(c = cigar_string[offset])) {
			offset++;
			opl *= 10;
			opl += c - '0';
		}
		/* Extract *OP */
		if (!(*OP = cigar_string[offset])) {
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
static int prev_cigar_OP(const char *cigar_string, int offset,
		char *OP, int *OPL)
{
	char c;
	int offset0, opl, powof10;

	if (offset == 0)
		return 0;
	offset0 = offset;
	do {
		/* Extract *OP */
		offset--;
		*OP = cigar_string[offset];
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
		while (offset >= 0 && isdigit(c = cigar_string[offset])) {
			opl += (c - '0') * powof10;
			powof10 *= 10;
			offset--;
		}
		offset++;
	} while (opl == 0);
	*OPL = opl;
	return offset0 - offset;
}

static int ops_lkup_table[256];

static void init_ops_lkup_table(SEXP ops)
{
	int ops_len, i;
	SEXP ops_elt;
	char OP;

	if (ops == R_NilValue) {
		for (i = 0; i < 256; i++)
			ops_lkup_table[i] = 1;
		return;
	}
	for (i = 0; i < 256; i++)
		ops_lkup_table[i] = 0;
	ops_len = LENGTH(ops);
	for (i = 0; i < ops_len; i++) {
		ops_elt = STRING_ELT(ops, i);
		if (ops_elt == NA_STRING || LENGTH(ops_elt) == 0)
			error("'ops' contains NAs and/or empty strings");
		OP = CHAR(ops_elt)[0];
		ops_lkup_table[(unsigned char) OP] = 1;
	}
	return;
}

static int is_in_ops(char OP)
{
	return ops_lkup_table[(unsigned char) OP];
}

static int is_visible_in_space(char OP, int space)
{
	if (OP == 'M')
		return 1;
	switch (space) {
	case QUERY_BEFORE_HARD_CLIPPING:
		if (OP == 'H')
			return 1;
		/* fall through */
	case QUERY:
		if (OP == 'S')
			return 1;
		/* fall through */
	case QUERY_AFTER_SOFT_CLIPPING:
		if (OP == 'I')
			return 1;
		break;
	case PAIRWISE:
		if (OP == 'I')
			return 1;
		/* fall through */
	case REFERENCE:
		if (OP == 'D' || OP == 'N')
			return 1;
		break;
	case PAIRWISE_N_REGIONS_REMOVED:
		if (OP == 'I')
			return 1;
		/* fall through */
	case REFERENCE_N_REGIONS_REMOVED:
		if (OP == 'D')
			return 1;
	}
	if (OP == '=' || OP == 'X')
		return 1;
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

/* Make sure init_ops_lkup_table() is called before parse_cigar_ranges(). */
static const char *parse_cigar_ranges(const char *cigar_string,
		int space, int pos,
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
		width = is_visible_in_space(OP, space) ? OPL : 0;
		if (is_in_ops(OP))
			drop_or_append_or_merge_range(start, width,
						      drop_empty_ranges,
						      reduce_ranges, buf_nelt0,
						      range_buf, &OP, OP_buf);
		start += width;
		cigar_offset += n;
	}
	return NULL;
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
		if (is_visible_in_space(OP, space))
			*width += OPL;
		cigar_offset += n;
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
 * explode_cigar_ops() and explode_cigar_op_lengths()
 */

/* Make sure init_ops_lkup_table() is called before split_cigar_string(). */
static const char *split_cigar_string(const char *cigar_string,
		CharAE *OPbuf, IntAE *OPLbuf)
{
	int offset, n, OPL /* Operation Length */;
	char OP /* Operation */;

	offset = 0;
	while ((n = next_cigar_OP(cigar_string, offset, &OP, &OPL))) {
		if (n == -1)
			return errmsg_buf;
		if (is_in_ops(OP)) {
			if (OPbuf != NULL)
				CharAE_insert_at(OPbuf,
					CharAE_get_nelt(OPbuf), OP);
			if (OPLbuf != NULL)
				IntAE_insert_at(OPLbuf,
					IntAE_get_nelt(OPLbuf), OPL);
		}
		offset += n;
	}
	return NULL;
}

/* --- .Call ENTRY POINTS ---
 *   - explode_cigar_ops()
 *   - explode_cigar_op_lengths()
 * Args:
 *   cigar: character vector containing the extended CIGAR strings to
 *          explode.
 *   ops:   NULL or a character vector containing the CIGAR operations to
 *          actually consider. If NULL, then all CIGAR operations are
 *          considered.
 * Both functions return a list of the same length as 'cigar' where each
 * list element is a character vector (for explode_cigar_ops()) or an integer
 * vector (for explode_cigar_op_lengths()). The 2 lists have the same shape,
 * that is, same length() and same elementLengths(). The i-th character vector
 * in the list returned by explode_cigar_ops() contains one single-letter
 * string per CIGAR operation in 'cigar[i]'. The i-th integer vector in the
 * list returned by explode_cigar_op_lengths() contains the corresponding
 * CIGAR operation lengths. Zero-length operations or operations not listed
 * in 'ops' are ignored.
 */
SEXP explode_cigar_ops(SEXP cigar, SEXP ops)
{
	SEXP ans, cigar_elt, ans_elt, ans_elt_elt;
	int cigar_len, ans_elt_len, i, j;
	CharAE OPbuf;
	const char *cigar_string, *errmsg;

	cigar_len = LENGTH(cigar);
	init_ops_lkup_table(ops);
	PROTECT(ans = NEW_LIST(cigar_len));
	OPbuf = new_CharAE(0);
	for (i = 0; i < cigar_len; i++) {
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
		CharAE_set_nelt(&OPbuf, 0);
		errmsg = split_cigar_string(cigar_string, &OPbuf, NULL);
		if (errmsg != NULL) {
			UNPROTECT(1);
			error("in 'cigar[%d]': %s", i + 1, errmsg);
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

SEXP explode_cigar_op_lengths(SEXP cigar, SEXP ops)
{
	SEXP ans, cigar_elt, ans_elt;
	int cigar_len, i;
	IntAE OPLbuf;
	const char *cigar_string, *errmsg;

	cigar_len = LENGTH(cigar);
	init_ops_lkup_table(ops);
	PROTECT(ans = NEW_LIST(cigar_len));
	OPLbuf = new_IntAE(0, 0, 0);
	for (i = 0; i < cigar_len; i++) {
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
		IntAE_set_nelt(&OPLbuf, 0);
		errmsg = split_cigar_string(cigar_string, NULL, &OPLbuf);
		if (errmsg != NULL) {
			UNPROTECT(1);
			error("in 'cigar[%d]': %s", i + 1, errmsg);
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
 *
 * TODO: Used by splitCigar() which is now deprecated. Remove this .Call
 * ENTRY POINT once splitCigar() is defunct.
 */
SEXP split_cigar(SEXP cigar)
{
	SEXP ans, cigar_elt, ans_elt, ans_elt_elt0, ans_elt_elt1;
	int cigar_len, i;
	CharAE OPbuf;
	IntAE OPLbuf;
	const char *cigar_string, *errmsg;

	cigar_len = LENGTH(cigar);
	init_ops_lkup_table(R_NilValue);
	PROTECT(ans = NEW_LIST(cigar_len));
	OPbuf = new_CharAE(0);
	OPLbuf = new_IntAE(0, 0, 0);
	for (i = 0; i < cigar_len; i++) {
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
		CharAE_set_nelt(&OPbuf, 0);
		IntAE_set_nelt(&OPLbuf, 0);
		errmsg = split_cigar_string(cigar_string, &OPbuf, &OPLbuf);
		if (errmsg != NULL) {
			UNPROTECT(1);
			error("in 'cigar[%d]': %s", i + 1, errmsg);
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
 * cigar_op_table()
 */

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

/* --- .Call ENTRY POINT ---
 * Args:
 *   cigar: character vector containing the extended CIGAR string for each
 *          read;
 * Return an integer matrix with the number of rows equal to the length of
 * 'cigar' and 9 columns, one for each extended CIGAR operation containing
 * a frequency count for the operations for each element of 'cigar'.
 */
SEXP cigar_op_table(SEXP cigar)
{
	SEXP cigar_string, ans, ans_dimnames, ans_colnames;
	int cigar_len, allOPs_len, i, j, *ans_row;
	const char *allOPs = "MIDNSHP=X", *errmsg;
	char OPstrbuf[2];

	cigar_len = LENGTH(cigar);
	allOPs_len = strlen(allOPs);
	PROTECT(ans = allocMatrix(INTSXP, cigar_len, allOPs_len));
	memset(INTEGER(ans), 0, LENGTH(ans) * sizeof(int));
	ans_row = INTEGER(ans);
	for (i = 0, ans_row = INTEGER(ans); i < cigar_len; i++, ans_row++) {
		cigar_string = STRING_ELT(cigar, i);
		if (cigar_string == NA_STRING) {
			INTEGER(ans)[i] = NA_INTEGER;
			continue;
		}
		errmsg = cigar_string_op_table(cigar_string, allOPs,
				ans_row, cigar_len);
		if (errmsg != NULL) {
			UNPROTECT(1);
			error("in 'cigar[%d]': %s", i + 1, errmsg);
		}
	}

	PROTECT(ans_colnames = NEW_CHARACTER(allOPs_len));
	OPstrbuf[1] = '\0';
	for (j = 0; j < allOPs_len; j++) {
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
 * cigar_ranges()
 */

static SEXP make_list_of_IRanges(const RangeAEAE *range_buf, SEXP names)
{
	SEXP ans, ans_names;

	PROTECT(ans = new_list_of_IRanges_from_RangeAEAE("IRanges", range_buf));
	PROTECT(ans_names = duplicate(names));
	SET_NAMES(ans, ans_names);
	UNPROTECT(2);
	return ans;
}

static SEXP make_CompressedIRangesList(const RangeAE *range_buf,
		const CharAEAE *OP_buf, SEXP breakpoints)
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
					      breakpoints, NULL));
	PROTECT(ans = new_CompressedList(
				"CompressedIRangesList",
				ans_unlistData, ans_partitioning));
	UNPROTECT(3);
	return ans;
}

/* --- .Call ENTRY POINT ---
 * Args:
 *   cigar: character vector containing extended CIGAR strings.
 *   flag:  NULL or an integer vector of the same length as 'cigar'
 *          containing the SAM flag for each read. Serves only as a way to
 *          indicate whether a read is mapped or not. According to the SAM
 *          Spec v1.4, flag bit 0x4 is the only reliable place to tell
 *          whether a segment (or read) is mapped (bit is 0) or not (bit is 1).
 *   space: single integer indicating one of the 8 supported spaces (defined
 *          at the top of this file).
 *   pos:   integer vector of the same length as 'cigar' (or of length 1)
 *          containing the 1-based leftmost position/coordinate of the
 *          clipped read sequences.
 *   f:     NULL or a factor of length 'cigar'. If NULL, then the ranges are
 *          grouped by alignment and stored in a CompressedIRangesList object
 *          with 1 list element per element in 'cigar'. If a factor, then they
 *          are grouped by factor level and stored in an ordinary list of
 *          IRanges objects with 1 list element per level in 'f' and named
 *          with those levels.
 *   ops:   NULL or a character vector containing the CIGAR operations to
 *          translate to ranges. If NULL, then all CIGAR operations are
 *          translated.
 *   drop_empty_ranges: TRUE or FALSE.
 *   reduce_ranges: TRUE or FALSE.
 *   with_ops: TRUE or FALSE indicating whether the returned ranges should be
 *          named with their corresponding CIGAR operation.
 *
 * Returns either a CompressedIRangesList object of the same length as 'cigar'
 * (if 'f' is NULL) or an ordinary list of IRanges objects with 1 list element
 * per level in 'f' (if 'f' is a factor). This list is then turned into a
 * SimpleIRangesList object in R.
 */
SEXP cigar_ranges(SEXP cigar, SEXP flag, SEXP space, SEXP pos, SEXP f,
		  SEXP ops, SEXP drop_empty_ranges, SEXP reduce_ranges,
		  SEXP with_ops)
{
	SEXP ans, ans_breakpoints, f_levels, cigar_elt;
	int cigar_len, space0, pos_len, f_is_NULL, ans_len, *breakpoint,
	    drop_empty_ranges0, reduce_ranges0, with_ops0, i;
	RangeAE range_buf1, *range_buf_p;
	RangeAEAE range_buf2;
	CharAEAE OP_buf, *OP_buf_p;
	const int *flag_elt, *pos_elt, *f_elt;
	const char *cigar_string, *errmsg;

	cigar_len = LENGTH(cigar);
	if (flag != R_NilValue)
		flag_elt = INTEGER(flag);
	init_ops_lkup_table(ops);
	space0 = INTEGER(space)[0];
	pos_len = LENGTH(pos);
	pos_elt = INTEGER(pos);
	f_is_NULL = f == R_NilValue;
	if (f_is_NULL) {
		ans_len = cigar_len;
		/* We will typically generate at least 'cigar_len' ranges. */
		range_buf1 = new_RangeAE(ans_len, 0);
		range_buf_p = &range_buf1;
		PROTECT(ans_breakpoints = NEW_INTEGER(ans_len));
		breakpoint = INTEGER(ans_breakpoints);
	} else {
		f_levels = GET_LEVELS(f);
		ans_len = LENGTH(f_levels);
		range_buf2 = new_RangeAEAE(ans_len, ans_len);
		f_elt = INTEGER(f);
	}
	drop_empty_ranges0 = LOGICAL(drop_empty_ranges)[0];
	reduce_ranges0 = LOGICAL(reduce_ranges)[0];
	with_ops0 = LOGICAL(with_ops)[0];
	if (with_ops0 && f_is_NULL) {
		OP_buf = new_CharAEAE(cigar_len, 0);
		OP_buf_p = &OP_buf;
	} else {
		OP_buf_p = NULL;
	}
	for (i = 0; i < cigar_len; i++) {
		if (flag != R_NilValue) {
			if (*flag_elt == NA_INTEGER) {
				if (f_is_NULL)
					UNPROTECT(1);
				error("'flag' contains NAs");
			}
			if (*flag_elt & 0x004) {
				/* The CIGAR of an unmapped read doesn't
				   produce any range i.e. it's treated as an
				   empty CIGAR. */
				goto for_tail;
			}
		}
		cigar_elt = STRING_ELT(cigar, i);
		if (cigar_elt == NA_STRING) {
			if (f_is_NULL)
				UNPROTECT(1);
			error("'cigar[%d]' is NA", i + 1);
		}
		cigar_string = CHAR(cigar_elt);
		if (strcmp(cigar_string, "*") == 0) {
			if (f_is_NULL)
				UNPROTECT(1);
			error("'cigar[%d]' is \"*\"", i + 1);
		}
		if (*pos_elt == NA_INTEGER || *pos_elt == 0) {
			if (f_is_NULL)
				UNPROTECT(1);
			error("'pos[%d]' is NA or 0", i + 1);
		}
		if (!f_is_NULL) {
			if (*f_elt == NA_INTEGER)
				error("'f[%d]' is NA", i + 1);
			range_buf_p = range_buf2.elts + *f_elt - 1;
		}
		errmsg = parse_cigar_ranges(cigar_string, space0, *pos_elt,
				drop_empty_ranges0, reduce_ranges0,
				range_buf_p, OP_buf_p);
		if (errmsg != NULL) {
			if (f_is_NULL)
				UNPROTECT(1);
			error("in 'cigar[%d]': %s", i + 1, errmsg);
		}
for_tail:
		if (flag != R_NilValue)
			flag_elt++;
		if (pos_len != 1)
			pos_elt++;
		if (f_is_NULL)
			*(breakpoint++) = RangeAE_get_nelt(range_buf_p);
		else
			f_elt++;
	}
	if (!f_is_NULL)
		return make_list_of_IRanges(&range_buf2, f_levels);
	PROTECT(ans = make_CompressedIRangesList(range_buf_p, OP_buf_p,
						 ans_breakpoints));
	UNPROTECT(2);
	return ans;
}


/****************************************************************************
 * --- .Call ENTRY POINT ---
 * Args:
 *   cigar, flag, space: see cigar_ranges() function above.
 * Return an integer vector of the same length as 'cigar' containing the
 * widths of the alignments as inferred from the cigar information.
 */
SEXP cigar_width(SEXP cigar, SEXP flag, SEXP space)
{
	SEXP ans, cigar_elt;
	int cigar_len, space0, i, *ans_elt;
	const int *flag_elt;
	const char *cigar_string, *errmsg;

	cigar_len = LENGTH(cigar);
	if (flag != R_NilValue)
		flag_elt = INTEGER(flag);
	space0 = INTEGER(space)[0];
	PROTECT(ans = NEW_INTEGER(cigar_len));
	for (i = 0, ans_elt = INTEGER(ans); i < cigar_len; i++, ans_elt++) {
		if (flag != R_NilValue) {
			if (*flag_elt == NA_INTEGER) {
				UNPROTECT(1);
				error("'flag' contains NAs");
			}
			if (*flag_elt & 0x004) {
				*ans_elt = NA_INTEGER;
				goto for_tail;
			}
		}
		cigar_elt = STRING_ELT(cigar, i);
		if (cigar_elt == NA_STRING) {
			*ans_elt = NA_INTEGER;
			goto for_tail;
		}
		cigar_string = CHAR(cigar_elt);
		if (strcmp(cigar_string, "*") == 0) {
			*ans_elt = NA_INTEGER;
			goto for_tail;
		}
		errmsg = parse_cigar_width(cigar_string, space0, ans_elt);
		if (errmsg != NULL) {
			UNPROTECT(1);
			error("in 'cigar[%d]': %s", i + 1, errmsg);
		}
for_tail:
		if (flag != R_NilValue)
			flag_elt++;
	}
	UNPROTECT(1);
	return ans;
}


/****************************************************************************
 * cigar_narrow()
 */

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

/* --- .Call ENTRY POINT --- */
SEXP cigar_narrow(SEXP cigar, SEXP left_width, SEXP right_width)
{
	SEXP ans, ans_cigar, ans_cigar_string, ans_rshift, cigar_string;
	int cigar_len, i;
	static char cigar_buf[1024];
	const char *errmsg;

	cigar_len = LENGTH(cigar);
	PROTECT(ans_cigar = NEW_CHARACTER(cigar_len));
	PROTECT(ans_rshift = NEW_INTEGER(cigar_len));
	for (i = 0; i < cigar_len; i++) {
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
			error("in 'cigar[%d]': %s", i + 1, errmsg);
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
 * cigar_qnarrow()
 */

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

/* --- .Call ENTRY POINT ---
 * Return a list of 2 elements: 1st elt is the narrowed cigar vector, 2nd elt
 * is the 'rshift' vector i.e. the integer vector of the same length as 'cigar'
 * that would need to be added to the 'pos' field of a SAM/BAM file as a
 * consequence of this narrowing.
 */
SEXP cigar_qnarrow(SEXP cigar, SEXP left_qwidth, SEXP right_qwidth)
{
	SEXP ans, ans_cigar, ans_cigar_string, ans_rshift, cigar_string;
	int cigar_len, i;
	static char cigar_buf[1024];
	const char *errmsg;

	cigar_len = LENGTH(cigar);
	PROTECT(ans_cigar = NEW_CHARACTER(cigar_len));
	PROTECT(ans_rshift = NEW_INTEGER(cigar_len));
	for (i = 0; i < cigar_len; i++) {
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
			error("in 'cigar[%d]': %s", i + 1, errmsg);
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

