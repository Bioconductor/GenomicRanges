#include "GenomicRanges.h"
#include "IRanges_interface.h"

#include <ctype.h> /* for isdigit() */

static char errmsg_buf[200];

/* Return the number of chars that was read, or 0 if there is no more char
   to read (i.e. cig0[offset] is '\0'), or -1 in case of a parse error.
   Zero-length operations are ignored. */
static int get_next_cigar_OP(const char *cig0, int offset,
		int *OPL, char *OP)
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
				 "unexpected CIGAR end at char %d",
				 offset + 1);
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
static int get_prev_cigar_OP(const char *cig0, int offset,
		int *OPL, char *OP)
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
				 "no CIGAR operation length at char %d",
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
	while ((n = get_next_cigar_OP(cig0, offset, &OPL, &OP))) {
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

static const char *cigar_string_to_qwidth(SEXP cigar_string, int clip_reads,
		int *qwidth)
{
	const char *cig0;
	int offset, n, OPL /* Operation Length */;
	char OP /* Operation */;

	if (cigar_string == NA_STRING)
		return "CIGAR string is NA";
	if (LENGTH(cigar_string) == 0)
		return "CIGAR string is empty";
	cig0 = CHAR(cigar_string);
	*qwidth = offset = 0;
	while ((n = get_next_cigar_OP(cig0, offset, &OPL, &OP))) {
		if (n == -1)
			return errmsg_buf;
		switch (OP) {
		/* Alignment match (can be a sequence match or mismatch) */
		    case 'M': case '=': case 'X': *qwidth += OPL; break;
		/* Insertion to the reference */
		    case 'I': *qwidth += OPL; break;
		/* Deletion (or skipped region) from the reference,
		   or silent deletion from the padded reference */
		    case 'D': case 'N': case 'P': break;
		/* Soft clip on the read */
		    case 'S': *qwidth += OPL; break;
		/* Hard clip on the read */
		    case 'H': if (!clip_reads) *qwidth += OPL; break;
		    default:
			snprintf(errmsg_buf, sizeof(errmsg_buf),
				 "unknown CIGAR operation '%c' at char %d",
				 OP, offset + 1);
			return errmsg_buf;
		}
		offset += n;
	}
	return NULL;
}

static const char *cigar_string_to_width(SEXP cigar_string, int *width)
{
	const char *cig0;
	int offset, n, OPL /* Operation Length */;
	char OP /* Operation */;

	if (cigar_string == NA_STRING)
		return "CIGAR string is NA";
	if (LENGTH(cigar_string) == 0)
		return "CIGAR string is empty";
	cig0 = CHAR(cigar_string);
	*width = offset = 0;
	while ((n = get_next_cigar_OP(cig0, offset, &OPL, &OP))) {
		if (n == -1)
			return errmsg_buf;
		switch (OP) {
		/* Alignment match (can be a sequence match or mismatch) */
		    case 'M': case '=': case 'X': *width += OPL; break;
		/* Insertion to the reference */
		    case 'I': break;
		/* Deletion (or skipped region) from the reference */
		    case 'D': case 'N': *width += OPL; break;
		/* Soft/Hard clip on the read */
		    case 'S': case 'H': break;
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
	while ((n = get_next_cigar_OP(cig0, offset, &OPL, &OP))) {
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
	while ((n = get_prev_cigar_OP(cig0, offset, &OPL, &OP))) {
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
		/* Deletion (or skipped region) from the reference */
		    case 'D': case 'N':
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
		n = get_next_cigar_OP(cig0, offset, &OPL, &OP);
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
	while ((n = get_next_cigar_OP(cig0, offset, &OPL, &OP))) {
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
	while ((n = get_prev_cigar_OP(cig0, offset, &OPL, &OP))) {
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
		n = get_next_cigar_OP(cig0, offset, &OPL, &OP);
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

static const char *split_cigar_string(SEXP cigar_string,
		CharAE *OPbuf, IntAE *OPLbuf)
{
	const char *cig0;
	int offset, n, OPL /* Operation Length */;
	char OP /* Operation */;

	cig0 = CHAR(cigar_string);
	offset = 0;
	while ((n = get_next_cigar_OP(cig0, offset, &OPL, &OP))) {
		if (n == -1)
			return errmsg_buf;
		CharAE_insert_at(OPbuf, CharAE_get_nelt(OPbuf), OP);
		IntAE_insert_at(OPLbuf, IntAE_get_nelt(OPLbuf), OPL);
		offset += n;
	}
	return NULL;
}

static void append_range(RangeAE *range_ae, int start, int width)
{
	RangeAE_insert_at(range_ae, RangeAE_get_nelt(range_ae), start, width);
}

static const char *cigar_string_to_ranges(SEXP cigar_string, int pos_elt,
		int Ds_as_Ns, RangeAE *range_ae)
{
	const char *cig0;
	int offset, n, OPL /* Operation Length */, start, width;
	char OP /* Operation */;

	cig0 = CHAR(cigar_string);
	offset = 0;
	start = pos_elt;
	while ((n = get_next_cigar_OP(cig0, offset, &OPL, &OP))) {
		if (n == -1)
			return errmsg_buf;
		width = 0;
		switch (OP) {
		/* Alignment match (can be a sequence match or mismatch) */
		    case 'M': case '=': case 'X': width = OPL; break;
		/* Insertion to the reference */
		    case 'I': break;
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
		if (width) {
			append_range(range_ae, start, width);
			start += width;
		}
		offset += n;
	}
	return NULL;
}

/* Unlike cigar_string_to_ranges(), cigar_string_to_ranges2() merges adjacent
   ranges. */
static const char *cigar_string_to_ranges2(SEXP cigar_string, int pos_elt,
		int Ds_as_Ns, RangeAE *range_ae)
{
	const char *cig0;
	int offset, n, OPL /* Operation Length */, start, width;
	char OP /* Operation */;

	cig0 = CHAR(cigar_string);
	offset = 0;
	start = pos_elt;
	width = 0;
	while ((n = get_next_cigar_OP(cig0, offset, &OPL, &OP))) {
		if (n == -1)
			return errmsg_buf;
		switch (OP) {
		/* Alignment match (can be a sequence match or mismatch) */
		    case 'M': case '=': case 'X': width += OPL; break;
		/* Insertion to the reference */
		    case 'I': break;
		/* Deletion from the reference */
		    case 'D':
			if (Ds_as_Ns) {
				if (width)
					append_range(range_ae, start, width);
				start += width + OPL;
				width = 0;
			} else {
				width += OPL;
			}
		    break;
		/* Skipped region from the reference */
		    case 'N':
			if (width)
				append_range(range_ae, start, width);
			start += width + OPL;
			width = 0;
		    break;
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
		offset += n;
	}
	if (width)
		append_range(range_ae, start, width);
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
	SEXP ans, cigar_string;
	int cigar_length, ans_type0, i, qwidth;
	const char *errmsg;
	char string_buf[200];

	cigar_length = LENGTH(cigar);
	ans_type0 = INTEGER(ans_type)[0];
	if (ans_type0 == 1)
		PROTECT(ans = NEW_LOGICAL(cigar_length));
	else
		ans = R_NilValue;
	for (i = 0; i < cigar_length; i++) {
		cigar_string = STRING_ELT(cigar, i);
		/* we use cigar_string_to_qwidth() here just for its ability
                   to parse and detect ill-formed CIGAR strings */
		errmsg = cigar_string_to_qwidth(cigar_string, 1, &qwidth);
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
 * --- .Call ENTRY POINT ---
 * Args:
 *   cigar: character vector containing the extended CIGAR string for each
 *          read.
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
 *   cigar: character vector containing the extended CIGAR string for each
 *          read;
 *   before_hard_clipping: TRUE or FALSE indicating whether the returned
 *          widths should be those of the reads before or after "hard
 *          clipping".
 * Return an integer vector of the same length as 'cigar' containing the
 * lengths of the query sequences as inferred from the cigar information.
 */
SEXP cigar_to_qwidth(SEXP cigar, SEXP before_hard_clipping)
{
	SEXP ans, cigar_string;
	int clip_reads, cigar_length, i, qwidth;
	const char *errmsg;

	clip_reads = !LOGICAL(before_hard_clipping)[0];
	cigar_length = LENGTH(cigar);
	PROTECT(ans = NEW_INTEGER(cigar_length));
	for (i = 0; i < cigar_length; i++) {
		cigar_string = STRING_ELT(cigar, i);
		if (cigar_string == NA_STRING) {
			INTEGER(ans)[i] = NA_INTEGER;
			continue;
		}
		errmsg = cigar_string_to_qwidth(cigar_string, clip_reads,
				&qwidth);
		if (errmsg != NULL) {
			UNPROTECT(1);
			error("in 'cigar' element %d: %s", i + 1, errmsg);
		}
		INTEGER(ans)[i] = qwidth;
	}
	UNPROTECT(1);
	return ans;
}


/****************************************************************************
 * --- .Call ENTRY POINT ---
 * Args:
 *   cigar: character vector containing the extended CIGAR string for each
 *          read.
 * Return an integer vector of the same length as 'cigar' containing the
 * widths of the alignments as inferred from the cigar information.
 */
SEXP cigar_to_width(SEXP cigar)
{
	SEXP ans, cigar_string;
	int cigar_length, i, width;
	const char *errmsg;

	cigar_length = LENGTH(cigar);
	PROTECT(ans = NEW_INTEGER(cigar_length));
	for (i = 0; i < cigar_length; i++) {
		cigar_string = STRING_ELT(cigar, i);
		if (cigar_string == NA_STRING) {
			INTEGER(ans)[i] = NA_INTEGER;
			continue;
		}
		errmsg = cigar_string_to_width(cigar_string, &width);
		if (errmsg != NULL) {
			UNPROTECT(1);
			error("in 'cigar' element %d: %s", i + 1, errmsg);
		}
		INTEGER(ans)[i] = width;
	}
	UNPROTECT(1);
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
 * Args:
 *   cigar: character string containing the extended CIGAR;
 *   drop_D_ranges: TRUE or FALSE indicating whether Ds should be treated
 *          like Ns or not;
 *   merge_ranges: TRUE or FALSE indicating whether adjacent ranges coming
 *          from the same cigar should be merged or not.
 * Return an IRanges object describing the alignment.
 */
SEXP cigar_to_IRanges(SEXP cigar, SEXP drop_D_ranges, SEXP merge_ranges)
{
	RangeAE range_ae;
	SEXP cigar_string;
	int Ds_as_Ns, merge_ranges0;
	const char *errmsg;

	cigar_string = STRING_ELT(cigar, 0);
	if (cigar_string == NA_STRING)
		error("'cigar' is NA");
	Ds_as_Ns = LOGICAL(drop_D_ranges)[0];
	merge_ranges0 = LOGICAL(merge_ranges)[0];
	range_ae = new_RangeAE(0, 0);
	errmsg = merge_ranges0 ?
			cigar_string_to_ranges2(cigar_string, 1,
				Ds_as_Ns, &range_ae) :
			cigar_string_to_ranges(cigar_string, 1,
				Ds_as_Ns, &range_ae);
	if (errmsg != NULL)
		error("%s", errmsg);
	return new_IRanges_from_RangeAE("IRanges", &range_ae);
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
 * Returns a CompressedNormalIRangesList object of the same length as the
 * input.
 * NOTE: See note for cigar_to_list_of_IRanges_by_rname() below about the
 * strand.
 * TODO: Support character factor 'cigar' in addition to current character
 *       vector format.
 */
SEXP cigar_to_list_of_IRanges_by_alignment(SEXP cigar, SEXP pos, SEXP flag,
		SEXP drop_D_ranges)
{
	SEXP cigar_string;
	SEXP ans, ans_unlistData, ans_partitioning, ans_partitioning_end;
	int cigar_length, Ds_as_Ns, i, pos_elt, flag_elt;
	RangeAE range_ae;
	const char *errmsg;

	cigar_length = LENGTH(cigar);
	Ds_as_Ns = LOGICAL(drop_D_ranges)[0];
	/* we will generate at least 'cigar_length' ranges, and possibly more */
	range_ae = new_RangeAE(cigar_length, 0);
	PROTECT(ans_partitioning_end = NEW_INTEGER(cigar_length));
	for (i = 0; i < cigar_length; i++) {
		if (flag != R_NilValue) {
			flag_elt = INTEGER(flag)[i];
			if (flag_elt == NA_INTEGER) {
				UNPROTECT(1);
				error("'flag' contains NAs");
			}
			if (flag_elt & 0x004)
				continue;
		}
		cigar_string = STRING_ELT(cigar, i);
		if (cigar_string == NA_STRING) {
			UNPROTECT(1);
			error("'cigar' contains %sNAs",
			      flag != R_NilValue ? "unexpected " : "");
		}
		pos_elt = INTEGER(pos)[i];
		if (pos_elt == NA_INTEGER) {
			UNPROTECT(1);
			error("'pos' contains %sNAs",
			      flag != R_NilValue ? "unexpected " : "");
		}
		errmsg = cigar_string_to_ranges2(cigar_string, pos_elt,
				Ds_as_Ns, &range_ae);
		if (errmsg != NULL) {
			UNPROTECT(1);
			error("in 'cigar' element %d: %s", i + 1, errmsg);
		}
		INTEGER(ans_partitioning_end)[i] = RangeAE_get_nelt(&range_ae);
	}
	PROTECT(ans_unlistData = new_IRanges_from_RangeAE(
			"IRanges", &range_ae));
	PROTECT(ans_partitioning = new_PartitioningByEnd(
			"PartitioningByEnd",
			ans_partitioning_end, NULL));
	PROTECT(ans = new_CompressedList(
			"CompressedNormalIRangesList",
			ans_unlistData, ans_partitioning));
	UNPROTECT(4);
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
 *   merge_ranges: TRUE or FALSE indicating whether adjacent ranges coming
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
		SEXP flag, SEXP drop_D_ranges, SEXP merge_ranges)
{
	SEXP rname_levels, cigar_string, ans, ans_names;
	int ans_length, nreads, Ds_as_Ns, merge_ranges0,
	    i, level, pos_elt, flag_elt;
	RangeAEAE range_aeae;
	const char *errmsg;

	rname_levels = GET_LEVELS(rname);
	ans_length = LENGTH(rname_levels);
	range_aeae = new_RangeAEAE(ans_length, ans_length);
	nreads = LENGTH(pos);
	Ds_as_Ns = LOGICAL(drop_D_ranges)[0];
	merge_ranges0 = LOGICAL(merge_ranges)[0];
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
		pos_elt = INTEGER(pos)[i];
		if (pos_elt == NA_INTEGER)
			error("'pos' contains %sNAs",
			      flag != R_NilValue ? "unexpected " : "");
		errmsg = merge_ranges0 ?
			cigar_string_to_ranges2(cigar_string, pos_elt,
				Ds_as_Ns, range_aeae.elts + level - 1) :
			cigar_string_to_ranges(cigar_string, pos_elt,
				Ds_as_Ns, range_aeae.elts + level - 1);
		if (errmsg != NULL)
			error("in 'cigar' element %d: %s", i + 1, errmsg);
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
 *   ref_locs: global positions in the reference that we will map
 *   cigar: character string containing the extended CIGAR;
 *   pos: reference position at which the query alignment begins
 *        (after clip)
 * Returns the local query positions. This assumes that the reference
 * positions actually occur in the read alignment region, outside of
 * any deletions or insertions. 
 */
SEXP ref_locs_to_query_locs(SEXP ref_locs, SEXP cigar, SEXP pos)
{
  int nlocs, i;
  SEXP query_locs;
  
  nlocs = LENGTH(ref_locs);
  PROTECT(query_locs = allocVector(INTSXP, nlocs));
  
  for (i = 0; i < nlocs; i++) {
    int query_loc = INTEGER(ref_locs)[i] - INTEGER(pos)[i] + 1;
    int n, offset = 0, OPL, query_consumed = 0;
    char OP;
    const char *cig0 = CHAR(STRING_ELT(cigar, i));
    while (query_consumed < query_loc &&
           (n = get_next_cigar_OP(cig0, offset, &OPL, &OP)))
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
      case 'N': /* Skipped region from the reference */
        query_loc -= OPL;
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
