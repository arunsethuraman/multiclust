/**
 * @file read_file.c
 * @author Arun Sethuraman
 * @author Karin Dorman, kdorman@iastate.edu
 * @date Thu Dec  6 08:21:03 CST 2012
 *
 * Functions for reading in data.
 *
 * TODO
 * - Are alleles in STRUCTURE files always stored as integers?
 */

#include "multiclust.h"
#define MAKE_1ARRAY MAKE_1ARRAY_RETURN	/* return on memory allocation error */

int summarize_alleles(options *opt, data *dat);
int sufficient_statistics(data *dat);

/* perhaps the functions below belong in a utility file */
int count_columns(FILE *fp);
int count_lines(FILE *fp);
void skip_line(FILE *fp);
int read_next_word(FILE *fp, char **word);
int add_to_string_set(char ***list, int *len, char *str, int *loc);
void bubbleSort(int *numbers, int array_size);

/**
 * Read the input data file.  This function reads STRUCTURE (Pritchard et al.,
 * 2000) formatted input files and stores the data in a data object.
 *
 * Note [BUG?], it assumes all alleles are integers.
 *
 * @param datafile name of the input file
 * @param dat data object
 * @return error status
 */
int read_file(options *opt, data *dat)
{
	FILE *f1;		/* file handle */
	int n_info_col = 2;	/* no. non-loci columns */
	int i, j, l, idv;	/* indices */
	int n_haplotypes;	/* total number of observed haplotypes */
	int debug = 0;		/* debugging level */
	char *locale = NULL;	/* read-in locale */
	char *name1 = NULL, *name2 = NULL;
	int err = NO_ERROR;
	char a;
	
	/* open file or abort */
	if ((f1 = fopen(opt->filename, "r")) == NULL)
		return message(stderr, __FILE__, __func__, __LINE__,
			ERROR_MSG, FILE_OPEN_ERROR, opt->filename);

	/* count number of loci (or loci * ploidy)*/
	dat->L = count_columns(f1);	/* f1 on first char of next line */
	if (opt->R_format)
		dat->L -= 2; /* KLUDGE: uses R-formatted Structure file format */

	if (opt->verbosity >= TALKATIVE)
		mmessage(INFO_MSG, NO_ERROR, "Input file '%s' has %u "
				"columns.\n", opt->filename, dat->L + 2);

	if ((err = read_next_word(f1, &name1)))
		return message(stderr, __FILE__, __func__, __LINE__, ERROR_MSG,
			err, "Failed to read first individual's name in file "
						"'%s'\n", opt->filename);
	skip_line(f1);
	if ((err = read_next_word(f1, &name2)))
		return message(stderr, __FILE__, __func__, __LINE__, ERROR_MSG,
			err, "Failed to read second individual's name in file "
						"'%s'\n", opt->filename);

	if (strcmp(name1, name2))
		opt->interleaved = 1;

	free(name1);
	free(name2);

	/* count number of individuals sampled */
	dat->I = count_lines(f1) + 1;	/* f1 now at EOF */

	if (opt->verbosity >= TALKATIVE)
		mmessage(INFO_MSG, NO_ERROR, "Input file '%s' has %u "
					"lines.\n", opt->filename, dat->I);

	/* self-consistency check on amount of data */
	if (!opt->interleaved && (dat->I % dat->ploidy))
		return message(stderr, __FILE__, __func__, __LINE__,
			ERROR_MSG, FILE_FORMAT_ERROR, "number of lines (%d) in "
			"'%s' is not a multiple of ploidy (%d)\n", dat->I,
						opt->filename, dat->ploidy);
	else if (opt->interleaved && (dat->L % dat->ploidy))
		return message(stderr, __FILE__, __func__, __LINE__,
			ERROR_MSG, FILE_FORMAT_ERROR, "number of columns (%u) "
			"in '%s' is not a multiple of ploidy (%d)\n", dat->L,
						opt->filename, dat->ploidy);

	/* allocate memory for raw data */
	if (opt->interleaved) {
		n_haplotypes = dat->I * dat->ploidy;
		dat->L /= dat->ploidy;
	} else {
		n_haplotypes = dat->I;
		dat->I /= dat->ploidy;
	}

	if (opt->verbosity >= TALKATIVE)
		mmessage(INFO_MSG, NO_ERROR, "Input file '%s' has %u "
				"haplotypes, %u individuals.\n", opt->filename,
							n_haplotypes, dat->I);


	MAKE_2ARRAY(dat->IL, n_haplotypes, dat->L);

	if (debug)
		message(stderr, __FILE__, __func__, __LINE__, ERROR_MSG,
			FILE_FORMAT_ERROR, "read %d individuals\n", dat->I);

	/* and individual names & locales */
	MAKE_1ARRAY(dat->idv, dat->I);
	MAKE_1ARRAY(dat->I_K, dat->I);
	if (opt->admixture)
		MAKE_2ARRAY(dat->IL_K, dat->I * dat->ploidy, dat->L);
	
	for (i=0; i<dat->I; i++)
		dat->idv[i].name = NULL;
	idv = 0;
	dat->numpops = 0;

	/* read alleles */
	fseek(f1, 0, SEEK_SET);
	skip_line(f1);

	for (i=0; i<n_haplotypes; i += opt->interleaved ? 1 : dat->ploidy) {

		/* first haplotype of individual idv */
		if (opt->interleaved || !(i % dat->ploidy)) {

			/* read name of individual */
			if ((err = read_next_word(f1, &(dat->idv[idv].name)))) {
				if (feof(f1))
					break;
				return err;
			}

			/* read locale of individual */
			if ((err = read_next_word(f1, &locale)))
				return err;

			/* add (new) locale to locale set */
			err = add_to_string_set(&(dat->pops), &dat->numpops,
				locale, &(dat->idv[idv].locale));
			free(locale);
			if (err)
				return err;

			idv++;

		/* other haplotypes of idv: first cols. are duplicate */
		} else {

			/* skip non-data columns */
			for (j=0; j<n_info_col; j++) {
				/* skip white space before next col. */
				do {
					a = fgetc(f1);
				} while (!feof(f1) && (a == ' ' || a == '\t'
					|| a == '\n'));
	
				/* skip to end of next column */
				while (!feof(f1) && a != ' ' && a != '\t'
					&& a != '\n') {
					a = fgetc(f1);
				}
			}
		}
		if (feof(f1))
			break;

		for (l=0; l<dat->L; l++) {
			for (j = 0; j < (opt->interleaved ? dat->ploidy : 1); ++j)
				if (fscanf(f1, "%d", &(dat->IL[i + j][l])) != 1)
					return message(stderr, __FILE__,
						__func__, __LINE__, ERROR_MSG,
						FILE_FORMAT_ERROR, "failed to "
						"read locus %d of haplotype %d "
						"in file '%s'.  Check option "
						"-R.\n", l + 1, i + j + 1,
								opt->filename);
			fgetc(f1);	/* clear new line */
		}
	
		if ((opt->interleaved || !(i % dat->ploidy))
						&& opt->verbosity >= VERBOSE)
			mmessage(INFO_MSG, NO_ERROR, "Read individual %u.\n",
				i / (opt->interleaved ? 1 : dat->ploidy) + 1);

		if (feof(f1))
			break;

	}

	fclose(f1);

	if (opt->verbosity >= TALKATIVE)
		mmessage(INFO_MSG, NO_ERROR, "Finished reading %u populations "
				"from '%s'.\n", dat->numpops, opt->filename);


	/* count no. individuals observed at each locale */
	CMAKE_1ARRAY(dat->i_p, dat->numpops);
	
	for (i=0; i<dat->I; i++)
		dat->i_p[dat->idv[i].locale]++;

	if (opt->verbosity >= TALKATIVE)
		mmessage(INFO_MSG, NO_ERROR, "Finished computing population "
							"sample sizes.\n");

	if (debug)
		for (i=0; i<dat->I; i++)
			fprintf(stderr, "%s: individual %s sampled at locale "
				"%s\n", __func__, dat->idv[i].name,
				dat->pops[dat->idv[i].locale]);

//At this point, I have an (I*ploidy)xL matrix which contains all the alleles
//from all individuals.  I need to extract the number of unique alleles at each
//locus, then build an IxLxM_l jagged 3D matrix, where m_l is the number of
//distinct alleles at locus l.  Each element of this matrix should contain the
//count of the number of that unique allele at locus L in that individual I.

	if ((err = summarize_alleles(opt, dat)))
                return err;

	if (opt->verbosity >= TALKATIVE)
		mmessage(INFO_MSG, NO_ERROR, "Finished summarizing alleles: "
			"maximum number of alleles at a single locus is %u.\n",
									dat->M);

	if ((err = sufficient_statistics(dat)))
		return err;

	if (opt->verbosity >= TALKATIVE)
		mmessage(INFO_MSG, NO_ERROR, "Finished computing sufficient "
							"statistics.\n");

	return err;

} /* End of read_file(). */


int make_ila(data *dat) {
	int i, l, a, m, m_start;

	MAKE_3ARRAY(dat->ila, dat->I, dat->L, dat->ploidy);
	if (!dat->ila)
		return message(stderr, __FILE__, __func__, __LINE__,
			ERROR_MSG, MEMORY_ALLOCATION,
			"dat->ila[%d][%d][%d]\n", dat->I, dat->L, dat->ploidy);

	for (i = 0; i < dat->I; i++) {
		for (l = 0; l < dat->L; l++) {
			m_start = dat->L_alleles && dat->L_alleles[l][0] == MISSING;
			for (a = 0; a < dat->ploidy; a++) {
				for (m = m_start; m < dat->uniquealleles[l];
					m++) {
					if (dat->L_alleles[l][m]
						== dat->IL[i*2 + a][l]) {
						dat->ila[i][l][a] = m;
						break;
					}
				}
			}
		}
	}

	return NO_ERROR;
} /* make_ila */


/**
 * Allocate and populate matrices for allele information.  Allocates and
 * populates uniquealleles, which stores number of unique alleles at each locus.
 * Allocates and populates L_alleles, a jagged array that lists number of
 * unique alleles at each locus.
 *
 * NOTE: hidden alleles coded as -9; hence sort first in the allele lists
 *
 * @param dat data object
 * @return error status
 */
int summarize_alleles(options *opt, data *dat)
{
	int i, k, l, m;
	int n_haplotypes
		= dat->I * dat->ploidy;	/* no. of haplotypes */
	int *locusgeno = NULL;		/* one column of data from dat->IL */

	/* count unique alleles at each locus */
	CMAKE_1ARRAY(dat->uniquealleles, dat->L);

	if (!opt->alleles_are_indices)
		MAKE_1ARRAY(locusgeno, n_haplotypes);

	if (opt->verbosity >= TALKATIVE)
		mmessage(INFO_MSG, NO_ERROR, "Processing loci ");

	dat->M = 0;
	for (l = 0; l < dat->L; l++) {
		if (opt->alleles_are_indices) {
			for (i=0; i < n_haplotypes; i++)
				dat->uniquealleles[l] 
					= MAX(dat->uniquealleles[l],
						dat->IL[i][l] + 1);
			if (opt->verbosity > MINIMAL)
				fprintf(stderr, "Site %d: %4d alleles\n",
					l+1, dat->uniquealleles[l]);
		} else {
			/* copy locus l observations to locusgeno */
			for (i = 0; i < n_haplotypes; i++)
				locusgeno[i] = dat->IL[i][l];
	
			bubbleSort(locusgeno, n_haplotypes);
	
			/* count unique alleles, minus missing (was BUG) */
			k = 0;
			while (k < n_haplotypes && locusgeno[k] == MISSING)
				k++;
			dat->uniquealleles[l] = 0;
			if (k == n_haplotypes)
				continue;
			if (locusgeno[0] == MISSING) {
				dat->missing_data = 1;
				dat->uniquealleles[l] = 2;
			} else
				dat->uniquealleles[l] = 1;
			for (k++; k < n_haplotypes; k++)
				if (locusgeno[k] != locusgeno[k - 1])
					dat->uniquealleles[l]++;
	
		}
		if (dat->uniquealleles[l] > dat->M)
			dat->M = dat->uniquealleles[l];
		if (opt->verbosity >= TALKATIVE && !(l % 100))
			fprintf(stderr, ".");
	}
	if (opt->verbosity >= TALKATIVE)
		fprintf(stderr, "\n");

	if (!opt->alleles_are_indices) {
		if (opt->verbosity >= TALKATIVE)
			mmessage(INFO_MSG, NO_ERROR, "Sorting alleles ");

		/* allocate jagged array listing unique alleles at each locus */
		MAKE_2JAGGED_ARRAY(dat->L_alleles, dat->L, dat->uniquealleles);
	
		for (l = 0; l < dat->L; l++) {
			/* copy locus l observations to locusgeno */
			for (i = 0; i < n_haplotypes; i++)
				locusgeno[i] = dat->IL[i][l];
	
			bubbleSort(locusgeno, n_haplotypes);
	
			/* copy unique alleles into L_alleles */
			m = 0;
			dat->L_alleles[l][m++] = locusgeno[0];
			for (k = 1; k < n_haplotypes; k++)
				if (locusgeno[k] != locusgeno[k - 1])
					dat->L_alleles[l][m++] = locusgeno[k];
			if (opt->verbosity >= TALKATIVE && !(l % 100))
				fprintf(stderr, ".");
		}
		if (opt->verbosity >= TALKATIVE)
			fprintf(stderr, "\n");
	}

	/* could now rewrite IL to be indices of L_alleles */
	
	/* free up all space */
	if (!opt->alleles_are_indices)
		FREE_1ARRAY(locusgeno);

	return NO_ERROR;
} /* End of summarize_alleles(). */

void print_unique_alleles(data *dat)
{
	int l, j;
	fprintf(stdout, "print_unique_alleles():\n");
	for (l = 0; l < dat->L; l++) {
		for (j=0; j < dat->uniquealleles[l]; j++)
			fprintf(stdout, "%d\t", dat->L_alleles[l][j]);
		fprintf(stdout, "\n");
	}
} /* print_unique_alleles */

void print_number_unique_alleles(data *dat)
{
	int l;
	fprintf(stdout, "print_number_unique_alleles():\n");
	for (l = 0; l < dat->L; l++) {
		fprintf(stdout, "%d\t", dat->uniquealleles[l]);
	}
	fprintf(stdout, "\n");
} /* print_number_unique_alleles */


/**
 * Allocate and populate the ILM matrix.  The ILM matrix stores the count of
 * all unique alleles at each locus in all individuals, the sufficient 
 * statistics.
 *
 * @param dat data object
 * @return void
 */
int sufficient_statistics(data *dat)
{
	int i, j = 0, l = 0, m = 0, z_start, z_end;

	CMAKE_3JAGGED_ARRAY(dat->ILM, dat->I, dat->L, dat->uniquealleles);

	for (l = 0; l < dat->L; l++) {				/* locus */
		z_end = 0;
		for (i = 0; i < dat->I; i++) {			/* individual */
			z_start = z_end;
			z_end += dat->ploidy;
			if (!dat->L_alleles)
				for (j = z_start; j < z_end; j++)
					dat->ILM[i][l][dat->IL[j][l]]++;
			else {
				for (m = 0; m < dat->uniquealleles[l]; m++)
								/* allele */
					for (j = z_start; j < z_end; j++)
								/* haplotype */
						if (dat->IL[j][l]
							== dat->L_alleles[l][m])
							dat->ILM[i][l][m]++;
			}
		}
	}

	return NO_ERROR;
} /* End of sufficient_statistics(). */

void print_sufficient_statistics(data *dat)
{
	int i, l, m;
	fprintf(stdout, "print_sufficient_statistics():\n");
	for (i = 0; i < dat->I; i++) {
		fprintf(stdout, "%d:\n", i);
		for (l = 0; l < dat->L; l++) {
			for (m = 0; m < dat->uniquealleles[l]; m++)
				fprintf(stdout, "\t%d", dat->ILM[i][l][m]);
			fprintf(stdout, "\n");
		}
		fprintf(stdout, "\n");
	}
} /* print_sufficient_statistics */


/**
 * Count number of columns in a file.  This function assumes the file pointer
 * is set to the beginning of the first line.  It counts columns by assuming
 * any amount of white space (space or tab) separates columns.  Leading white
 * space is ignored.  It leaves the file pointer at the beginning of the next
 * line (or at EOF).
 *
 * @param fp file handle
 * @return number of columns
 */
int count_columns(FILE *fp)
{
	char a;
	int ncol = 0;

	/* skip any white space to first column header */
	do {
		a = fgetc(fp);
	} while (!feof(fp) && (a == ' ' || a == '\t'));

	while (!feof(fp)) {
		/* skip to end of column header */
		while (!feof(fp) && a != ' ' && a != '\t' && a != '\n')
			a = fgetc(fp);
		ncol++;
		if (a == '\n')	/* last column */
			break;
		/* skip through white space to next column  header */
		while (!feof(fp) && (a == ' ' || a == '\t'))
			a = fgetc(fp);
	}
	return ncol;
} /* count_columns */

/**
 * Count number of non-empty lines before end of file.  This function counts
 * how many non-empty lines (lines with more than just white space) before EOF.
 * Make sure your file pointer is at the beginning of the file if you want to
 * the total number of lines in the file.  If the file pointer is on any 
 * character of non-empty line other than the terminating newline, it will
 * include the current line in the count.
 *
 * @param fp file handle
 * @return number of lines
 */
int count_lines(FILE *fp)
{
	char a;			/* store characters off stream */
	int nlines = 0;

	while (!feof(fp)) {
		/* skip white space until content  
		 * (note: last line need not end w/ newline) */
		do {
			a = fgetc(fp);
		} while (!feof(fp) && (a == ' ' || a == '\t' || a == '\n'
			|| a == '\r'));	/* handle \n\r<EOF> */

		/* found another line of content */
		if (!feof(fp)) {
			nlines++;
			/* skip to end of line */
			while (!feof(fp) && a != '\n')
				a = fgetc(fp);
		}
	}

	return nlines;
} /* count_lines */

/**
 * Advances file pointer to next newline.
 *
 * @param fp file pointer
 * @return void
 */
void skip_line(FILE *fp)
{
	char a;

	/* skip first line */
	do {
		a = fgetc(fp);
		if (a == '\n')
			break;
	} while (!feof(fp));
} /* skip_line */

/**
 * Read next word from file stream.  This function allocates just enough
 * memory in its second argument to store the result.  Initial white space
 * is skipped and not included in the word.  The file pointer is left
 * pointing at the character immediately after the word.
 *
 * @param fp file handle
 * @param word unallocated pointer to char to store word
 */
int read_next_word(FILE *fp, char **word)
{
	char a;
	int len = 0;
	long int fpos;	/* file position */

	/* skip through white space */
	do {
		fpos = ftell(fp);
		a = fgetc(fp);
	} while (!feof(fp) && (a == '\n' || a == ' ' || a == '\t'));

	/* count length of word */
	do {
		len++;
		a = fgetc(fp);
	} while (!feof(fp) && a != '\n' && a != ' ' && a != '\t');

	/* rewind to beginning of word */
	fseek(fp, fpos, SEEK_SET);

	/* allocate just enough memory */
	*word = malloc((len + 1) * sizeof **word);
	if (!*word)
		return message(stderr, __FILE__, __func__, __LINE__,
			ERROR_MSG, MEMORY_ALLOCATION,
			"word of length %d\n", len);

	/* read word */
	if (fgets(*word, len+1, fp) == NULL) {
		free(*word);
		return FILE_FORMAT_ERROR;
	}

	return NO_ERROR;
} /* read_next_word */

/**
 * Add new string to a set.
 *
 * @param list current set of strings
 * @param len current size of set
 * @param str string to add
 * @param loc location of string in final set
 * @return error status
 */
int add_to_string_set(char ***list, int *len, char *str, int *loc)
{
	int i;
	char **nlist;

	for (i=0; i<*len; i++)
		if (!strcmp((*list)[i], str)) {
			*loc = i;
			return NO_ERROR;
		}

	if (*len) {
		nlist = realloc(*list, (*len + 1) * sizeof **list);
		if (!nlist)
			return message(stderr, __FILE__, __func__, __LINE__,
				ERROR_MSG, MEMORY_ALLOCATION, "list");
	} else {
		nlist = malloc((*len + 1) * sizeof *nlist);
		if (!nlist)
			return message(stderr, __FILE__, __func__, __LINE__,
				ERROR_MSG, MEMORY_ALLOCATION, "list");
	}
	*list = nlist;
	(*list)[*len] = malloc((strlen(str)+1) * sizeof *str);
	if (!(*list)[*len]) 
		return message(stderr, __FILE__, __func__, __LINE__,
			ERROR_MSG, MEMORY_ALLOCATION, "list[%d]", *len);
	strcpy((*list)[*len], str);
	*loc = *len;
	(*len)++;
	return NO_ERROR;
} /* add_to_string_set */

/**
 * Bubble sort.  Sort an array of integers and using the bubble sort algorithm.
 *
 * @param numbers array of numbers to be sorted
 * @param array_size size of the array that is to be sorted
 * @return void
 */
void bubbleSort(int *numbers, int array_size)
{
	int i, j, temp;

	for (i = (array_size - 1); i >= 0; i--)
		for (j = 1; j <= i; j++)
			if (numbers[j-1] > numbers[j]) {
				temp = numbers[j-1];
				numbers[j-1] = numbers[j];
				numbers[j] = temp;
			}
} /* End of bubbleSort(). */

int read_qfile(options *opt, data *dat, model *mod)
{
	int i, k;
	int err = NO_ERROR;
	FILE *fp = fopen(opt->qfile, "r");

	if (!fp) {
		err = message(stderr, __FILE__, __func__, __LINE__, ERROR_MSG,
			FILE_OPEN_ERROR, opt->qfile);
		return err;
	}

	if (opt->admixture && !opt->eta_constrained)
		for (i = 0; i < dat->I; i++) {
			for (k = 0; k < mod->K; k++)
#ifndef OLDWAY
				if (fscanf(fp, "%lf", &mod->vetaik[mod->tindex][i][k]) != 1) {
#else
				if (fscanf(fp, "%lf", &mod->etaik[i][k]) != 1) {
#endif
					err = message(stderr, __FILE__,
						__func__, __LINE__, ERROR_MSG,
						FILE_FORMAT_ERROR, opt->qfile);
					goto RETURN;
				}
		}
	else
		for (k = 0; k < mod->K; k++)
#ifndef OLDWAY
			if (fscanf(fp, "%lf", &mod->vetak[mod->tindex][k]) != 1) {
#else
			if (fscanf(fp, "%lf", &mod->etak[k]) != 1) {
#endif
				err = message(stderr, __FILE__, __func__,
					__LINE__, ERROR_MSG, FILE_FORMAT_ERROR,
					opt->qfile);
				goto RETURN;
			}

RETURN:
	fclose(fp);
	return err;
} /* read_qfile */

int read_pfile(options *opt, data *dat, model *mod)
{
	int l, k;
	int err = NO_ERROR;
	FILE *fp = fopen(opt->pfile, "r");

	if (!fp) {
		err = message(stderr, __FILE__, __func__, __LINE__, ERROR_MSG,
			FILE_OPEN_ERROR, opt->qfile);
		return err;
	}

	/* assumes biallelic locus */
	for (l = 0; l < dat->L; l++)
		for (k = 0; k < mod->K; k++) {
#ifndef OLDWAY
			if (fscanf(fp, "%lf", &mod->vpklm[mod->tindex][k][l][0]) != 1) {
#else
			if (fscanf(fp, "%lf", &mod->pKLM[k][l][0]) != 1) {
#endif
				err = message(stderr, __FILE__,
					__func__, __LINE__, ERROR_MSG,
					FILE_FORMAT_ERROR, opt->qfile);
				goto RETURN;
			}
#ifndef OLDWAY
			mod->vpklm[mod->tindex][k][l][1] = 1 - mod->vpklm[mod->tindex][k][l][0];
#else
			mod->pKLM[k][l][1] = 1 - mod->pKLM[k][l][0];
#endif
		}

RETURN:
	fclose(fp);
	return err;
} /* read_pfile */

/**
 * Read assignment file.  Assumes contiguous integers {1,2,...,K} are used so
 * the maximum is the number of clusters.
 *
 * @param opt options
 * @param dat data
 * @param mod model
 * @return error state
 */
int read_afile(options *opt, data *dat)
{
	int i, err = NO_ERROR;
	FILE *fp = fopen(opt->afile, "r");

	if (!fp) {
		err = message(stderr, __FILE__, __func__, __LINE__, ERROR_MSG,
			FILE_OPEN_ERROR, opt->afile);
		return err;
	}

	MAKE_VECTOR(opt->partition_from_file, dat->I);
	opt->pK = 0;

	for (i = 0; i < dat->I; i++) {
		if (fscanf(fp, "%d", &opt->partition_from_file[i]) != 1) {
			err = message(stderr, __FILE__, __func__, __LINE__,
				ERROR_MSG, FILE_FORMAT_ERROR, opt->afile);
			goto RETURN;
		}
		opt->partition_from_file[i]--;
		if (opt->partition_from_file[i] > opt->pK)
			opt->pK = opt->partition_from_file[i];
	}
	opt->pK++;

RETURN:
	fclose(fp);
	return err;
} /* read_afile */
