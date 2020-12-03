/**
 * @file simplex.c
 * @author Arun Sethuraman
 * @author Karin S. Dorman
 * @date Mon Jan  5 12:39:51 CST 2015
 * 
 * This file contains functions for simplex projections.
 */

#include "multiclust.h"
#define MAKE_1ARRAY MAKE_1ARRAY_RETURN /* return on memory allocation error */

/* projection onto the simplex */
void michelot_project(double *, int, double, double);
double chenye_project(double *, int);
void bubbleSort_doubles(double *numbers, int array_size);

/**
 * Simplex projection function. This function imposes the condition that all
 * the mixing proportions or the pklas for a particular individual, or an
 * allele m at a locus l should sum to 1. If not, it projects it to the nearest
 * positive point on the simplex.  
 *
 * @param mod model object
 * @param opt options object
 * @return none
 * TODO: use opt->lower_bound instead of hard-coded 1e-5
 */
int simplex_project_eta(model *mod, options *opt, int i)
{
	if (opt->admixture && !opt->eta_constrained)
#ifndef OLDWAY
		michelot_project(mod->vetaik[mod->tindex][i], mod->K, 1.0, opt->eta_lower_bound);
#else
		michelot_project(mod->etaik[i], mod->K, 1.0, opt->eta_lower_bound);
#endif
	else
#ifndef OLDWAY
		michelot_project(mod->vetak[mod->tindex], mod->K, 1.0, opt->eta_lower_bound);
#else
		michelot_project(mod->etak, mod->K, 1.0, opt->eta_lower_bound);
#endif

	return NO_ERROR;
} /* simplex_project_eta */

int simplex_project_pklm(model *mod, data *dat, options *opt, int k, int l)
{
	double *params;
	int m, m_start;

	m_start = dat->L_alleles && dat->L_alleles[l][0] == MISSING;
	MAKE_1ARRAY(params, dat->uniquealleles[l] - m_start);
	for (m = m_start; m < dat->uniquealleles[l]; m++)
#ifndef OLDWAY
		params[m - m_start] = mod->vpklm[mod->tindex][k][l][m];
#else
		params[m - m_start] = mod->pKLM[k][l][m];
#endif
	michelot_project(params, dat->uniquealleles[l] - m_start, 1.0, opt->p_lower_bound);
	for (m = m_start; m < dat->uniquealleles[l]; m++)
#ifndef OLDWAY
		mod->vpklm[mod->tindex][k][l][m] = params[m - m_start];
#else
		mod->pKLM[k][l][m] = params[m - m_start];
#endif
	FREE_1ARRAY(params);
	return NO_ERROR;
} /* simplex_project_pklm */

/**
 * Simplex projection algorithm. This method uses the algorithm proposed in
 * the following reference: Yunmei Chen, Xiaojing Ye. "Projection Onto A
 * Simplex" http://arxiv.org/abs/1101.6081
 *
 */
double chenye_project(double *params, int len)
{
	int index, j;
	double numer;
	int oloop_index = 0;
	double t,ti;
	double sum = 0;	/* WARNING: not properly set up (this code now broken) */

	bubbleSort_doubles(params, len);

	oloop_index = len;
	ti = 0.0;

	/* simplex projection method */
	t = (sum - 1) / oloop_index;
	for (index = oloop_index - 1; index > 0; index--) {
		numer = 0.0;
		for (j = index; j < oloop_index; j++)	/* [KSD] was BUG for pkla */
			numer += params[j];
		ti = (numer - 1)/(oloop_index - index);
		if (ti >= params[index - 1]) {
			t = ti;
			break;
		}/* else if (index < 0) {
			t = (sum - 1)/(oloop_index-1);
			break;
		}*/
	}

	return t;
} /* End of chenye_project function */

void michelot_project(double *params, int len, double sum, double min)
{
	int debug = 0;
	double csum;
	double aty_c_a2;
	int fixed[len];
	int i, n = len;
	int can_terminate;

	for (i = 0; i < len; i++)
		fixed[i] = 0;
	
	while (n) {
		csum = 0.0;
		for (i = 0; i < len; i++)
			csum += params[i];
		aty_c_a2 = (csum - sum) / n;
		can_terminate = 1;
		for (i = 0; i < len; i++)
			if (!fixed[i]) {
				params[i] -= aty_c_a2;
				if (params[i] < min) {
					if (debug)
						fprintf(stderr, "%s::%s(%d): resetting to minimum %e\n", __FILE__, __func__, __LINE__, min);
					params[i] = min;
					fixed[i] = 1;
					n--;
					can_terminate = 0;
				}
			}
		if (can_terminate)
			break;
	}

} /* michelot_project */

	
/**
 * Bubble sort.  Sort an array of doubles and using the bubble sort algorithm.
 *
 * @param numbers array of numbers to be sorted
 * @param array_size size of the array that is to be sorted
 * @return none
 */
void bubbleSort_doubles(double *numbers, int array_size)
{
	int i, j, temp;

	for (i = array_size - 1; i >= 0; i--)
		for (j = 1; j <= i; j++)
			if (numbers[j-1] > numbers[j]) {
				temp = numbers[j-1];
				numbers[j-1] = numbers[j];
				numbers[j] = temp;
			}
} /* End of bubbleSort_doubles(). */
