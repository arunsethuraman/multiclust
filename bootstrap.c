/**
 * @file bootstrap.c
 * @author Arun Sethuraman
 * @author Karin Dorman, kdorman@iastate.edu
 * @date Sun Dec  9 23:44:32 CST 2012
 *
 * Functions for parametric bootstrap.
 *
 * TODO
 * - need a serious runif() RNG instead of rand() in this file
 */

#include "multiclust.h"

void parametric_bootstrap_admixture(options *opt, data *dat, model *mod);
void parametric_bootstrap_mixture(data *dat, model *mod);

/**
 * Generate one parametric bootstrap data set.  This function allocates the
 * memory it needs to perform bootstrap.  Remember to call 
 * cleanup_parametric_bootstrap() when done bootstrapping to release this data.
 * If one were really cramped for space, it could be possible to bootstrap in
 * place, so long as the number of missing alleles were stored for each 
 * individual/locus.
 *
 * @param opt options object
 * @param dat data object
 * @param mod model object
 * @return error status
 */
int parametric_bootstrap(options *opt, data *dat, model *mod)
{

	/* allocate bootstrapped data, saving orig. data in dat->bs_ILM */
	if (!dat->bs_ILM) {
		int ***tmp_ptr = dat->ILM;
	        MAKE_3JAGGED_ARRAY(dat->bs_ILM, dat->I, dat->L,
			dat->uniquealleles);
		dat->ILM = dat->bs_ILM;
		dat->bs_ILM = tmp_ptr;
	}

	if (opt->admixture)
		parametric_bootstrap_admixture(opt, dat, mod);
	else
		parametric_bootstrap_mixture(dat, mod);
	
	/* CAUTION: cleanup_parametric_bootstrap() to release bootstrap data
	 * and restore original data */
	
	return NO_ERROR;
} /* End parameteric_bootstrap(). */

/**
 * Release data allocated for bootstrapping.
 *
 * @param dat data object
 * @return none
 */
void cleanup_parametric_bootstrap(data *dat)
{
	int ***tmp_ptr = dat->ILM;
	dat->ILM = dat->bs_ILM;
	FREE_3ARRAY(tmp_ptr);
	dat->bs_ILM = NULL;
} /* End cleanup_parametric_bootstrap() */

/**
 * Generates parametric bootstrap dataset under admixture model.
 *
 * @param opt options object
 * @param dat data object
 * @param mod model object
 * @return error status
 */
void parametric_bootstrap_admixture(options *opt, data *dat, model *mod)
{
	int i, j, l, m, n, n_end, m_start;
	double r, sum;

	for (i=0; i<dat->I; i++) {
		for (l=0; l<dat->L; l++) {
			/* missing data must stay missing */
			dat->ILM[i][l][0] = dat->bs_ILM[i][l][0];

			/* initialize counts of observed alleles */
			m_start = dat->L_alleles
				&& dat->L_alleles[l][0] == MISSING ? 1 : 0;
			for (m=m_start; m<dat->uniquealleles[l]; m++)
				dat->ILM[i][l][m] = 0;

			/* number of observed alleles */
			n_end = dat->ploidy - dat->ILM[i][l][0];

			for (n=0; n<n_end; n++) {
				/* choose subpopulation source: store in j */
				r = (double) rand() / RAND_MAX;
				j = 0;
				sum = 0;
				while (j < mod->K && r > sum)
					if (opt->eta_constrained)
						sum += mod->mle_etak[j++];
					else
						sum += mod->mle_etaik[i][j++];
				if (j)
					j--;

				/* choose allele: store in m */
				r = (double) rand() / RAND_MAX;
				m = dat->L_alleles
					&& dat->L_alleles[l][0] == MISSING ? 1 : 0;
				sum = 0;
				while (m < dat->uniquealleles[l] && r > sum)
					sum += mod->mle_pKLM[j][l][m++];
				if (m)
					m--;

				/* add chosen allele to data */
				dat->ILM[i][l][m]++;
			}
		}
	}
} /* End parametric_bootstrap_admixture(). */

/**
 * Generates parametric bootstrap dataset under mixture model.
 *
 * @param dat data object
 * @param mod model object
 * @return error status
 */
void parametric_bootstrap_mixture(data *dat, model *mod)
{
	int i, k, l, m, n, n_end, m_start;
	double r, sum;

	/* for each individual */
	for (i=0; i<dat->I; i++) {
		/* choose subpopulation source of individual */
		r = (double) rand() / RAND_MAX;	/* TODO need a serious rng */
		k = 0;
		sum = 0;
		while (k < mod->K && r > sum)
			sum += mod->mle_etak[k++];
		if (k)
			k--;

		/* for each locus */
		for (l=0; l<dat->L; l++) {
			/* missing data (cannot create data where was none) */
			dat->ILM[i][l][0] = dat->bs_ILM[i][l][0];

			/* initialize allele counts */
			m_start = dat->L_alleles
				&& dat->L_alleles[l][0] == MISSING ? 1 : 0;
			for (m=m_start; m<dat->uniquealleles[l]; m++)
				dat->ILM[i][l][m] = 0;

			/* for all non-missing observations, draw an allele */
			n_end = dat->ploidy - dat->ILM[i][l][0];
			for (n=0; n<n_end; n++) {
				r = (double) rand() / RAND_MAX;
				m = dat->L_alleles
					&& dat->L_alleles[l][0] == MISSING
					? 1 : 0;
				sum = 0;
				while (m < dat->uniquealleles[l] && r > sum)
					sum += mod->mle_pKLM[k][l][m++];
				if (m)
					m--;
				dat->ILM[i][l][m]++;
			}
		}
	}
} /* End parametric_bootstrap_mixture(). */
