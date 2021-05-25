/**
 * @file log_likelihood.c
 * @author Arun Sethuraman
 * @author Karin Dorman
 * @date Mon Jan  5 12:17:44 CST 2015
 * 
 * This file contains functions that compute the log likelihood.
 */

#include <math.h>

#include "multiclust.h"
#define MAKE_1ARRAY MAKE_1ARRAY_RETURN /* return on memory allocation error */

/* log likelihood calculation and tests */
int loglikelihood_decrease(options *opt, double pll, double ll);
double logL_admixture(options *opt, data *dat, model *mod, int);
double logL_mixture(data *dat, model *mod, int);


/**
 * Determine if the log likelihood has decreased.  Very small decreases in log
 * likelihood could be caused by numerical error and tend to happen near
 * convergence, slowing down the algorithm by forcing a full-blown EM step.
 * This function checks for a decrease in log likelihood, but requires the
 * decrease to be of non-trivial size.  It is a bit tricky to decide what is
 * non-trivial.
 *
 * @param opt options object
 * @param pll previous log likelihood
 * @param ll current log likelihood
 * @return non-zero to indicate log likelihood has descreased
 */
int loglikelihood_decrease(options *opt, double pll, double ll)
{
	/* if abs_error is set, use it to determine a "true" decrease */
	if (pll > ll && opt->abs_error && pll - ll > opt->abs_error)
		return 1;
	/* else use something very small for double precision */
	else if (pll > ll && !opt->abs_error && pll - ll > 1e-15)
		return 1;

	return 0;
} /* loglikelihood_decrease */

				
/**
 * Compute log likelihood given current parameters.
 *
 * @param opt options object
 * @param dat data object
 * @param mod model object
 * @param which which parameters to use
 * @return log likelihood
 */
double log_likelihood(options *opt, data *dat, model *mod, int which)
{
	if (opt->admixture)	
		return logL_admixture(opt, dat, mod, which);	/* admixture model */
	else
		return logL_mixture(dat, mod, which);	/* mixtude model */
} /* log_likelihood */

/**
 * AIC
 *
 * @param mod model object
 * @return double
 */
double aic(model *mod)
{
	return(-2 * mod->max_logL + 2 * mod->no_parameters);
}/* aic */

/**
 * AIC
 *
 * @param dat data object
 * @param mod model object
 * @return double
 */
double bic(data *dat, model *mod)
{
	return(-2 * mod->max_logL + mod->no_parameters * log((double) dat->I));
} /* bic */

/**
 * Compute log likelihood given current parameter estimates given admixture model.
 * 
 * @param opt options object
 * @param dat data object
 * @param mod model object
 * @param which which parameters to use
 * @return double
 */
double logL_admixture(options *opt, data *dat, model *mod, int which)
{
	int i, k, l, m, m_start;
	double temp = 0.0, loglt1 = 0.0;
	double ***pKLM, **etaik, *etak;

#ifndef OLDWAY
	pKLM = mod->vpklm[which];
	if (!opt->eta_constrained)
		etaik = mod->vetaik[which];
	else
		etak = mod->vetak[which];
#else
	if (!which) {
		pKLM = mod->pKLM;
		etaik = mod->etaik;
		etak = mod->etak;
	} else if (which == 1) {
		pKLM = mod->init_pKLM;
		etaik = mod->init_etaik;
		etak = mod->init_etak;
	} else if (which == 2) {
		pKLM = mod->iter1_pKLM;
		etaik = mod->iter1_etaik;
		etak = mod->iter1_etak;
	} else {
		pKLM = mod->iter2_pKLM;
		etaik = mod->iter2_etaik;
		etak = mod->iter2_etak;
	}
#endif

	for (i = 0; i < dat->I; i++)
		for (l = 0; l < dat->L; l++) {
			m_start = dat->L_alleles
				&& dat->L_alleles[l][0] == MISSING ? 1 : 0;
			for (m = m_start; m < dat->uniquealleles[l]; m++) {
				if (dat->ILM[i][l][m] == 0)
					continue;

				temp = 0.0;
				for (k = 0; k < mod->K; k++)
					temp += (opt->eta_constrained
						? etak[k]
						: etaik[i][k])
						* pKLM[k][l][m];
				loglt1 += dat->ILM[i][l][m] * log(temp);
			}
		}

        return loglt1;
} /* End of logL_admixture(). */


/**
 * Compute log likelihood given current parameter estimates under mixture model.
 *
 * @param mod model object
 * @param dat data object
 * @return double 
 */
double logL_mixture(data *dat, model *mod, int which)
{
	int i, k, l, m, flag_out_range;
	double max_exp = 0.0, temp_exp = 0.0, scale_exp = 0.0, loglt1 = 0.0;
	double temp_vik[mod->K], log_etak[mod->K];
	double ***pKLM, *etak;

#ifndef OLDWAY
	pKLM = mod->vpklm[which];
	etak = mod->vetak[which];
#else
	if (!which) {
		pKLM = mod->pKLM;
		etak = mod->etak;
	} else if (which == 1) {
		pKLM = mod->init_pKLM;
		etak = mod->init_etak;
	} else if (which == 2) {
		pKLM = mod->iter1_pKLM;
		etak = mod->iter1_etak;
	} else {
		pKLM = mod->iter2_pKLM;
		etak = mod->iter2_etak;
	}
#endif

	
	for (k = 0; k < mod->K; k++)
		log_etak[k] = log(etak[k]);
	
	for (i = 0; i < dat->I; i++) {
		max_exp = -INFINITY;
		for (k = 0; k < mod->K; k++) {
			temp_vik[k] = 0.0;
			for (l = 0; l < dat->L; l++) {
				for (m = (dat->L_alleles &&
						dat->L_alleles[l][0] == MISSING
						? 1 : 0);
					m < dat->uniquealleles[l];
					m++) {
					if (dat->ILM[i][l][m] == 0)
						continue;
					temp_vik[k] += dat->ILM[i][l][m]
						* log(pKLM[k][l][m]);
				}
			}
			temp_vik[k] += log_etak[k];
			if (temp_vik[k] > max_exp)
				max_exp = temp_vik[k];
		}

		/* optionally scale term in log likelihood */
		temp_exp = exp(max_exp);
		scale_exp = 0.0;
		flag_out_range = 0.0;
		if (temp_exp == 0.0 || temp_exp == HUGE_VAL) {
			flag_out_range = 1;
			scale_exp = (temp_exp == HUGE_VAL) ? max_exp : -max_exp;
			do {
				scale_exp *= 0.5;
				temp_exp = exp(scale_exp);
			} while (temp_exp == HUGE_VAL);
			scale_exp = max_exp - scale_exp;
		}
	
		if (flag_out_range)
			for (k = 0; k < mod->K; k++)
				temp_vik[k] -= scale_exp;
		temp_exp = 0.0;
		for (k = 0; k < mod->K ; k++)
			temp_exp = temp_exp + exp(temp_vik[k]);
		loglt1 = loglt1 + log(temp_exp) + scale_exp;
	}

	return loglt1;
} /* End of logL_mixture() function */
