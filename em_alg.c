/**
 * @file em_alg.c
 * @author Wei-Chen Chen
 * @author Arun Sethuraman
 * @author Karin Dorman, kdorman@iastate.edu
 * @date Wed Dec  5 19:05:51 CST 2012
 *
 * This file contains the functions for running iterations of the EM algorithm
 * under the (ad)mixture model.  Public functions are em(), which runs the EM to
 * convergence, and em_step(), which runs a single step.
 *
 * TODO
 * - There is code to deal with numerical errors in E step (admixture), but it
 *   is commented out for time efficiency.  We may want to check for numerical
 *   problems somehow and trigger this code when appropriate.
 */

#include "multiclust.h"

//#define DEBUG
#ifdef DEBUG
int print_bits(size_t size, void *v);
#endif

double e_step_admixture_new(options *opt, data *dat, model *mod);
double e_step_admixture_orig(options *opt, data *dat, model *mod);
void m_step_admixture_new(options *opt, data *dat, model *mod);
void m_step_admixture_orig(options *opt, data *dat, model *mod);
double e_step_mixture(data *dat, model *mod);
void m_step_mixture(options *opt, data *dat, model *mod);
double scale_log_sum(double *v, int n, double max_v);
void internal_print_param(options *opt, data *dat, model *mod, double ***pklm, double **etaik, double *etak);
int stop_condition(options *opt, model *mod, double loglik);


#define e_step_admixture(A, B, C) e_step_admixture_orig(A, B, C)
#define m_step_admixture(A, B, C) m_step_admixture_orig(A, B, C)

/**
 * [TODO] Comment this function.
 *
 * This is all code for block relaxation using Zhou's algorithm.
 */
void em(options *opt, data *dat, model *mod)
{
	int i, stop = 0;

	/* no iteration required when _model::K is 1 */
	if (mod->K == 1) {
		em_step(opt, dat, mod);
		/* log likelihood calculation is one step late */
#ifndef OLDWAY
		mod->logL = log_likelihood(opt, dat, mod, mod->tindex);
#else
		mod->logL = log_likelihood(opt, dat, mod, 0);
#endif
		return;
	}

	/* run a few (stable) EM iterations to get close to maximum */
	while (mod->n_iter < opt->n_init_iter && !stop) {
		/* EM runs in place: no storage of iterates/increments */
		stop = em_step(opt, dat, mod);
	}

	/* at this point, mod->pindex contains previous iterate x^{(t)} */
	/* collect all but last secant conditions */
#ifndef OLDWAY
	for (i = 1; i < opt->q; i++) {
	        em_2_steps(mod, dat, opt);
		mod->pindex = mod->findex;
	}
	if (mod->converged)
		return;
#endif


	do {
		
		/* [KSD TODO: no block relax implemented] */
		if (!opt->accel_scheme)
			stop = em_step(opt, dat, mod); 
		else
			stop = accelerated_em_step(opt, dat, mod); 
//		if (mod->n_iter > 5) exit(0);
//exit(0);

	} while (!stop);
//	print_param(opt, dat, mod, mod->tindex);
} /* em */


/**
 * Check if iterations should stop.
 *
 * @param opt options object
 * @param mod model object
 * @param loglik most recent computed log likelihood; previous in mod->logL
 * @return non-zero to stop
 */
int stop(options *opt, model *mod, double loglik)
{
	mod->n_iter++;

	/* some fatal error checks */
	if (isnan(loglik)) {
		message(stderr, __FILE__, __func__, __LINE__, ERROR_MSG,
			CUSTOM_ERROR, "nan\n");
		exit(0);
	}

	mod->stopped = stop_condition(opt, mod, loglik);

	/* a small decrease in log likelihood is OK */
	if (loglik < mod->logL && !mod->stopped) {
		message(stderr, __FILE__, __func__, __LINE__, ERROR_MSG,
			CUSTOM_ERROR, "log likelihood decrease (%f < %f; %e)\n",
			loglik, mod->logL, (loglik - mod->logL)/loglik);
		exit(0);
	}


	if (opt->verbosity > MINIMAL)
		fprintf(stderr, "%4d (%s): %.2f (delta): %.5g\n",
			mod->n_iter,
			mod->accel_step?accel_method_abbreviations[opt->accel_scheme]:"EM",
			loglik, loglik - mod->logL);


	mod->accel_step = 0;
	mod->logL = loglik;

	return mod->stopped;
} /* stop */

int stop_condition(options *opt, model *mod, double loglik)
{
	mod->seconds_run = ((double) clock() - mod->start)/CLOCKS_PER_SEC;

	/* stop by hard constraints */
	if (opt->max_iter && mod->n_iter > opt->max_iter) {
		mod->iter_stop = 1;
		return 1;
	}

	if (opt->n_seconds && mod->seconds_run > opt->n_seconds) {
		mod->time_stop = 1;
		return 1;
	}
	
	return converged(opt, mod, loglik);
} /* stop_condition */

int converged(options *opt, model *mod, double loglik)
{
	int stop = 1;
	double abs_diff = 0;	/* absolute difference */
	double rel_diff = 0;	/* relative difference */

	if (opt->abs_error)
		abs_diff = fabs(loglik - mod->logL);
	if (opt->rel_error)
		rel_diff = abs_diff / fabs(mod->logL);

	/* any one of these triggers can make the iterations continue */
	if (opt->abs_error && abs_diff > opt->abs_error)
		stop &= 0;
	if (opt->rel_error && rel_diff > opt->rel_error)
		stop &= 0;
	if (stop)
		mod->converged = 1;
	return stop;
} /* converged */

/**
 * One iteration of EM.  This function does one iteration of the EM algorithm,
 * using previous parameter estimates in _model::findex and placing updated
 * parameter estimates in _model::tindex, which is in place by default and
 * if #define OLDWAY.
 *
 * @param opt options object
 * @param dat data object
 * @param mod mod object
 * @return converged status
 */
int em_step(options *opt, data *dat, model *mod)
{
	double ll;

	if (opt->admixture) {
		ll = e_step_admixture(opt, dat, mod);
		m_step_admixture(opt, dat, mod);
	} else {
		ll = e_step_mixture(dat, mod);
		m_step_mixture(opt, dat, mod);
	}
	return stop(opt, mod, ll);
} /* End of em_step(). */


/**
 * One iteration of EM followed by one E step.  The second E step is needed to
 * get the updated log likelihood following one iteration of EM.
 *
 * @param opt options object
 * @param dat data object
 * @param mod mod object
 * @return log likelihood from previous iteration
 */
double em_e_step(options *opt, data *dat, model *mod)
{
	double ll;

	if (opt->admixture) {
		e_step_admixture(opt, dat, mod);
		m_step_admixture(opt, dat, mod);
		ll = e_step_admixture(opt, dat, mod);
	} else {
		e_step_mixture(dat, mod);
		m_step_mixture(opt, dat, mod);
		ll = e_step_mixture(dat, mod);
	}
	return ll;
} /* End of em_e_step(). */

double e_step_admixture_new(options *opt, data *dat, model *mod)
{
	int i, k, l, a, m;
	double tmp;
	double loglik = 0;

	if (!dat->ila && make_ila(dat))
		return 0;

	for (i = 0; i < dat->I; i++) {
		for (l = 0; l < dat->L; l++) {
			for (a = 0; a < dat->ploidy; a++) {
				m = dat->ila[i][l][a];

				/* [TODO] [KSD] cannot handle missing alleles */
				if (m == MISSING) {
					message(stderr, __FILE__, __func__,
						__LINE__, ERROR_MSG,
						CUSTOM_ERROR,
						"e_step_admixture() does not "
						"handle missing data");
					exit(0);
				}
				tmp = 0;
				for (k = 0; k < mod->K; k++) {
					mod->diklm[i][k][l][m] =
						(opt->eta_constrained
#ifndef OLDWAY
						? mod->vetak[mod->findex][k]
						: mod->vetaik[mod->findex][i][k])
						* mod->vpklm[mod->findex][k][l][m];
#else
						? mod->etak[k] : mod->etaik[i][k])
						* mod->pKLM[k][l][m];
#endif
					tmp += mod->diklm[i][k][l][m];
				}
				for (k = 0; k < mod->K; k++)
					mod->diklm[i][k][l][m] /= tmp;
				loglik += log(tmp);
			}
		}
	}

	return loglik;
} /* e_step_admixture_new */


/**
 * E step.
 *
 * @param opt options object
 * @param dat data object
 * @param mod model object
 * @return previous step log likelihood
 */
double e_step_admixture_orig(options *opt, data *dat, model *mod)
{
	int debug = 0;
	int i, k, l, m, m_start;
	double ldilmk[mod->K];
/*	
	double scale, max_ldilmk;
*/
	double tmp, loglik = 0;
/*	
	double logeta[mod->K];
*/

	/* TIME COMPLEXITY (let T = A_1 + ... + A_L) */
	/* current version: I*K*T*(2*K+3) \propto I*K^2*T */
	/* faster version: I*K*L*M*(2*K + 2) \propto I*K^2*L*M and L*M>T*/
	/* compare for full EM cycle):		current	vs. faster
		5*I*K*T + 2*I*K^2*T + 2*I*K + 2*K*T	vs. 2*I*K^2*L*M + 4*I*K*L*M + 2*I*K + I*K*M + 2*K*T
		5*I*K*T + 2*I*K^2*T			vs. 2*I*K^2*L*M + 4*I*K*L*M + I*K*M
		5*T + 2*K*T				vs. 2*K*L*M + 5*L*M*(4/5 + 1/(4*L))
		T*(5 + 2*K)				vs. L*M*[2*K + 5*(4/5 + 1/(4*L))]
		Suppose T = alpha*L*M, then the fast version is faster if
		alpha 					> 4/5 + 1/(4*L)
		which is almost certainly true except for SNP data on diploids
	 */
/*
	if (opt->eta_constrained)
		for (k = 0; k < mod->K; k++)
#ifndef OLDWAY
			logeta[k] = log(mod->vetak[mod->findex][k]);
#else
			logeta[k] = log(mod->etak[k]);
#endif
*/
	for (i = 0; i < dat->I; i++) {
/*
		if (!opt->eta_constrained)
			for (k = 0; k < mod->K; k++)
#ifndef OLDWAY
				logeta[k] = log(mod->etaik[mod->findex][i][k]);
#else
				logeta[k] = log(mod->etaik[i][k]);
#endif
*/
		for (l = 0; l < dat->L; l++) {
			m_start = dat->L_alleles && dat->L_alleles[l][0] == MISSING;
			for (m = m_start; m < dat->uniquealleles[l]; m++) {
				if (dat->ILM[i][l][m] == 0) {
					for (k = 0; k < mod->K; k++)
						mod->diklm[i][k][l][m] = 0;
					continue;
				}
/*
				max_ldilmk = -INFINITY;
*/
				tmp = 0;
				if (debug)
					fprintf(stderr, "%d, %d, %d:", i, l, m);
				for (k = 0; k < mod->K; k++) {
					ldilmk[k] = (opt->eta_constrained
#ifndef OLDWAY
						? mod->vetak[mod->findex][k]
						: mod->vetaik[mod->findex][i][k])
						* mod->vpklm[mod->findex][k][l][m];
#else
						? mod->etak[k] : mod->etaik[i][k])
						* mod->pKLM[k][l][m];
#endif
					if (debug)
						fprintf(stderr, " %f*%f",
#ifndef OLDWAY
							mod->vpklm[mod->findex][k][l][m],
							mod->vetaik[mod->findex][i][k]);
#else
							mod->pKLM[k][l][m],
							mod->etaik[i][k]);
#endif
/*
					ldilmk[k] = logeta[k]
#ifndef OLDWAY
						+ log(mod->vpklm[mod->findex][k][l][m]);
#else
						+ log(mod->pKLM[k][l][m]);
#endif
					if (ldilmk[k] > max_ldilmk)
						max_ldilmk = ldilmk[k];
					tmp += exp(ldilmk[k]);
*/
					tmp += ldilmk[k];
				}
				if (debug)
					fprintf(stderr, " = %f\n", tmp);
/*
				scale = scale_log_sum(ldilmk, mod->K, max_ldilmk);
*/
				for (k = 0; k < mod->K; k++)
					mod->diklm[i][k][l][m]
						= dat->ILM[i][l][m]
						* ldilmk[k]/tmp;
/*
					mod->diklm[i][k][l][m]
						= dat->ILM[i][l][m]
						* exp(ldilmk[k])/tmp;
*/
				loglik += dat->ILM[i][l][m] * log(tmp);
/*
				+ scale;
*/
			}
			if ((!dat->L_alleles
				|| dat->L_alleles[l][0] == MISSING)
				&& dat->ILM[i][l][0] > 0) {
				for (m = m_start;
					m < dat->uniquealleles[l]; m++){
					for (k=0; k < mod->K; k++) {
						mod->diklm[i][k][l][m] += 
							(opt->eta_constrained
#ifndef OLDWAY
							? mod->vetak[mod->findex][k]
							: mod->vetaik[mod->findex][i][k])
							* mod->vpklm[mod->findex][k][l][m]
							* dat->ILM[i][l][0];
#else
							? mod->etak[k]
							: mod->etaik[i][k])
							* mod->pKLM[k][l][m]
							* dat->ILM[i][l][0];
#endif
					}}

				/* TODO: FASTER?
				for (m = 0; m < M; m++) {
					temp = 0.0;
					for (j = 0; j < K; j+=)
						temp += etaik[i][j]
							* KLM[j][l][X[i][l][m]];
					diklm[i][k][l][X[i][l][m]] = etaik[i][k]
						* KLM[k][l][X[i][l][m]] / temp;
				}
				(see M-Step for matching changes) */
			}
		}
	}

	/* ALTERNATIVE: E+M combined uses less memory as I->infinity and less time as I,L->infinity
	 * Let T = A_1+...+A_L and U = max_l A_l
	 * MEMORY USAGE: NEW vs. OLD
	 *	2*I*K + 2*K*T + U*(1+K)	vs. I*K*[1 + L*M + T/I]
	 *	2*I*K*[1+T/I] + U*(1+K)	vs. I*K*[1 + L*M + T/I]
	 *	I*K*[1 + T/I] + U*(1+K)	vs. I*K*L*M
	 *	K*T + U*(1+K)		vs. I*K*(L*M-1)
	 * TIME COMPLEXITY: NEW vs. OLD
	 *	I*(1 + L*(2*calloc + 5*K*M) + 6*K*T + K) + 2*K*T	vs. I*K*T*(2*K+3) + I*K*(2+T) + K*T*(I+2)
	 *	I + I*L*(2*calloc) + 5*I*L*K*M + 6*I*K*T + I*K + 2*K*T	vs. 5*I*K*T + 2*I*K + 2*I*T*K^2 + 2*K*T
	 *	I + I*L*(2*calloc) + 5*I*L*K*M + I*K*T			vs. I*K + 2*I*T*K^2
	 */
	/*
	double *denom;
	double *dilka;
	double sum, tmp;
	int itmp;
	for (i = 0; i < I; i++) {
		sum = 0;
		for (l = 0; l < L; l++) {
			denom = calloc(uniquealleles[l], sizeof *denom);
			dilka = calloc(K*uniquealleles[l], sizeof *dilka);
			for (k = 0; k < K; k++)
				for (m = 0; m < M; m++) {
					tmp = etaik[i][k] * pkla[k][l][X[i][l][m]]
					dilka[K*X[i][l][m] + k] += tmp;
					denom[X[i][l][m]] += tmp;
				}
			for (m = 0; m < uniquealleles[l]; m++)
				for (k = 0; k < K; k++) {
					itmp = K*m + k;
					dilka[itmp] /= denom[m];
					new_pkla[k][l][m] += dilka[itmp];
					next_etaik[i][k] += dilka[itmp];
					sum += dilka[itmp];
				}
		}
		for (k = 0; k < K; k++)
			next_etaik[i][k] /= sum;
	}
	for (k = 0; k < K; k++)
		for (l = 0; l < L; l++) {
			sum = 0;
			for (m = 0; m < uniquealleles[l]; m++)
				sum += new_pkla[k][l][m];
			for (m = 0; m < uniquealleles[l]; m++)
				new_pkla[k][l][m] /= sum;
		}
	*/

	return loglik;
} /* End of e_step_admixture_orig(). */


void m_step_admixture_new(options *opt, data *dat, model *mod)
{
	int i, m, l, k, a, nhaplotypes;
	double tmp;

	for (l = 0; l < dat->L; l++) {
		if (dat->L_alleles && dat->L_alleles[l][0] == MISSING) {
			message(stderr, __FILE__, __func__, __LINE__, ERROR_MSG,
				CUSTOM_ERROR, "m_step_admixture can not handle "
				"missing data");
			exit(0);
		}
		for (m = 0; m < dat->uniquealleles[l]; m++)
			for (k = 0; k < mod->K; k++)
#ifndef OLDWAY
				mod->vpklm[mod->tindex][k][l][m] = 0;
#else
				mod->pKLM[k][l][m] = 0;
#endif
		for (k = 0; k < mod->K; k++) {
			tmp = 0;
			for (i = 0; i < dat->I; i++)
				for (a = 0; a < dat->ploidy; a++) {
					m = dat->ila[i][l][a];
#ifndef OLDWAY
					mod->vpklm[mod->tindex][k][l][m] +=
#else
					mod->pKLM[k][l][m] +=
#endif
						mod->diklm[i][k][l][m];
					tmp += mod->diklm[i][k][l][m];
				}

			for (m = 0; m < dat->uniquealleles[l]; m++)
#ifndef OLDWAY
				mod->vpklm[mod->tindex][k][l][m] /= tmp;
#else
				mod->pKLM[k][l][m] /= tmp;
#endif
		}
	}

	if (!opt->eta_constrained) {
		nhaplotypes = dat->ploidy * dat->L;
		for (i = 0; i < dat->I; i++) {
			for (k = 0; k < mod->K; k++) {
#ifndef OLDWAY
				mod->vetaik[mod->tindex][i][k] = 0;
#else
				mod->etaik[i][k] = 0;
#endif
				for (l = 0; l < dat->L; l++)
					for (a = 0; a < dat->ploidy; a++) {
						m = dat->ila[i][l][a];
#ifndef OLDWAY
						mod->vetaik[mod->tindex][i][k] +=
#else
						mod->etaik[i][k] +=
#endif
							mod->diklm[i][k][l][m];
					}
#ifndef OLDWAY
				mod->vetaik[mod->tindex][i][k] /= nhaplotypes;
#else
				mod->etaik[i][k] /= nhaplotypes; 
#endif
			}
		}
	} else {
		nhaplotypes = dat->ploidy * dat->L * dat->I;
		for (k = 0; k < mod->K; k++)
#ifndef OLDWAY
			mod->vetak[mod->tindex][k] = 0;
#else
			mod->etak[k] = 0;
#endif
		for (i = 0; i < dat->I; i++)
			for (a = 0; a < dat->ploidy; a++) {
				m = dat->ila[i][l][a];
#ifndef OLDWAY
				mod->vetak[mod->tindex][k] += mod->diklm[i][k][l][m];
#else
				mod->etak[k] += mod->diklm[i][k][l][m];
#endif
			}
		for (k = 0; k < mod->K; k++)
#ifndef OLDWAY
			mod->vetak[mod->tindex][k] /= nhaplotypes;
#else
			mod->etak[k] /= nhaplotypes;
#endif
	}
} /* m_step_admixture_new */


/**
 * M step for admixture model, without accounting for missing data.
 *
 * @param opt options object
 * @param dat data object
 * @param mod model object
 * @return void
 */
void m_step_admixture_orig(options *opt, data *dat, model *mod)
{
	int debug = 0;
	int i, k, l, m, m_start;
	double temp;

	/* TIME COMPLEXITY: 2*I*K*(1 + T)         + 2*K*T */
	/* FASTER VERSION:  2*I*K*(1 + L*M + M/2) + 2*K*T */
	/* T <= L*M with equality if same number of alleles per site */
	/* M << T if L is large */

	/* estimate new eta_k or eta_{ik} */
	if (!opt->admixture || opt->eta_constrained) {
		temp = 0.0;
		for (k = 0; k < mod->K; k++) {
#ifndef OLDWAY
			mod->vetak[mod->tindex][k] = 0; //opt->eta_lower_bound;
#else
			mod->etak[k] = 0; //opt->eta_lower_bound;
#endif
			for (l = 0; l < dat->L; l++) {
				for (m = dat->L_alleles && dat->L_alleles[l][0] == MISSING;
					m < dat->uniquealleles[l];
					m++)
					for (i = 0; i < dat->I; i++)
#ifndef OLDWAY
						mod->vetak[mod->tindex][k] +=
#else
						mod->etak[k] +=
#endif
							mod->diklm[i][k][l][m];
			}
#ifndef OLDWAY
			temp += mod->vetak[mod->tindex][k];
#else
			temp += mod->etak[k];
#endif
		}
		if (debug)
			fprintf(stderr, "etak:");
		for (k = 0; k < mod->K; k++) {
#ifndef OLDWAY
			mod->vetak[mod->tindex][k] /= temp;
#else
			mod->etak[k] /= temp;
#endif
			if (debug)
#ifndef OLDWAY
				fprintf(stderr, " %f", mod->vetak[mod->tindex][k]);
#else
				fprintf(stderr, " %f", mod->etak[k]);
#endif
		}
		if (debug)
			fprintf(stderr, "\n");
		if (opt->do_projection)
			simplex_project_eta(mod, opt, 0);
	} else {
		for (i = 0; i < dat->I; i++) {
			temp = 0.0;
			for (k = 0; k < mod->K; k++) {
#ifndef OLDWAY
				mod->vetaik[mod->tindex][i][k] = 0;//opt->eta_lower_bound;
#else
				mod->etaik[i][k] = 0;//opt->eta_lower_bound;
#endif
				for (l = 0; l < dat->L; l++) {
					for (m = dat->L_alleles && dat->L_alleles[l][0] == MISSING;
						m < dat->uniquealleles[l]; m++) {
#ifndef OLDWAY
						mod->vetaik[mod->tindex][i][k]
#else
						mod->etaik[i][k]
#endif
							+= mod->diklm[i][k][l][m];
						if (debug & 1)
							fprintf(stderr, "diklm[%d][%d][%d][%d]: %f (%f)\n",
								i, k, l, m, mod->diklm[i][k][l][m],
#ifndef OLDWAY
								mod->vetaik[mod->tindex][i][k]);
#else
								mod->etaik[i][k]);
#endif
					}
				}
#ifndef OLDWAY
				temp += mod->vetaik[mod->tindex][i][k];
#else
				temp += mod->etaik[i][k];
#endif
			}
			if (debug & 1)
				fprintf(stderr, "etaik[%d]:", i);
			for (k = 0; k < mod->K; k++) {
#ifndef OLDWAY
				mod->vetaik[mod->tindex][i][k] /= temp;
#else
				mod->etaik[i][k] /= temp;
#endif
				if (debug & 1)
#ifndef OLDWAY
					fprintf(stderr, " %f", mod->vetaik[mod->tindex][i][k]);
#else
					fprintf(stderr, " %f", mod->etaik[i][k]);
#endif
			}
			if (debug & 1)
				fprintf(stderr, "\n");
			if (opt->do_projection)
				simplex_project_eta(mod, opt, i);
		}
	}

	// estimate new p_KLMs 
	for (k = 0; k < mod->K; k++) {
		for (l = 0; l < dat->L; l++) {
			temp = 0.0;
			m_start = dat->L_alleles
				&& dat->L_alleles[l][0] == MISSING;
			for (m = m_start; m < dat->uniquealleles[l]; m++) {
#ifndef OLDWAY
				mod->vpklm[mod->tindex][k][l][m] = 0;//opt->p_lower_bound;
#else
				mod->pKLM[k][l][m] = 0;//opt->p_lower_bound;
#endif
				for (i = 0; i < dat->I; i++) {
					if (debug & 4)
						fprintf(stderr, "diklm[%d][%d][%d][%d]: %f\n", i, k, l, m, mod->diklm[i][k][l][m]);
#ifndef OLDWAY
					mod->vpklm[mod->tindex][k][l][m] += mod->diklm[i][k][l][m];
#else
					mod->pKLM[k][l][m] += mod->diklm[i][k][l][m];
#endif
				}
#ifndef OLDWAY
				temp += mod->vpklm[mod->tindex][k][l][m];
#else
				temp += mod->pKLM[k][l][m];
#endif
			}
			if (debug & 2)
				fprintf(stderr, "pKLM[%d][%d]:", k, l);
			for (m = m_start; m < dat->uniquealleles[l]; m++) {
#ifndef OLDWAY
				mod->vpklm[mod->tindex][k][l][m] /= temp; 
#else
				mod->pKLM[k][l][m] /= temp; 
#endif
				if (debug & 2)
#ifndef OLDWAY
					fprintf(stderr, " %f", mod->vpklm[mod->tindex][k][l][m]);
#else
					fprintf(stderr, " %f", mod->pKLM[k][l][m]);
#endif
			}
			if (debug & 2)
				fprintf(stderr, "\n");
			if (opt->do_projection)
				simplex_project_pklm(mod, dat, opt, k, l);
		}
	}
	//print_param(opt, dat, mod, mod->tindex);
} /* End of m_step_admixture_orig */
	
/**
 * E step for mixture model.
 *
 * @param dat data object
 * @param mod model object
 * @return previous step log likelihood
 */
double e_step_mixture(data *dat, model *mod)
{
	int i, l, m, k;
	double log_etak[mod->K];
	double temp;
	double max_ll, loglik = 0;

	for (k = 0; k < mod->K; k++)
#ifndef OLDWAY
		log_etak[k] = log(mod->vetak[mod->findex][k]);
#else
		log_etak[k] = log(mod->etak[k]);
#endif

	for (i = 0; i < dat->I; i++) {
		max_ll = -INFINITY;

		/* numerators */
		for (k = 0; k < mod->K; k++) {
			mod->vik[i][k] = log_etak[k];
#ifdef DEBUG
if (!isfinite(mod->vik[i][k])) {
	mmessage(DEBUG_MSG, NO_ERROR, "log(vik[%d,%d]) is not finite after initialization!\n", i, k);
	exit(0);
}
#endif
			for (l = 0; l < dat->L; l++) {
				for (m = dat->L_alleles && dat->L_alleles[l][0] == MISSING;
					m < dat->uniquealleles[l]; m++) {
					if (dat->ILM[i][l][m] == 0 ||
#ifndef OLDWAY
						mod->vpklm[mod->findex][k][l][m] == 0.0
#else
						mod->pKLM[k][l][m] == 0.0
#endif
					)
						continue;
					mod->vik[i][k] += dat->ILM[i][l][m]
#ifndef OLDWAY
						* log(mod->vpklm[mod->findex][k][l][m]);
#else
						* log(mod->pKLM[k][l][m]);
#endif
#ifdef DEBUG
if (!isfinite(mod->vik[i][k])) {
	mmessage(DEBUG_MSG, NO_ERROR, "log(vik[%d, %d]) is not finite after locus %d, allele %d (ILM=%d, pklm=%f)!\n", i, k, l, m, dat->ILM[i][l][m], mod->vpklm[mod->findex][k][l][m]);
	exit(0);
}
#endif
				}
			}
			if (mod->vik[i][k] > max_ll)
				max_ll = mod->vik[i][k];
		}

		/* normalize (possibly scaling) */
/*
		double scale = scale_log_sum(mod->vik[i], mod->K, max_ll);
 */
		temp = 0;
		for(k = 0; k < mod->K; k++) {
#ifdef DEBUG
fprintf(stderr, "i=%d, k=%d: (vik=%f; %d %d)", i, k, mod->vik[i][k], isnan(mod->vik[i][k]), mod->vik[i][k] == mod->vik[i][k]);
#endif
			mod->vik[i][k] = exp(mod->vik[i][k] - max_ll);
			temp += mod->vik[i][k];
#ifdef DEBUG
fprintf(stderr, " (vik=%e; %d %d) (temp=%e; %d) (max_ll=%f; %d)\n", mod->vik[i][k], isnan(mod->vik[i][k]), mod->vik[i][k] == mod->vik[i][k], temp, isnan(temp), max_ll, isnan(max_ll));
fprintf(stderr, "mod->vik[%d][%d]: ", i, k);
if (print_bits(sizeof(mod->vik[i][k]), (void *) &mod->vik[i][k]))
	exit(0);
fprintf(stderr, "temp: ");
if (print_bits(sizeof(temp), (void *) &temp))
	exit(0);
fprintf(stderr, "max_ll: ");
if (print_bits(sizeof(max_ll), (void *) &max_ll))
	exit(0);
if (isnan(temp) || isnan(mod->vik[i][k]) || isnan(max_ll))
	exit(0);
#endif
		}
#ifdef DEBUG
if (isnan(temp) || temp == 0.0) {
	mmessage(DEBUG_MSG, NO_ERROR, "temp is not finite for individual %d!\n", i);
	exit(0);
}
#endif
		for(k = 0; k < mod->K; k++)
			mod->vik[i][k] /= temp;

		/* restore scaling to log likelihood */
		loglik += log(temp) + max_ll;
#ifdef DEBUG
fprintf(stderr, "%d: %f + %f = %f\n", i, log(temp), max_ll, loglik);
if (isnan(loglik)) {
	mmessage(DEBUG_MSG, NO_ERROR, "loglik is not finite for individual %d (temp=%e, max_ll=%f)!\n", i, temp, max_ll);
	exit(0);
}
#endif
	}

	return loglik;
} /* End of e_step_mixture(). */

/**
 * M step for admixture model, without accounting for missing data.
 *
 * @param opt options object
 * @param dat data object
 * @param mod model object
 * @return void
 */
void m_step_mixture(options *opt, data *dat, model *mod)
{
	int i, l, m, k, m_start;
	double temp, maxl = -INFINITY;

	/* estimate eta */
	temp = 0.0;
	for (k = 0; k < mod->K; k++) {
#ifndef OLDWAY
		mod->vetak[mod->tindex][k] = 0;//opt->eta_lower_bound;
#else
		mod->etak[k] = 0;//opt->eta_lower_bound;
#endif
		for (i = 0; i < dat->I; i++) {
#ifndef OLDWAY
			mod->vetak[mod->tindex][k] += mod->vik[i][k];
#else
			mod->etak[k] += mod->vik[i][k];
#endif
		}
#ifndef OLDWAY
		temp += mod->vetak[mod->tindex][k];
		if (maxl < mod->vetak[mod->tindex][k])
			maxl = mod->vetak[mod->tindex][k];
#else
		temp += mod->etak[k];
		if (maxl < mod->etak[k])
			maxl = mod->etak[k];
#endif
#ifdef DEBUG
		if (isnan(temp)) {
			mmessage(DEBUG_MSG, NO_ERROR, "m_step_mixture, computing eta[%d]: temp=%f\n", k, temp);
			exit(0);
		}
#endif
	}
	for (k = 0; k < mod->K; k++) {
#ifndef OLDWAY
		mod->vetak[mod->tindex][k] /= temp;
#ifdef DEBUG
		fprintf(stderr, "mod->vetak[%d][%d] = %e (temp = %e; maxl = %e)\n", mod->tindex, k, mod->vetak[mod->tindex][k], temp, maxl);
#endif
#else
		mod->etak[k] /= temp;
#endif
	}
	if (opt->do_projection)
		simplex_project_eta(mod, opt, 0);

	/* estimate p_{klm} */
	for (k = 0; k < mod->K; k++)
		for(l = 0; l < dat->L; l++) {
			m_start = dat->L_alleles
				&& dat->L_alleles[l][0] == MISSING;
			temp = 0.0;
			for (m = m_start; m < dat->uniquealleles[l]; m++) {
#ifndef OLDWAY
				mod->vpklm[mod->tindex][k][l][m] = opt->p_lower_bound;
#else
				mod->pKLM[k][l][m] = opt->p_lower_bound;
#endif
				for (i = 0; i < dat->I; i++) {
					if (dat->ILM[i][l][m])
#ifndef OLDWAY
						mod->vpklm[mod->tindex][k][l][m]
#else
						mod->pKLM[k][l][m]
#endif
							+= mod->vik[i][k]
							* dat->ILM[i][l][m];
				}
#ifndef OLDWAY
				temp += mod->vpklm[mod->tindex][k][l][m];
#else
				temp += mod->pKLM[k][l][m];
#endif
#ifdef DEBUG
				if (isnan(temp)) {
					mmessage(DEBUG_MSG, NO_ERROR, "m_step_mixture, compute pklm[%d][%d][%d]: temp=%f\n", k, l, m, temp);
					exit(0);
				}
#endif
			}
			for (m = m_start; m < dat->uniquealleles[l]; m++)
#ifndef OLDWAY
				mod->vpklm[mod->tindex][k][l][m] /= temp;
#else
				mod->pKLM[k][l][m] /= temp;
#endif
			if (opt->do_projection)
				simplex_project_pklm(mod, dat, opt, k, l);
		}
} /* End of m_step_mixture(). */


/**
 * Scale vector elements.  In the EM algorithm there comes a time when
 * probabilities need to be computed by renormalizing values computed on an
 * arbitrary scale.  During normalization, we must return to the exp scale and
 * consequently there is the risk of numerical underflow (rounding to 0 of very
 * small numbers) or overflow (rounding to the biggest representable number for
 * very large numbers).  In the first case, if all numbers round to 0, then
 * there will be a divide by zero error.  In the second case, if the largest
 * number is truncated, then it will renormalize to 1.  To alleviate the
 * problem, we can translate the logged numbers so they fall within range when
 * exponentiated.  The resulting scaling factor on the non-logged scale cancels
 * during normalization.  This code finds an appropriate translation and
 * applies it so that subsequent normalization on the non-log scale will be
 * stable.
 *
 * @param v vector of logged values to normalize
 * @param n length of vector
 * @param max_v largest element of vector v
 * @return scaling parameter used
 */
double scale_log_sum(double *v, int n, double max_v)
{
	double max_ev = exp(max_v);
	double scale = 0;
	int k;

	if (max_ev == 0.0 || max_ev == HUGE_VAL) {
		/* find scaling factor that brings max_ev in range */
		scale = (max_ev == HUGE_VAL) ? max_v : -max_v;
		do {
			scale *= 0.5;
			max_ev = exp(scale);
		} while (max_ev == HUGE_VAL);

		/* diff on log scale that puts exp() in range */
		scale = max_v - scale;

		/* apply scale to vector */
		for (k=0; k<n; k++)
			v[k] -= scale;
	}
	return scale;
} /* scale_log_sum */

/**
 * Take three EM steps and store parameter estimates.  Accelerated steps
 * use secant conditions to improve the update.  These secant conditions
 * require storage of previous three iterations or two sets of differences.
 * This function takes these steps and stores the parameters/differences.
 *
 * In the process of updating...
 *
 * @param mod model object
 * @param dat data object 
 * @param opt options object
 * @param n number of iterations
 * @return double
 */
int em_2_steps(model *mod, data *dat, options *opt)
{
	int debug = 0;

#ifndef OLDWAY
	int i, k, l, m, j;

	if (debug>1) print_param(opt, dat, mod, mod->pindex);
	if (debug) fprintf(stderr, "%s: filling %d\n", __func__, mod->delta_index);

	/* the playground where we take the iterates */
	mod->findex = mod->pindex;
	mod->tindex = (mod->findex + 1) % 3;

	for (j = 0; j < 2; j++) {
		/* read _model::findex, update _model::tindex */
		if (em_step(opt, dat, mod))
			return 1;
		if (debug>1) print_param(opt, dat, mod, mod->tindex);
		for (k = 0; k < mod->K; k++)
			for (l = 0; l < dat->L; l++) {
				for (m = dat->L_alleles && dat->L_alleles[l][0] == MISSING;
					m < dat->uniquealleles[l]; m++) {
					if (!j)
						mod->u_pklm[mod->delta_index][k][l][m]
							= mod->vpklm[mod->tindex][k][l][m]
							- mod->vpklm[mod->findex][k][l][m];
					else
						mod->v_pklm[mod->delta_index][k][l][m]
							= mod->vpklm[mod->tindex][k][l][m]
							- mod->vpklm[mod->findex][k][l][m];
				}
			}
		if (opt->admixture && !opt->eta_constrained)
			for (i = 0; i < dat->I; i++)
				for (k = 0; k < mod->K; k++) {
					if (!j)
						mod->u_etaik[mod->delta_index][i][k]
							= mod->vetaik[mod->tindex][i][k]
							- mod->vetaik[mod->findex][i][k];
					else
						mod->v_etaik[mod->delta_index][i][k]
							= mod->vetaik[mod->tindex][i][k]
							- mod->vetaik[mod->findex][i][k];
				}
		else
			for (k = 0; k < mod->K; k++) {
				if (!j)
					mod->u_etak[mod->delta_index][k]
						= mod->vetak[mod->tindex][k]
						- mod->vetak[mod->findex][k];
				else
					mod->v_etak[mod->delta_index][k]
						= mod->vetak[mod->tindex][k]
						- mod->vetak[mod->findex][k];
			}

		/* update playground indices */
		mod->findex = mod->tindex;
		mod->tindex = (mod->findex + 1) % 3;
		if (mod->tindex == mod->pindex)
			mod->tindex = (mod->tindex + 1) % 3;
	}

	/* prepare for next secant condition */
	mod->delta_index = (mod->delta_index + 1) % opt->q;

#else
	COPY_3JAGGED_ARRAY(mod->init_pKLM, mod->pKLM, dat->uniquealleles);
	if (opt->admixture && !opt->eta_constrained)
		COPY_2ARRAY(mod->init_etaik, mod->etaik, mod->K);
	else
		COPY_1ARRAY(mod->init_etak, mod->etak, mod->K);

	if (em_step(opt, dat, mod))
		return 1;
	COPY_3JAGGED_ARRAY(mod->iter1_pKLM, mod->pKLM, dat->uniquealleles);
	if (opt->admixture && !opt->eta_constrained)
		COPY_2ARRAY(mod->iter1_etaik, mod->etaik, mod->K);
	else
		COPY_1ARRAY(mod->iter1_etak, mod->etak, mod->K);
	
	if (debug>1) {
		print_param(opt, dat, mod, 1);
		print_param(opt, dat, mod, 2);
	}

	if (em_step(opt, dat, mod))
		return 1;

	COPY_3JAGGED_ARRAY(mod->iter2_pKLM, mod->pKLM, dat->uniquealleles);
	if (opt->admixture && !opt->eta_constrained)
		COPY_2ARRAY(mod->iter2_etaik, mod->etaik, mod->K);
	else
		COPY_1ARRAY(mod->iter2_etak, mod->etak, mod->K);

	if (debug)
		print_param(opt, dat, mod, 3);

#endif
	return 0;
} /* em_2_steps */

/**
 * Print current parameter estimates.
 *
 * @param mod model object
 * @param dat data object
 * @param opt options object
 * @return null
 */
void print_param(options *opt, data *dat, model *mod, int which)
{
	double **etaik = NULL;
	double *etak = NULL;
	double ***pklm = NULL;

#ifndef OLDWAY
	pklm = mod->vpklm[which];
	if (opt->admixture && !opt->eta_constrained)
		etaik = mod->vetaik[which];
	else
		etak = mod->vetak[which];
#else
	if (!which) {
		pklm = mod->pKLM;
		etaik = mod->etaik;
		etak = mod->etak;
	} else if (which == 1) {
		pklm = mod->init_pKLM;
		etaik = mod->init_etaik;
		etak = mod->init_etak;
	} else if (which == 2) {
		pklm = mod->iter1_pKLM;
		etaik = mod->iter1_etaik;
		etak = mod->iter1_etak;
	} else if (which == 3) {
		pklm = mod->iter2_pKLM;
		etaik = mod->iter2_etaik;
		etak = mod->iter2_etak;
	}
#endif
	internal_print_param(opt, dat, mod, pklm, etaik, etak);
} /* print_param */

#ifndef OLDWAY
void print_udiff(options *opt, data *dat, model *mod, int which)
{
	internal_print_param(opt, dat, mod, mod->u_pklm[which], 
		opt->admixture && !opt->eta_constrained ? mod->u_etaik[which] : NULL,
		opt->admixture && !opt->eta_constrained ? NULL : mod->u_etak[which]);
} /* print_udiff */

void print_vdiff(options *opt, data *dat, model *mod, int which)
{
	internal_print_param(opt, dat, mod, mod->v_pklm[which], 
		opt->admixture && !opt->eta_constrained ? mod->v_etaik[which] : NULL,
		opt->admixture && !opt->eta_constrained ? NULL : mod->v_etak[which]);
} /* print_vdiff */
#endif

#ifdef OLDWAY
int print_param_diff(options *opt, data *dat, model *mod, int first, int second)
{
	int i, k, l, m;
	double **etaik = NULL;
	double *etak = NULL;
	double ***pklm = NULL;
	double **etaik1 = NULL, **etaik2 = NULL;
	double *etak1 = NULL, *etak2 = NULL;
	double ***pklm1 = NULL, ***pklm2 = NULL;

	switch (first) {
		case 0:
			if (opt->admixture && !opt->eta_constrained)
				etaik1 = mod->init_etaik;
			else
				etak1 = mod->init_etak;
			pklm1 = mod->init_pKLM;
			break;
		case 1:
			if (opt->admixture && !opt->eta_constrained)
				etaik1 = mod->iter1_etaik;
			else
				etak1 = mod->iter1_etak;
			pklm1 = mod->iter1_pKLM;
			break;
		case 2:
			if (opt->admixture && !opt->eta_constrained)
				etaik1 = mod->iter2_etaik;
			else
				etak1 = mod->iter2_etak;
			pklm1 = mod->iter2_pKLM;
			break;
	}
	switch (second) {
		case 0:
			if (opt->admixture && !opt->eta_constrained)
				etaik2 = mod->init_etaik;
			else
				etak2 = mod->init_etak;
			pklm2 = mod->init_pKLM;
			break;
		case 1:
			if (opt->admixture && !opt->eta_constrained)
				etaik2 = mod->iter1_etaik;
			else
				etak2 = mod->iter1_etak;
			pklm2 = mod->iter1_pKLM;
			break;
		case 2:
			if (opt->admixture && !opt->eta_constrained)
				etaik2 = mod->iter2_etaik;
			else
				etak2 = mod->iter2_etak;
			pklm2 = mod->iter2_pKLM;
			break;
	}

	if (opt->admixture && !opt->eta_constrained)
		MAKE_2ARRAY(etaik, dat->I, mod->K);
	else
		MAKE_1ARRAY(etak, mod->K);
	MAKE_3JAGGED_ARRAY(pklm, mod->K, dat->L, dat->uniquealleles);

	if (opt->admixture && !opt->eta_constrained) {
		for (i = 0; i < dat->I; i++) {
			for (k = 0; k < mod->K; k++) {
				etaik[i][k] = etaik1[i][k] - etaik2[i][k];
			}
		}
	} else {
		for (k = 0; k < mod->K; k++)
			etak[k] = etak1[k] - etak2[k];
	}

	for (l = 0; l < dat->L; l++)
		for (k = 0; k < mod->K; k++)
			for (m = 0; m < dat->uniquealleles[l]; m++)
				pklm[k][l][m] = pklm1[k][l][m] - pklm2[k][l][m];
	internal_print_param(opt, dat, mod, pklm, etaik, etak);
	
	if (opt->admixture && !opt->eta_constrained)
		FREE_2ARRAY(etaik);
	else
		FREE_1ARRAY(etak);
	FREE_3ARRAY(pklm);
	return NO_ERROR;
} /* print_param_diff */
#endif


void internal_print_param(options *opt, data *dat, model *mod, double ***pklm, double **etaik, double *etak)
{
	int i, k, l, m;
	if (opt->admixture && !opt->eta_constrained) {
		for (i = 0; i < dat->I; i++) {
			fprintf(stderr, "etaik[%d]:", i);
			for (k = 0; k < mod->K; k++) {
				fprintf(stderr, " %f", etaik[i][k]);
			}
			if (!((i+1) % 4))
				fprintf(stderr, "\n");
			else
				fprintf(stderr, " ");
		}
	} else {
		fprintf(stderr, "etaik:");
		for (k = 0; k < mod->K; k++)
			fprintf(stderr, " %f", etak[k]);
		fprintf(stderr, "\n");
	}

	for (l = 0; l < dat->L; l++)
		for (k = 0; k < mod->K; k++) {
			fprintf(stderr, "Population %d, locus %d:", k, l);
			for (m = dat->L_alleles && dat->L_alleles[l][0] == MISSING;
				m < dat->uniquealleles[l]; m++) {
				fprintf(stderr, " %f", pklm[k][l][m]);
			}
			fprintf(stderr, "\n");
		}
} /* internal_print_param */

#ifdef DEBUG
int print_bits(size_t size, void *v)
{
        char *cv = (char *) v;
	int is_nan = 1;

        for (size_t i = size - 1; ; --i) {
                for (int j = 7; j >= 0; --j) {
                        printf("%d", cv[i] >> j & 1);
			if (i == size - 1 && !(cv[i] >> j & 1))
				is_nan = 0;
			else if (i == size - 2 && j == 7 && !(cv[i] >> j & 1))
				is_nan = 0;
		}
                if (i)
                        printf(" ");
		else
			break;
        }
        printf("\n");

	return is_nan;
} /* print_bits */
#endif
