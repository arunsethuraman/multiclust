/**
 * @file accel_em.c
 * @author Arun Sethuraman
 * @author Karin S. Dorman
 * @date Sun Jan  4 22:33:34 CST 2015
 * 
 * This file contains functions that will perform accelerated EM.
 */

#include "multiclust.h"

double step_size(options *opt, data *dat, model *mod);
double accelerated_update(options *opt, data *dat, model *mod, double s);
#ifndef OLDWAY
double qn_accelerated_update(options *opt, data *dat, model *mod);
#endif

#ifdef LAPACK
	/* LU decomposition of A */
	void dgetrf_(int *M, int *N, double *A, int *lda, int *IPIV, int *INFO);

	/* inverse of A given LU decomposition */
	void dgetri_(int *N, double *A, int *lda, int *IPIV, double *WORK,
		int *lwork, int *INFO);
#endif

/**
 * Accelerated update. TODO: add comments
 *
 * @param dat data object
 * @param mod mod object
 * @param opt options object
 * @return converged status
 */
int accelerated_em_step(options *opt, data *dat, model *mod)
{
	int n_adjust = 0;
	double emll, ll=0, s = 0;

	/* get another secant condition */
	em_2_steps(mod, dat, opt);

	/* converged during EM iterations; no need for acceleration */
	if (mod->stopped)
		return 1;
	
	/* mod->findex contains latest EM iterate */
	/* mod->tindex is ready to take aEM iterate */
	/* mod->pindex contains last (attempted) aEM iterate */

	/* mod->logL is penultimate EM log likelihood (one step late) */
#ifndef OLDWAY
	emll = log_likelihood(opt, dat, mod, mod->findex);
#else
	emll = log_likelihood(opt, dat, mod, 3);
#endif

	/* calculate step size */
	if (opt->accel_scheme <= QN) {
		s = step_size(opt, dat, mod);
		/* invalid step size: fall back on EM step */
		if (isnan(s) || isinf(s))
			goto EM_EXIT;
	}

	/* attempt acceleration step, possibly with backtracking */
	do {
		if (opt->accel_scheme <= QN)
			ll = accelerated_update(opt, dat, mod, s);
#ifndef OLDWAY
		else
			/* note: you can test this function by running it for QN q=1 */
			ll = qn_accelerated_update(opt, dat, mod);
#endif

		if (opt->adjust_step && ll < emll) {
			if (opt->verbosity > MINIMAL)
				fprintf(stderr, "after attempt %d (of %d) of accel ll is %f, EM ll is %f, step is %f\n", n_adjust, opt->adjust_step, ll, emll, s);
			s = (s - 1)/2;
		}

	} while (n_adjust++ < opt->adjust_step && ll < emll && s < -1);

	if (opt->verbosity > MINIMAL)
		fprintf(stderr, "accelerated_em_step (%s): aEM %f %s EM %f (step size: %f)\n", 
			ll > emll ? opt->accel_abbreviation : "EM",
			ll, ll > emll ? ">" : "<", emll, s);

	/* accelerated step improved log likelihood */
	if (ll > emll) {

		/* aEM step: current iterate is what was just computed in tindex */
#ifndef OLDWAY
		mod->pindex = mod->tindex;
#endif
		mod->accel_step = 1;
		return 0;
	}

EM_EXIT:
	/* fall back on EM iterate */
#ifdef OLDWAY
	COPY_3JAGGED_ARRAY(mod->pKLM, mod->iter2_pKLM, dat->uniquealleles);
	if (opt->admixture && !opt->eta_constrained)
		COPY_2ARRAY(mod->etaik, mod->iter2_etaik, mod->K);
	else
		COPY_1ARRAY(mod->etak, mod->iter2_etak, mod->K);
#else
	/* EM step: current iterate is what was just previously computed in findex */
	mod->pindex = mod->findex;
#endif

	return 0;
} /* accelerated_em_step */


/**
 * Compute step size for SQUAREM and QN1 acceleration methods.  All these
 * methods compute a step size to use in the update equation that uses three
 * regular iterates of EM, or more precisely the last two increments from these
 * three iterates.  These last two increments are stored in the _model::u_*
 * and _model::v_* arrays.  The next position to be written is stored in
 * _model::tindex.
 *
 * @param opt options object
 * @param dat data object
 * @param mod model object
 * @return computed step size
 */
double step_size(options *opt, data *dat, model *mod)
{
	int debug = 0;
	int i, l, k, m;
	int m_start;
	double utu = 0;
	double utvu = 0;
	double vutvu = 0;
	double s;

#ifndef OLDWAY

	if (opt->admixture && !opt->eta_constrained) {
		for (i = 0; i < dat->I; i++)
			for (k = 0; k < mod->K; k++) {
				utu += mod->u_etaik[mod->delta_index][i][k]
					* mod->u_etaik[mod->delta_index][i][k];
				utvu += mod->u_etaik[mod->delta_index][i][k]
					* (mod->v_etaik[mod->delta_index][i][k]
					- mod->u_etaik[mod->delta_index][i][k]);
				vutvu += (mod->v_etaik[mod->delta_index][i][k]
					- mod->u_etaik[mod->delta_index][i][k])
					* (mod->v_etaik[mod->delta_index][i][k]
					- mod->u_etaik[mod->delta_index][i][k]);
if (debug>1) fprintf(stderr, "(%d,%d): %13e %13e : %13e %13e %13e\n", i, k, mod->u_etaik[mod->delta_index][i][k], mod->v_etaik[mod->delta_index][i][k], utu, utvu, vutvu);
			}
	} else {
		for (k = 0; k < mod->K; k++) {
			utu += mod->u_etak[mod->delta_index][k]
				* mod->u_etak[mod->delta_index][k];
			utvu += mod->u_etak[mod->delta_index][k]
				* (mod->v_etak[mod->delta_index][k]
				- mod->u_etak[mod->delta_index][k]);
			vutvu += (mod->v_etak[mod->delta_index][k]
				- mod->u_etak[mod->delta_index][k])
				* (mod->v_etak[mod->delta_index][k]
				- mod->u_etak[mod->delta_index][k]);
		}
	}

	for (k = 0; k < mod->K; k++)
		for (l = 0; l < dat->L; l++) {
			m_start = dat->L_alleles && dat->L_alleles[l][0] == MISSING;
			for (m = m_start; m < dat->uniquealleles[l]; m++) {
				utu += mod->u_pklm[mod->delta_index][k][l][m]
					* mod->u_pklm[mod->delta_index][k][l][m];
				utvu += mod->u_pklm[mod->delta_index][k][l][m]
					* (mod->v_pklm[mod->delta_index][k][l][m]
					- mod->u_pklm[mod->delta_index][k][l][m]);
				vutvu += (mod->v_pklm[mod->delta_index][k][l][m]
					- mod->u_pklm[mod->delta_index][k][l][m])
					* (mod->v_pklm[mod->delta_index][k][l][m]
					- mod->u_pklm[mod->delta_index][k][l][m]);
			}
		}

#else
	if (opt->admixture && !opt->eta_constrained) {
		for (i = 0; i < dat->I; i++) {
			mod->current_i = i;
			initialize_etaiks(mod);
			for (k = 0; k < mod->K; k++) {	
				utu += mod->U[k] * mod->U[k];
				utvu += mod->U[k] * (mod->V[k] - mod->U[k]);
				vutvu += (mod->V[k] - mod->U[k]) * (mod->V[k] - mod->U[k]);
if (debug>1) fprintf(stderr, "(%d,%d): %13e %13e : %13e %13e %13e\n", i, k, mod->U[k], mod->V[k], utu, utvu, vutvu);
			}
		}
	} else {
		initialize_etaks(mod);
		for (k = 0; k < mod->K; k++) {	
			utu += mod->U[k] * mod->U[k];
			utvu += mod->U[k] * (mod->V[k] - mod->U[k]);
			vutvu += (mod->V[k] - mod->U[k]) * (mod->V[k] - mod->U[k]);
		}
	}

	for (k = 0; k < mod->K; k++) {
		for (l = 0; l < dat->L; l++) {
			mod->current_k = k;
			mod->current_l = l;
			initialize_pklas(mod, dat);
			m_start = dat->L_alleles
				&& dat->L_alleles[l][0] == MISSING;
			for (m = m_start; m < dat->uniquealleles[l]; m++) {
				utu += mod->U[m] * mod->U[m];
				utvu += mod->U[m] * (mod->V[m] - mod->U[m]);
				vutvu += (mod->V[m] - mod->U[m]) * (mod->V[m] - mod->U[m]);
			}
		}
	}
#endif
	if (opt->accel_scheme == SQS1)
		s = utu / utvu;
	else if (opt->accel_scheme == SQS2)
		s = utvu / vutvu;
	else if (opt->accel_scheme == SQS3) {
		if (sqrt(utu) < 1e-8)
			return NAN;
		s = -sqrt( utu / vutvu );
	} else if (opt->accel_scheme == QN)
		s = - utu / utvu;
	else s = -1;

	if (opt->accel_scheme < QN && s > -1)
		s = -1;

	if (debug)
		fprintf(stderr, "step_size: %f, %f, %f -> %f\n", utu, utvu, vutvu, s);

	return s;
} /* step_size */

#ifndef OLDWAY
/**
 * Quasi-Newton acceleration for q>1 as per Zhou2011.  This code is currently
 * not optimized in any sense.  In particular, most of the entries in matrix A
 * need not be recomputed each call; only the row and column corresponding to
 * the newest secant condition needs to be updated.  Also, this code only
 * implements the cases q=2 and q=3 using explicit formulae for the matrix
 * inverse.
 *
 * KSD TODO: I'm working on avoiding extra calculations, but it seems to be
 * broken.  See variable first_time, set now not to trigger anything special.
 *
 * @param dat data object
 * @param mod mod object
 * @param opt options object
 * @return log likelihood after acceleration step
 */
double qn_accelerated_update(options *opt, data *dat, model *mod)
{
	int debug = 0;
	int q1, q2;
	int i, j, k, l, m, n;
	int vindex = mod->delta_index ? mod->delta_index - 1 : opt->q - 1;
	int uindex = vindex ? vindex - 1 : opt->q - 1;
	int first_time = 1;
	double utu, utv, det, ll;

	if (debug)
		fprintf(stderr, "%s: %d\n", __func__, mod->delta_index);

/*
	if (opt->q > 3) {
		message(stderr, __FILE__, __func__, __LINE__, ERROR_MSG,
			INTERNAL_ERROR, "QN for q>2 not implemented");
		exit(0);
	}
	*/

	/* invert matrix */
	q1 = first_time ? mod->delta_index : vindex;	// index of oldest U
	j = first_time ? 0 : opt->q - 1;
	do {
		// index of oldest U
		q2 = first_time ? mod->delta_index : vindex;
		n = first_time ? 0 : opt->q - 1;
		do {
			utu = 0;
			utv = 0;
			for (k=0; k<mod->K; k++) {
				if (opt->admixture && !opt->eta_constrained)
					for (i=0; i<dat->I; i++) {
						utu += mod->u_etaik[q1][i][k] * mod->u_etaik[q2][i][k];
						utv += mod->u_etaik[q1][i][k] * mod->v_etaik[q2][i][k];
					}
				else {
					utu += mod->u_etak[q1][k] * mod->u_etak[q2][k];
					utv += mod->u_etak[q1][k] * mod->v_etak[q2][k];
				}
				for (l = 0; l < dat->L; l++) {
					for (m = dat->L_alleles && dat->L_alleles[l][0] == MISSING;
						m < dat->uniquealleles[l]; m++) {
						utu += mod->u_pklm[q1][k][l][m] * mod->u_pklm[q2][k][l][m];
						utv += mod->u_pklm[q1][k][l][m] * mod->v_pklm[q2][k][l][m];
					}
				}
			}
			mod->cutu[n] = utu;	/* yes: it overwrites repeatedly, but ends saving latest update */
			if (opt->q > 3)
				mod->Ainv[j*opt->q + n] = utu - utv;
			else
				mod->A[j*opt->q + n] = utu - utv;
			n++;
			if (debug)
				fprintf(stderr, "A[%d][%d] (%d, %d): %e - %e\n", j, n-1, q1, q2, utu, utv);
			q2 = (q2 + 1) % opt->q;
		} while (q2 != mod->delta_index);
		q1 = (q1 + 1) % opt->q;
		j++;
	} while (q1 != mod->delta_index);


	/* compute inverse */
	if (opt->q == 1) {
		mod->Ainv[0] = 1/mod->A[0];
		if (debug) fprintf(stderr, "step_size: %f, %f -> %f\n", mod->A[0], mod->cutu[0], -mod->Ainv[0]*mod->cutu[0]);
	} else if (opt->q == 2) {
		if (debug) fprintf(stderr, "   A: %e %e %e %e\n", mod->A[0], mod->A[1], mod->A[2], mod->A[3]);
		det = mod->A[0]*mod->A[3] - mod->A[1]*mod->A[2];
		mod->Ainv[0] = mod->A[3]/det;
		mod->Ainv[3] = mod->A[0]/det;
		mod->Ainv[1] = -mod->A[1]/det;
		mod->Ainv[2] = -mod->A[2]/det;
		if (debug)
			fprintf(stderr, "Ainv: %e %e %e %e (%e)\n",
				mod->Ainv[0], mod->Ainv[1],
				mod->Ainv[2], mod->Ainv[3], det);
	} else if (opt->q == 3) {
		det = mod->A[0]*(mod->A[4]*mod->A[8] - mod->A[5]*mod->A[7])
			- mod->A[1]*(mod->A[8]*mod->A[3] - mod->A[5]*mod->A[6])
			+ mod->A[2]*(mod->A[3]*mod->A[7] - mod->A[4]*mod->A[6]);
		mod->Ainv[0] = (mod->A[4]*mod->A[8] - mod->A[5]*mod->A[7])/det;
		mod->Ainv[1] = (mod->A[2]*mod->A[7] - mod->A[1]*mod->A[8])/det;
		mod->Ainv[2] = (mod->A[1]*mod->A[5] - mod->A[2]*mod->A[4])/det;
		mod->Ainv[3] = (mod->A[5]*mod->A[6] - mod->A[3]*mod->A[8])/det;
		mod->Ainv[4] = (mod->A[0]*mod->A[8] - mod->A[2]*mod->A[6])/det;
		mod->Ainv[5] = (mod->A[2]*mod->A[3] - mod->A[0]*mod->A[5])/det;
		mod->Ainv[6] = (mod->A[3]*mod->A[7] - mod->A[4]*mod->A[6])/det;
		mod->Ainv[7] = (mod->A[1]*mod->A[6] - mod->A[0]*mod->A[7])/det;
		mod->Ainv[8] = (mod->A[0]*mod->A[4] - mod->A[1]*mod->A[3])/det;
#ifdef LAPACK
	} else {
		int lwork = opt->q * opt->q;
		int info;
		dgetrf_(&opt->q, &opt->q, mod->Ainv, &opt->q, mod->ipiv, &info);
		dgetri_(&opt->q, mod->Ainv, &opt->q, mod->ipiv, mod->work, &lwork, &info);
#endif
	}


	for (k=0; k<mod->K; k++) {
		if (opt->admixture && !opt->eta_constrained)
			for (i=0; i<dat->I; i++) {
				mod->vetaik[mod->tindex][i][k] = mod->vetaik[mod->pindex][i][k] + mod->u_etaik[uindex][i][k];
			}
		else
			mod->vetak[mod->tindex][k] = mod->vetak[mod->pindex][k] + mod->u_etak[uindex][k];
		for (l = 0; l < dat->L; l++) {
			for (m = dat->L_alleles && dat->L_alleles[l][0] == MISSING;
				m < dat->uniquealleles[l]; m++)
				mod->vpklm[mod->tindex][k][l][m] = mod->vpklm[mod->pindex][k][l][m]
					+ mod->u_pklm[uindex][k][l][m];
			if (debug>1) fprintf(stderr, "pklm[%d][%d][0]: %f -> %f\n", k, l, mod->vpklm[mod->tindex][k][l][0], mod->vpklm[mod->pindex][k][l][0]);
		}
	}
	q1 = mod->delta_index;
	j = 0;
	do {
		n = 0;
		q2 = mod->delta_index;
		do {
			for (k=0; k<mod->K; k++) {
				if (opt->admixture && !opt->eta_constrained)
					for (i=0; i<dat->I; i++)
						mod->vetaik[mod->tindex][i][k] += mod->v_etaik[q1][i][k] * mod->Ainv[j*opt->q + n] * mod->cutu[n];
				else
					mod->vetak[mod->tindex][k] += mod->v_etak[q1][k] * mod->Ainv[j*opt->q + n] * mod->cutu[n];
				for (l = 0; l < dat->L; l++) {
					for (m = dat->L_alleles && dat->L_alleles[l][0] == MISSING;
						m < dat->uniquealleles[l]; m++)
						mod->vpklm[mod->tindex][k][l][m] += mod->v_pklm[q1][k][l][m] * mod->Ainv[j*opt->q + n] * mod->cutu[n];
					if (debug>1) fprintf(stderr, "pklm[%d][%d][0]: %f -> %f\n", k, l, mod->vpklm[mod->tindex][k][l][0], mod->vpklm[mod->pindex][k][l][0]);
				}
			}
			q2 = (q2 + 1) % opt->q;
			n++;
		} while (q2 != mod->delta_index);
		q1 = (q1 + 1) % opt->q;
		j++;
	} while (q1 != mod->delta_index);
	if (opt->admixture && !opt->eta_constrained)
		for (i = 0; i < dat->I; i++) {
			simplex_project_eta(mod, opt, i);
			if (debug>1) fprintf(stderr, "etaik[%d][0]: %f -> %f\n", i, mod->vetaik[mod->tindex][i][0], mod->vetaik[mod->pindex][i][0]);
		}
	else
		simplex_project_eta(mod, opt, 0);
	for (k = 0; k < mod->K; k++)
		for (l = 0; l < dat->L; l++)
			simplex_project_pklm(mod, dat, opt, k, l);

	ll = log_likelihood(opt, dat, mod, mod->tindex);
	return ll;
} /* qn_accelerated_update */
#endif

double accelerated_update(options *opt, data *dat, model *mod, double s)
{
	int debug = 0;
	double ll;
	int i, l, k, m, m_start;

	if (debug>1) {
#ifndef OLDWAY
	print_param(opt, dat, mod, mod->pindex);
	print_udiff(opt, dat, mod, mod->delta_index);
	print_vdiff(opt, dat, mod, mod->delta_index);
#else
	print_param(opt, dat, mod, 1);
	print_param_diff(opt, dat, mod, 1, 0);
	print_param_diff(opt, dat, mod, 2, 1);
#endif
	}
#ifndef OLDWAY
	mod->delta_index = mod->delta_index ? mod->delta_index - 1 : opt->q - 1;
#endif

	/* update allele frequencies */
	for (l = 0; l < dat->L; l++)
		for (k = 0; k < mod->K; k++) {
			m_start = dat->L_alleles
				&& dat->L_alleles[l][0] == MISSING;
			for (m=m_start; m < dat->uniquealleles[l]; m++) {
				if (opt->accel_scheme == QN)
#ifndef OLDWAY
					mod->vpklm[mod->tindex][k][l][m] = 
						mod->vpklm[mod->pindex][k][l][m]
						+ mod->u_pklm[mod->delta_index][k][l][m]
						+ s * mod->v_pklm[mod->delta_index][k][l][m];
#else
					mod->pKLM[k][l][m] = 
						(1-s) * mod->iter1_pKLM[k][l][m]
						+ s * mod->iter2_pKLM[k][l][m];
#endif
				else
#ifndef OLDWAY
					mod->vpklm[mod->tindex][k][l][m] = 
						mod->vpklm[mod->pindex][k][l][m]
						- 2*s*mod->u_pklm[mod->delta_index][k][l][m]
						+ s*s*(mod->v_pklm[mod->delta_index][k][l][m]
						- mod->u_pklm[mod->delta_index][k][l][m]);
#else
					mod->pKLM[k][l][m] = 
						mod->init_pKLM[k][l][m]
						- 2*s*(mod->iter1_pKLM[k][l][m] - mod->init_pKLM[k][l][m])
						+ s*s*(mod->iter2_pKLM[k][l][m] - 2*mod->iter1_pKLM[k][l][m] + mod->init_pKLM[k][l][m]);
#endif
			}
			simplex_project_pklm(mod, dat, opt, k, l);
#ifndef OLDWAY
			if (debug)
				fprintf(stderr, "pklm[%d][%d][0]: %f -> %f\n", k, l, mod->vpklm[mod->tindex][k][l][0], mod->vpklm[mod->pindex][k][l][0]);
#endif
		}
		
	/* update mixing proportions */
	if (opt->admixture && !opt->eta_constrained) {
		for (i = 0; i < dat->I; i++) {
			for (k = 0; k < mod->K; k++)
				if (opt->accel_scheme == QN)
#ifndef OLDWAY
					mod->vetaik[mod->tindex][i][k] =
						mod->vetaik[mod->pindex][i][k]
						+ mod->u_etaik[mod->delta_index][i][k]
						+ s * mod->v_etaik[mod->delta_index][i][k];
#else
					mod->etaik[i][k] =
						(1-s) * mod->iter1_etaik[i][k]
						+ s * mod->iter2_etaik[i][k];
#endif
				else
#ifndef OLDWAY
					mod->vetaik[mod->tindex][i][k] =
						mod->vetaik[mod->pindex][i][k]
						- 2 * s * mod->u_etaik[mod->delta_index][i][k]
						+ s * s * (mod->v_etaik[mod->delta_index][i][k]
						- mod->u_etaik[mod->delta_index][i][k]);
#else
					mod->etaik[i][k] =
						mod->init_etaik[i][k]
						- 2 * s * (mod->iter1_etaik[i][k] - mod->init_etaik[i][k])
						+ s * s * (mod->iter2_etaik[i][k] - 2*mod->iter1_etaik[i][k] + mod->init_etaik[i][k]);
#endif
			simplex_project_eta(mod, opt, i);
		}
	} else {
		for (k = 0; k < mod->K; k++)
			if (opt->accel_scheme == QN)
#ifndef OLDWAY
				mod->vetak[mod->tindex][k] =
					mod->vetak[mod->pindex][k]
					+ mod->u_etak[mod->delta_index][k]
					+ s * mod->v_etak[mod->delta_index][k];
#else
				mod->etak[k] =
					(1-s) * mod->iter1_etak[k]
					+ s * mod->iter2_etak[k];
#endif
			else
#ifndef OLDWAY
				mod->vetak[mod->tindex][k] =
					mod->vetak[mod->pindex][k]
					- 2 * s * mod->u_etak[mod->delta_index][k]
					+ s * s * (mod->v_etak[mod->delta_index][k]
					- mod->u_etak[mod->delta_index][k]);
#else
				mod->etak[k] =
					mod->init_etak[k]
					- 2 * s * (mod->iter1_etak[k] - mod->init_etak[k])
					+ s * s * (mod->iter2_etak[k] - 2*mod->iter1_etak[k] + mod->init_etak[k]);
#endif
		simplex_project_eta(mod, opt, -1);
	}
#ifndef OLDWAY
	//print_param(opt, dat, mod, mod->tindex);
	ll = log_likelihood(opt, dat, mod, mod->tindex);
	mod->delta_index = (mod->delta_index + 1) % opt->q;
#else
	//print_param(opt, dat, mod, 0);
	ll = log_likelihood(opt, dat, mod, 0);
#endif
	return ll;
} /* accelerated_udpate */

#ifdef OLDWAY
/**
 * Compute last two finite differences in pkla block.
 *
 * @param mod model object
 * @param dat data object
 * @return none
 */
void initialize_pklas(model *mod, data *dat)
{
	int k, l, m, m_start;

	l = mod->current_l;
	k = mod->current_k;
	m_start = dat->L_alleles && dat->L_alleles[l][0] == MISSING;

	for (m = m_start; m < dat->uniquealleles[l] ; m++) {
		/* init_ -> iter1_ -> iter2_ */
		mod->U[m] = mod->iter1_pKLM[k][l][m]
			- mod->init_pKLM[k][l][m];
		mod->V[m] = mod->iter2_pKLM[k][l][m]
			- mod->iter1_pKLM[k][l][m];
	}				
}/* End of initialize_pklas */


/**
 * Compute last two finite differences in etaik block.
 *
 * @param mod model object
 * @return none
 */
void initialize_etaiks(model *mod)
{
	int k;

	/* current individual */
	int i = mod->current_i;

	for (k = 0; k < mod->K ; k++) {
		mod->U[k] = mod->iter1_etaik[i][k]
			- mod->init_etaik[i][k];
		mod->V[k] = mod->iter2_etaik[i][k]
			- mod->iter1_etaik[i][k];
	}
} /* End of initialize_etaiks function */


/**
 * Compute last two finite differences in the etak block.
 *
 * @param mod model object
 * @return none
 */
void initialize_etaks(model *mod)
{
	int k;

	for (k = 0; k < mod->K; k++) {
		mod->U[k] = mod->iter1_etak[k]
			- mod->init_etak[k];
		mod->V[k] = mod->iter2_etak[k]
			- mod->iter1_etak[k];
	}
} /* End of initialize_etaks function */
#endif
