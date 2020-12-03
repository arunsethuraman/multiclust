/**
 * @file rnd_init.c
 * @author Arun Sethuraman
 * @author Karin Dorman, kdorman@iastate.edu
 * @date Tue Dec  4 22:28:19 CST 2012
 *
 * Functions to initialize the model.  At conclusion, initial parameter values
 * have been chosen and placed in _model::vpklm[0], _model::vetaik[0] or
 * _model::vetak[0], and _model::cindex is set to 0.  The only public function
 * is initialize_model().  Initialization is done by randomly assigning
 * individuals (or alleles under admixture) to initial clusters, estimating the
 * parameters assuming the partition is correct
 * (initialize_parameters_mixture() or initialize_parameters_admixture()), and
 * then optionally doing further processing.  Random assignment of individuals
 * is done by a simple random partition (random_individual_partition() or
 * random_allele_partition()) or by selecting random centers and clustering by
 * distance (random_individual_center()) or identity (random_allele_center()).
 * The additional processing involves the RandEM algorithm
 * (randem_initialize_mixture() or randem_initialize_admixture()).
 *
 * TODO See Chen & Goodman "An empirical study of smoothing techniques for
 * language modeling" CAp 2006, Conférence d'Apprentissage, Harvard University,
 * Presses Universitaires de Grenoble, 1998 for more information about
 * the kind of smoothing used in this file.
 */

#include "multiclust.h"
#define MAKE_1ARRAY MAKE_1ARRAY_RETURN

void test_initialization(data *dat, model *mod);
void random_initialize_mixture(options *opt, data *dat, model *mod);
int randem_initialize_mixture(options *opt, data *dat, model *mod);
void random_initialize_admixture(options *opt, data *dat, model *mod);
int randem_initialize_admixture(options *opt, data *dat, model *mod);
void random_individual_partition(data *dat, model *mod);
void random_individual_center(data *dat, model *mod);
void random_allele_partition(data *dat, model *mod);
void random_allele_center(data *dat, model *mod);
void initialize_parameters_mixture(data *dat, model *mod);
void initialize_parameters_admixture(options *opt, data *dat, model *mod);

/**
 * Pick a random initialization.  Set _options::initialization_procedure to
 * RAND_EM to select Rand EM.  Otherwise, the random initialization is
 * unprocessed and sent back to the caller.
 *
 * @param opt options object
 * @param dat data object
 * @param mod model object
 * @return error status
 */
int initialize_model(options *opt, data *dat, model *mod)
{

	/* reset model state */
	mod->n_iter = 0;
	mod->logL = -Inf;
	mod->converged = 0;

	/* reset indices of parameter copies if using acceleration */
	if (opt->accel_scheme) {
#ifndef OLDWAY
		mod->pindex = 0;
		mod->tindex = 0;
		mod->findex = 0;
		mod->uindex = 0;	/* not necessary */
		mod->vindex = 0;	/* not necessary */
#endif
	}

	if (opt->admixture) {
		if (opt->pfile && opt->qfile) {
			read_qfile(opt, dat, mod);
			read_pfile(opt, dat, mod);
		} else if (opt->initialization_procedure == NOTHING)
			random_initialize_admixture(opt, dat, mod);
		else
			return randem_initialize_admixture(opt, dat, mod);
	} else {
		if (opt->initialization_procedure == NOTHING)
			random_initialize_mixture(opt, dat, mod);
		else
			return randem_initialize_mixture(opt, dat, mod);
	}

	return NO_ERROR;
} /* End of initialize_model(). */


/**
 * Initialize mixture model by random partition of individuals.  Setting
 * _options::initialization_method to RANDOM_PARTITION selects a completely
 * random partition.  Otherwise RANDOM_CENTERS selects _model::K random centers
 * and partitions based on proximity to the centers. 
 *
 * @param opt options object
 * @param dat data object
 * @param mod model object
 * @return void
 */
void random_initialize_mixture(options *opt, data *dat, model *mod)
{
	if (opt->initialization_method == RANDOM_PARTITION)
		random_individual_partition(dat, mod);
	else
		random_individual_center(dat, mod);
	initialize_parameters_mixture(dat, mod);
} /* End of random_initialize_mixture(). */


/**
 * Rand-EM initialization for mixture model. Run multiple initializations (by
 * choosing random centers) through one step of EM and initialize with the
 * initialization that gives the best log likelihood.
 *
 * @param opt options object
 * @param dat data object
 * @param mod model object
 * @return void
 */
int randem_initialize_mixture(options *opt, data *dat, model *mod)
{
	int i;
	int n_init = mod->K > 1 ? opt->n_rand_em_init : 1;
	int *I_K;
	double logL;
	double max_logL = -Inf;

	MAKE_1ARRAY(I_K, dat->I);

	for (i = 0; i < n_init; i++) {
		/* pick random partition */
		if (opt->initialization_method == RANDOM_PARTITION)
			random_individual_partition(dat, mod);
		else
			random_individual_center(dat, mod);

		/* use partition to initialize parameters */
		initialize_parameters_mixture(dat, mod);

		/* refine through one iteration of EM plus another E step
		 * since likelihood is computed one-step-late in E step */
		logL = em_e_step(opt, dat, mod);

		/* check for better log likelihood */
		if (logL > max_logL) {
			max_logL = logL;
			COPY_1ARRAY(I_K, dat->I_K, dat->I);
		}
	}
	COPY_1ARRAY(dat->I_K, I_K, dat->I);

	/* estimate etak and pKLM  */
	initialize_parameters_mixture(dat, mod);

	FREE_1ARRAY(I_K);

	return NO_ERROR;
} /* End of randem_initialize_mixture(). */


/**
 * Initialize the mixture model by a random partition of individuals.  This is
 * the simplest initialization method.  It randomly assigns each observation
 * (individual) to one of the K clusters.
 *
 * @param dat data object
 * @param mod model object
 * @return void
 */
void random_individual_partition(data *dat, model *mod)
{
	int i;

	for (i = 0; i < dat->I; i++)
		dat->I_K[i] = (int) (rand() % mod->K);
} /* End of random_individual_partition(). */


/**
 * Initialize the mixture model by picking random individuals.  This method
 * randomly picks K "center" individuals and then forms clusters by assigning
 * all other individuals to the closest center individual.  It populates
 * mod->I_K with the assignments.
 *
 * @param dat data object
 * @param mod model object
 * @return void
 */
void random_individual_center(data *dat, model *mod)
{
	int i, j, k, l, m, flag, m_start;
	int center[mod->K];
	double count_diff, min_count_diff;

	if (mod->K == 1) {
		for (i = 0; i < dat->I; i++)
			dat->I_K[i] = 0;
		return;
	}

	/* pick unique center individuals */
	for (k = 0; k < mod->K; k++) {
		center[k] = (int) (rand() % dat->I);
		do {
			flag = 0;
			for (j = 0; j < k; j++) {
				if (center[k] == center[j]) {
					center[k] = (int) (rand() % dat->I);
					flag = 1;
					break;
				}
			}
		} while (flag == 1);
	}

	/* assign individuals to closest center */
	for (i = 0; i < dat->I; i++) {
		dat->I_K[i] = 0;
		if(i == center[0])
			continue;

		/* find closest center */
		min_count_diff = Inf;
		for (k = 0; k < mod->K; k++) {
			/* this individual is a center */
			if (i == center[k]) {
				dat->I_K[i] = k;
				break;
			}

			/* compute distance to kth center */
			count_diff = 0;
			for (l = 0; l < dat->L; l++) {
				m_start = dat->L_alleles 
					&& dat->L_alleles[l][0] == MISSING
					? 1 : 0;
				for (m = m_start; m < dat->uniquealleles[l];
					m++)
					count_diff += abs(dat->ILM[i][l][m]
						- dat->ILM[center[k]][l][m]);
				/* under uniform prior model, the missing data
				 * has chance of matching */
				if (m_start)
					count_diff -= abs(dat->ILM[i][l][0]
						- dat->ILM[center[0]][l][0])
						/(dat->uniquealleles[l]-1.0);
			}

			/* record if this is closer center */
			if (count_diff < min_count_diff) {
				dat->I_K[i] = k;
				min_count_diff = count_diff;
			}
		}
	}
} /* End of random_individual_center(). */

/**
 * Initialize model parameters based on initial partition of individuals.
 *
 * @param dat data object
 * @param mod model object
 * @return void
 */
void initialize_parameters_mixture(data *dat, model *mod)
{
	int k, l, m, i, m_start;
	double temp;

	/* estimate estak */
	for (k = 0; k < mod->K; k++)
#ifndef OLDWAY
		mod->vetak[mod->tindex][k] = 1;	/* avoid 0 */
#else
		mod->etak[k] = 1;	/* avoid 0 */
#endif

	for (i = 0; i < dat->I; i++)
#ifndef OLDWAY
		mod->vetak[mod->tindex][dat->I_K[i]]++;
#else
		mod->etak[dat->I_K[i]]++;
#endif

	for (k = 0; k < mod->K; k++)
#ifndef OLDWAY
		mod->vetak[mod->tindex][k] /= dat->I + mod->K;
#else
		mod->etak[k] /= dat->I + mod->K;
#endif

	/* estimate pKLM */
	for (k = 0; k < mod->K; k++)
		for (l = 0; l < dat->L; l++) {
			/* ignore possible missing data */
			m_start = (dat->L_alleles
				&& dat->L_alleles[l][0] == MISSING);
			for (m = m_start; m < dat->uniquealleles[l]; m++) {
#ifndef OLDWAY
				mod->vpklm[mod->tindex][k][l][m] = 1.0;	/* avoid 0 */
#else
				mod->pKLM[k][l][m] = 1.0;	/* avoid 0 */
#endif
				for (i = 0; i < dat->I; i++) {
					if (dat->ILM[i][l][m] == 0)
						continue;
#ifndef OLDWAY
					mod->vpklm[mod->tindex][dat->I_K[i]][l][m]
#else
					mod->pKLM[dat->I_K[i]][l][m]
#endif
						+= dat->ILM[i][l][m];
				}
			}
		}

	/* normalize pKLM */
	for (k = 0; k < mod->K; k++)
		for (l = 0; l < dat->L; l++) {
			m_start = (dat->L_alleles
				&& dat->L_alleles[l][0] == MISSING);
			temp = 0.0;
			for (m = m_start; m < dat->uniquealleles[l]; m++)
#ifndef OLDWAY
				temp +=  mod->vpklm[mod->tindex][k][l][m];
#else
				temp +=  mod->pKLM[k][l][m];
#endif
			for(m = m_start; m < dat->uniquealleles[l]; m++)
#ifndef OLDWAY
				mod->vpklm[mod->tindex][k][l][m] /= temp;
#else
				mod->pKLM[k][l][m] /= temp;
#endif
		}
} /* End of initialize_parameters_mixture(). */

/**
 * Initialize parameters by randomly assigning alleles to populations.
 *
 * @param opt options object
 * @param dat data object
 * @param mod model object
 * @return void
 */
void random_initialize_admixture(options *opt, data *dat, model *mod)
{
	if (opt->initialization_method == TESTING)
		test_initialization(dat, mod);
	else
		random_allele_partition(dat, mod);

	m_step_admixture(opt, dat, mod);
} /* End of random_initialize_admixture(). */

void test_initialization(data *dat, model *mod)
{
	int debug = 0;
	int i, l, m, k, a, s = 0;

	for (i = 0; i < dat->I; i++) {
		for (l = 0; l < dat->L; l++) {
			int m_start = dat->L_alleles
				&& dat->L_alleles[l][0] == MISSING;
			for (m = m_start; m < dat->uniquealleles[l]; m++)
				for (k = 0; k < mod->K; k++)
					mod->diklm[i][k][l][m] = 0;
			for (a = 0; a < dat->ploidy; a++) {
				k = s++ % mod->K;
				/* alleles are not array indices */
				if (dat->L_alleles) {
					for (m = m_start;
						m < dat->uniquealleles[l]; m++)
						if (dat->IL[2*i + a][l]
							== dat->L_alleles[l][m]) {
							mod->diklm[i][k][l][m] = 1;
							if (debug)
								fprintf(stdout,
									"%d %d "
									"%d: %d*\n",
									i, l,
									dat->L_alleles[l][m],
									k);
						}
				/* alleles are array indices */
				} else {
					mod->diklm[i][k][l][dat->IL[2*i + a][l]]
						= 1;
					if (debug)
						fprintf(stdout, "%d %d %d: %d\n",
							i, l, dat->IL[2*i + a][l], k);
				}
			}
		}
	}
} /* test_initialization */

/**
 * Randomly initialize algorithm using RandEM.  For _options::n_init times, randomly initialize the parameters using
 * random_allele_center(), run one round of EM, and record 
 * the initialization with highest log likelihood.

 *
 * @param opt options object
 * @param dat data object
 * @param mod model object
 * @return error status
 */
int randem_initialize_admixture(options *opt, data *dat, model *mod)
{
	int i;
	int n_init = mod->K > 1 ? opt->n_rand_em_init : 1;
	int n_haplotypes = dat->I * dat->ploidy;
	int **IL_K_max;
	double logL;
	double max_logL = -Inf;

	MAKE_2ARRAY(IL_K_max, n_haplotypes, dat->L);

	for (i = 0; i < n_init; i++) {
		/* randomly initialize */
		random_allele_center(dat, mod);
		initialize_parameters_admixture(opt, dat, mod);

		/* one iteration of EM */
		logL = em_e_step(opt, dat, mod);

		if (logL > max_logL) {
			max_logL = logL;
			COPY_2ARRAY(IL_K_max, dat->IL_K, dat->L);
		}
	}
	COPY_2ARRAY(dat->IL_K, IL_K_max, dat->L);

	/* Update etak and KLM by the I_K which has the highest logL. */
	initialize_parameters_admixture(opt, dat, mod);

	FREE_2ARRAY(IL_K_max);

	return NO_ERROR;
} /* End of randem_initialize_admixture(). */


/**
 * Randomly initialize the admixture model with random partition.  This is the
 * simplest initialization method.  It randomly assigns each observed allele to
 * one of the K clusters.
 *
 * @param dat data object
 * @param mod model object
 * @return void
 */
void random_allele_partition(data *dat, model *mod)
{
	int i, l, k, a, m, m_start;

	for (i = 0; i < dat->I; i++)
		for (l = 0; l < dat->L; l++) {
			m_start = dat->L_alleles && dat->L_alleles[l][0] == MISSING;
			for (k = 0; k < mod->K; k++)
				for (m = m_start; m < dat->uniquealleles[l]; m++)
					mod->diklm[i][k][l][m] = 0;
			for (a = 0; a < dat->ploidy; a++) {
				k = (int) rand() % mod->K;
				/* alleles are not array indices */
				if (dat->L_alleles) {
					for (m = m_start; m < dat->uniquealleles[l]; m++)
						if (dat->IL[2*i + a][l] == dat->L_alleles[l][m]) {
							mod->diklm[i][k][l][m] = 1;
//fprintf(stderr, "%d %d %d: %d\n", i, l, dat->IL[2*i + a][l], k);
						}
				/* alleles are array indices */
				} else {
					mod->diklm[i][k][l][dat->IL[2*i + a][l]] = 1;
				}
//				dat->IL_K[2*i+a][l] = (int) rand() % mod->K;
			}
		}
} /* End of random_allele_partition(). */

/**
 * Randomly intialize the admixture model with random centers.  This method
 * picks a "center" allele from the uniquealleles and then randomly associates
 * that allele to subpopulation k.  Then, observed alleles are assigned to
 * clusters based on identity with the centers.
 *
 * TODO There is nothing particular optimal about this approach.
 *
 * @param dat data object
 * @param mod model object
 * @return void
 */
void random_allele_center(data *dat, model *mod)
{
	int i, j, l, k, z_start, z_end, i_l, flag;
	int m_start;
	int n_haplotypes = dat->I * dat->ploidy;
	int center[mod->K];

	/* one subpopulation: all alleles come from subpopulation 0 */
	if (mod->K == 1) {
		for (i=0; i < n_haplotypes; i++)
			for (l = 0; l < dat->L; l++)
				dat->IL_K[i][l] = 0;
		return;
	}

	/* for each locus */
	for (l = 0; l < dat->L; l++) {
		/* check for missing data at this locus */
		m_start = dat->L_alleles && dat->L_alleles[l][0] == MISSING
			? 1 : 0;

		/* pick center allele for each cluster from uniquealleles */
		/* fewer alleles than centers; use all alleles */
		if (dat->uniquealleles[l] - m_start < mod->K) {
			for (k = m_start; k < dat->uniquealleles[l]; k++)
				center[k] = k;
			for (k = dat->uniquealleles[l]; k < mod->K; k++)
				center[k] = -1;
		/* more alleles than centers; assign alleles randomly */
		} else {
			for (k = 0; k < mod->K; k++) {

				/* choose center allele at random */
				center[k] = m_start + (int) rand()
					% (dat->uniquealleles[l] - m_start);

				/* make sure it is unique */
				do {
					flag = 0;
					for (j = 0; j < k; j++) {
						if (center[k] == center[j]) {
							center[k] = m_start
								+ (int) rand()
								% (dat->uniquealleles[l] - m_start);
							flag = 1;
							break;
						}
					}
				} while (flag == 1);
			}
		}

		/* assign alleles by identity with centers, or randomly */
		z_end = 0;
		for (i = 0; i < dat->I; i++) {
			z_start = z_end;
			z_end += dat->ploidy;
			/* for each allele */
			for (i_l = z_start; i_l < z_end; i_l++) {
				flag = 0;
				/* assign to cluster when it matches */
				for (k = 0; k < mod->K; k++) {
					if (center[k] == -1)
						break;
					/* allele matches a center */
					if (dat->L_alleles && dat->IL[i_l][l]
						== dat->L_alleles[l][center[k]]) {

						dat->IL_K[i_l][l] = k;
						flag = 1;
						break;
					} else if (!dat->L_alleles &&
						dat->IL[i_l][l] == center[k]) {
						dat->IL_K[i_l][l] = k;
						flag = 1;
						break;
					}
				}
				/* no matching center allele: assign randomly */
				if (flag == 0)
					dat->IL_K[i_l][l] = (int) rand() % mod->K;
			}
		}
	}
} /* End of random_allele_center(). */

/**
 * Initialize parameters given initial assignment of alleles to populations.
 *
 * @param opt options object
 * @param dat data object
 * @param mod model object
 * @return void
 */
void initialize_parameters_admixture(options *opt, data *dat, model *mod)
{
	int k, l, m, i, i_l, m_start, z_start = 0, z_end = 0;
	int n_haplotypes = dat->I * dat->ploidy;
	double temp = 0.0;


	/* estimate etaks or etaiks */
	if (opt->eta_constrained) {
		temp = dat->ploidy*dat->L*dat->I + mod->K;	/* was BUG: normalized wrong */
		for (k=0; k<mod->K; k++)
#ifndef OLDWAY
			mod->vetak[mod->tindex][k] = 1;
#else
			mod->etak[k] = 1;
#endif
		for (i=0; i<n_haplotypes; i++)
			for (l=0; l<dat->L; l++)
#ifndef OLDWAY
				mod->vetak[mod->tindex][dat->IL_K[i][l]]++;
#else
				mod->etak[dat->IL_K[i][l]]++;
#endif
		for (k=0; k<mod->K; k++)
#ifndef OLDWAY
			mod->vetak[mod->tindex][k] /= temp;
#else
			mod->etak[k] /= temp;
#endif
	} else {
		temp = dat->ploidy*dat->L + mod->K;	/* was BUG: normalized wrong */
		for (i = 0; i < dat->I; i++) {	/* individual */
			z_start = z_end;
			z_end += dat->ploidy;

			for (k = 0; k < mod->K; k++)
#ifndef OLDWAY
				mod->vetaik[mod->tindex][i][k] = 1;	/* avoid 0 */
#else
				mod->etaik[i][k] = 1;	/* avoid 0 */
#endif

			for (i_l = z_start; i_l < z_end; i_l++)
				for (l = 0; l < dat->L; l++)
#ifndef OLDWAY
					mod->vetaik[mod->tindex][i][dat->IL_K[i_l][l]]++;
#else
					mod->etaik[i][dat->IL_K[i_l][l]]++;
#endif

			/* normalize */
			for (k = 0; k < mod->K; k++)
#ifndef OLDWAY
				mod->vetaik[mod->tindex][i][k] /= temp;
#else
				mod->etaik[i][k] /= temp;
#endif
		}
	}

	/* estimate pKLM */
	for (l = 0; l < dat->L; l++) {				/* locus */
		/* 0th allele may be missing */
		m_start = dat->L_alleles && dat->L_alleles[l][0] == MISSING;
		for (m = m_start; m < dat->uniquealleles[l]; m++)/* allele */
			for (k = 0; k < mod->K; k++)		/* population */
#ifndef OLDWAY
				mod->vpklm[mod->tindex][k][l][m] = 1.0;	/* avoid 0 */
#else
				mod->pKLM[k][l][m] = 1.0;	/* avoid 0 */
#endif
	}

	for (i = 0; i < n_haplotypes; i++)	/* haplotype */
		for (l = 0; l < dat->L; l++) {	/* locus */
			/* allele names are not indices */
			if (dat->L_alleles) {
				m_start = dat->L_alleles && dat->L_alleles[l][0] == MISSING;
				/* [TODO] Inefficient.  Store indices in IL and map to
				 * allele name as it appears in data file only when
				 * interfacing with user. */
				for (m = m_start; m < dat->uniquealleles[l]; m++)
					if (dat->IL[i][l] == dat->L_alleles[l][m])
#ifndef OLDWAY
						mod->vpklm[mod->tindex][dat->IL_K[i][l]][l][m]++;
#else
						mod->pKLM[dat->IL_K[i][l]][l][m]++;
#endif
			} else {
#ifndef OLDWAY
				mod->vpklm[mod->tindex][dat->IL_K[i][l]][l][dat->IL[i][l]]++;
#else
				mod->pKLM[dat->IL_K[i][l]][l][dat->IL[i][l]]++;
#endif
			}
		}

	/* normalize */
	for (k = 0; k < mod->K; k++)		/* population */
		for (l = 0; l < dat->L; l++) {	/* locus */
			m_start = dat->L_alleles && dat->L_alleles[l][0] == MISSING;
			temp = 0;
			for (m = m_start; m < dat->uniquealleles[l]; m++)
#ifndef OLDWAY
				temp += mod->vpklm[mod->tindex][k][l][m];
#else
				temp += mod->pKLM[k][l][m];
#endif
			for (m = m_start; m < dat->uniquealleles[l]; m++)
#ifndef OLDWAY
				mod->vpklm[mod->tindex][k][l][m] /= temp;
#else
				mod->pKLM[k][l][m] /= temp;
#endif
		}
} /* End of initialize_parameters_admixture(). */
