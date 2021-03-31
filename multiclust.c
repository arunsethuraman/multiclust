/**
 * @file multiclust.c
 * @author Karin Dorman, kdorman@iastate.edu
 * @date Wed Dec  5 09:20:48 CST 2012
 *
 * User interface to multiclust.
 *
 * TODOS
 * [TODO] write the results of parametric bootstrap to a permanent output file
 * [TODO] run matrix of hypothesis tests and control FDR (as per Maitra & 
 *        Melnykov) to auto-estimate K
 * [TODO] add hypothesis test for H0: mixture vs. HA: admixture
 * [TODO] add code to run variable number of initializations, perhaps until the
 *        current best log likelihood has been observed X times
 */


#include "cline.h"	/* command line */
#include "multiclust.h"
#define MAKE_1ARRAY MAKE_1ARRAY_RETURN

/* create structures for options, data, and model */
int make_options(options **opt);
int make_data(data **dat);
int make_model(model **);

/* set/parse options and model; verify integrity */
int parse_options(options *opt, data *dat, int argc, const char **argv);
int allocate_model_for_k(options *opt, model *mod, data *dat);
int synchronize(options *opt, data *dat, model *mod);

/* do the data fitting or bootstrap */
int timed_model_estimation(options *, data *, model *);
int estimate_model(options *opt, data *dat, model *mod, int bootstrap);
int maximize_likelihood(options *opt, data *dat, model *mod, int bootstrap);
int run_bootstrap(options *opt, data *dat, model *mod);
void print_model_state(options *opt, data *dat, model *mod, int diff, int newline);
double adj_rand(int n, int k1, int k2, int *cl1, int *cl2, int type);

/* cleanup all allocated memory */
void free_options(options *opt);
void free_data(data *dat);
void free_model(model *mod, options *opt);
void free_model_mles(model *mod);
void free_model_data(model *mod, options *opt);

const char *accel_method_abbreviations[NUM_ACCELERATION_METHODS] = {
	"EM",
	"S1",
	"S2",
	"S3",
	"Q",
};
const char *accel_method_names[NUM_ACCELERATION_METHODS] = {
	"No acceleration",
	"SQUAREM version 1",
	"SQUAREM version 2",
	"SQUAREM version 3",
	"Quasi Newton",
};




int main(int argc, const char **argv)
{
	options *opt = NULL;	/* run options */
	data *dat = NULL;	/* genetic data */
	model *mod = NULL;	/* model parameters */
	int err = NO_ERROR;	/* error code */

/*
	for (int i=0; i<argc; i++)
		fprintf(stderr, " %s", argv[i]);
	fprintf(stderr, "\n\n");
	exit(0);
*/

#ifdef OLDWAY
fprintf(stdout, "Using OLDWAY code\n");
#endif

	/* make various structures to store run information */
	if ((err = make_options(&opt)))
		goto FREE_AND_EXIT;

	/* set up data structure */
	if ((err = make_data(&dat)))
		goto FREE_AND_EXIT;

	/* set up model structure */
	if ((err = make_model(&mod)))
		goto FREE_AND_EXIT;

	/* parse command-line options */
	if ((err = parse_options(opt, dat, argc, argv)))
		goto FREE_AND_EXIT;

	/* read data */
	if ((err = read_file(opt, dat)))
		goto FREE_AND_EXIT;

	/* finalize settings that refer to data; check settings */
	if ((err = synchronize(opt, dat, mod)))
		goto FREE_AND_EXIT;

	/* estimate the model(s) using the observed data */
	if (opt->n_repeat > 1 && (err = timed_model_estimation(opt, dat, mod)))
		goto FREE_AND_EXIT;
	if (opt->n_repeat == 1 && (err = estimate_model(opt, dat, mod, 0)))
		goto FREE_AND_EXIT;

	/* set up parallelization mode */
	if (opt->parallel)
		err = estimate_model(opt, dat, mod, 0);
		if (err)
			goto FREE_AND_EXIT;
		printf("%f\n", mod->max_logL);
		goto FREE_AND_EXIT;


	/* optionally run a bootstrap */
	if (opt->n_bootstrap) {
		if ((err = run_bootstrap(opt, dat, mod)))
			goto FREE_AND_EXIT;

		/* TODO write this result to some file! */
		fprintf(stdout, "p-value to reject H0: K=%d is %f\n", 
			mod->null_K, mod->pvalue);
	}

FREE_AND_EXIT:

	free_model(mod,opt);
	free_options(opt);
	free_data(dat);

	return err;

} /* main */

/**
 * Time model estimation by repeating multiple times.  Repeat model-fitting
 * process chosen by user _options:n_repeat, compute the total time, and report
 * the average time.  The purpose of this code is to produce better timing
 * information by replication.  It is not useful for model fitting, as one
 * would normally not want to repeat the same model fit!  Currently, this
 * function cannot be used to time the bootstrap process.
 *
 * @param opt options object
 * @param dat data object
 * @param mod model object
 * @return error status
 */
int timed_model_estimation(options *opt, data *dat, model *mod)
{
	int err = NO_ERROR;
	clock_t start;
	int enough_time = opt->repeat_seconds ? 0 : 1;
	double esec = 0;
	int max_iter = 0;
	int max_init = 0;
	double sum_init = 0;
	double sum_init2 = 0;
	double sum_iter = 0;
	double sum_iter2 = 0;
	double sum_aic_K = 0;
	double sum_aic_K2 = 0;
	double sum_bic_K = 0;
	double sum_bic_K2 = 0;
	double sum_ll = 0;
	double sum_ll2 = 0;
	double sum_ar = 0;
	double sum_ar2 = 0;
	double max_ll = -Inf;
	double min_aic = 0, min_bic = 0;
	double max_ar = -1;
	double first_ll = -Inf;
	double max_ll_rand = 0;
	int first_hit_index = 0;
	int converged_repeats = 0;
	int n_repeats = 0;
	int target_reached = 0;

	/* time multiple runs */
	start = clock();

	for (; n_repeats < opt->n_repeat || !enough_time;) {
		if ((err = estimate_model(opt, dat, mod, 0)))
			return err;
		if (mod->n_init > max_init)
			max_init = mod->n_init;
		if (mod->n_max_iter > max_iter)
			max_iter = mod->n_max_iter;
		/* [TODO] assumes no repetitions inside estimate_mode() */
		if (mod->max_logL > max_ll) {
			max_ll = mod->max_logL;
			min_aic = mod->aic;
			min_bic = mod->bic;
			max_ll_rand = mod->arand;
			if (!converged(opt, mod, first_ll)) {
				first_ll = mod->max_logL;
				first_hit_index = n_repeats;
			}
		}
		if (mod->arand > max_ar)
			max_ar = mod->arand;
		sum_init += mod->n_init;
		sum_init2 += mod->n_init * mod->n_init;
		sum_iter += mod->n_total_iter;
		sum_iter2 += mod->n_total_iter * mod->n_total_iter;
		sum_aic_K += mod->aic_K;
		sum_aic_K2 += mod->aic_K * mod->aic_K;
		sum_bic_K += mod->bic_K;
		sum_bic_K2 += mod->bic_K * mod->bic_K;
		sum_ll += mod->max_logL;
		sum_ll2 += mod->max_logL * mod->max_logL;
		if (opt->afile) {
			sum_ar += mod->arand;
			sum_ar2 += mod->arand * mod->arand;
		}
		n_repeats++;
		if (mod->ever_converged)
			converged_repeats++;
		if (mod->n_targetll_times)
			target_reached++;

		esec = ((double) clock() - start) / CLOCKS_PER_SEC;
		if (opt->verbosity > SILENT) {
			print_model_state(opt, dat, mod, 
				((double) clock() - start)/CLOCKS_PER_SEC, 0);
			fprintf(stdout, " %f %f %d %d", esec, esec/n_repeats,
				target_reached, converged_repeats);
			fprintf(stdout, " %f", max_ll);
			if (opt->target_ll)
				fprintf(stdout, " %f", opt->desired_ll);
			else
				fprintf(stdout, " NA");
			fprintf(stdout, " %d %d %d %d", opt->target_revisit,
				n_repeats, opt->n_repeat, opt->repeat_seconds);
			fprintf(stdout, "\n");
		}

		/* check time spent against required min and max*/
		if (!enough_time || opt->max_repeat_seconds) {
			if (!enough_time && esec > opt->repeat_seconds)
				enough_time = 1;
			if (opt->max_repeat_seconds && esec > opt->max_repeat_seconds)
				break;

		}
	}

	if (opt->verbosity >= SILENT) {
		fprintf(stdout, "Data, Method, Model: %s, %s, %s\n", opt->filename,
			opt->accel_abbreviation,
			opt->admixture && opt->eta_constrained
			?"admix constrained":opt->admixture?"admix":"mix");
		fprintf(stdout, "Run: %e %e %e %e n=%d i=%d u=(%f,%d) w=(%d,%d)\n", 
			opt->abs_error, opt->rel_error, opt->eta_lower_bound,
			opt->p_lower_bound, opt->n_init, opt->n_init_iter, opt->desired_ll,
			opt->target_revisit, opt->n_repeat, opt->repeat_seconds);
		fprintf(stdout, "Number of repetitions: %d of %d requested, %d converged, %d reach target\n", 
			n_repeats, opt->n_repeat, converged_repeats, target_reached);
		fprintf(stdout, "Average time: %fs (total: %fs; target: %d)\n", 
			esec/n_repeats, esec, opt->repeat_seconds);
		fprintf(stdout, "Average log likelihood: %f (+/- %f)\n",
			sum_ll/n_repeats, sqrt((sum_ll2 - sum_ll*sum_ll/n_repeats)
			/(n_repeats-1)));
		fprintf(stdout, "Maximum log likelihood: %f first hit at run %d (AIC %f; BIC %f; RAND: %f)\n", 
			max_ll, first_hit_index, min_aic, min_bic, max_ll_rand);
		fprintf(stdout, "Adjusted RAND: avg = %f +/- %f; max = %f\n",
			sum_ar/n_repeats, sqrt((sum_ar2 - sum_ar*sum_ar/n_repeats)
			/(n_repeats-1)), max_ar);
		if (opt->max_K != opt->min_K) {
			fprintf(stdout, "Average K (AIC): %f (+/- %f)\n",
				sum_aic_K/n_repeats, sqrt((sum_aic_K2 - 
				sum_aic_K*sum_aic_K/n_repeats)/(n_repeats - 1)));
			fprintf(stdout, "Average K (BIC): %f (+/- %f)\n",
				sum_bic_K/n_repeats, sqrt((sum_bic_K2 - sum_bic_K
				* sum_bic_K / n_repeats)/(n_repeats - 1)));
		} else {
			fprintf(stdout, "Total initializations, iterations: %d, %d\n",
				(int)sum_init, (int)sum_iter);
			fprintf(stdout, "Average initializations: %f (+/- %f)"
				" [%e, %e]\n", sum_init/n_repeats, sqrt((sum_init2
				- sum_init*sum_init/n_repeats)/(n_repeats - 1)),
				sum_init2, sum_init);
			fprintf(stdout, "Average iterations: %f (+/- %f) [%e %e]\n",
				sum_iter/sum_init, sqrt((sum_iter2 - sum_iter
				* sum_iter/sum_init)/(sum_init - 1)),
				sum_iter2, sum_iter);
			fprintf(stdout, "Maximum initializations: %d\n",
				max_init);
			fprintf(stdout, "Maximum iterations: %d\n", max_iter);
		}
	}

	return err;
} /* timed_model_estimation */

/**
 * Estimate model(s) using maximum likelihood.  For each model, the function
 * allocates the model (parameters), maximizes the likelihood (probably using
 * multiple initializations), extracts relevant statistics (mles, log
 * likelihood), and frees the model.  Currently, the code either tests H0
 * against HA (and aborts if the HA log likelihood does not exceed the H0 log
 * likelihood) or fits a series of models from K=_options:min_K to K=_options::max_K.
 * Importantly, when hypothesis testing, it fits H0 //before// HA.  One could
 * modify this code to select other models for fitting.
 *
 * @param opt options object
 * @param dat data object
 * @param mod model object
 * @param bootstrap indicate if bootstrap run
 * @return error status
 */
int estimate_model(options *opt, data *dat, model *mod, int bootstrap)
{
	int total_iter = 0;
	int err = NO_ERROR;
	clock_t start;
	double min_bic, min_aic;

	start = clock();

	/* initialize: if bootstrap, comparing H0: K = mod->null_K vs.
	 * Ha: K = mod->alt_K; otherwise running for 1, 2, ..., opt->max_K
	 * so log likelihood will increase as progress through models */
	mod->max_logL = -Inf;
	mod->max_logL_H0 = -Inf;
	mod->K = opt->n_bootstrap ? mod->null_K : opt->min_K;
	dat->max_M = dat->M;
	min_aic = min_bic = Inf;

	//if(opt->block_relax == 1)
	//	mod->K = 2;

	do {
		/* Setting the largest size of the U and V arrays to be the larger of K and max_M */
		if (dat->max_M < mod->K)
			dat->max_M = mod->K;
		
		/* knowing K, can allocate space for parameters */
		if ((err = allocate_model_for_k(opt, mod, dat)))
			return err;

		/* maximize the likelihood; involves multiple initializations */
		/* also stores mles in mod->mle_* parameters */
		if ((err = maximize_likelihood(opt, dat, mod, bootstrap)))
			return err;

		/* final output for this model */
		if (opt->n_repeat == 1 && opt->verbosity) /* || opt->verbosity > SILENT) */
			print_model_state(opt, dat, mod, 
				((double) clock() - start)/CLOCKS_PER_SEC, 1);

		total_iter += mod->n_total_iter;


		/* possibly store maximum likelihood under H0 */
		if (opt->n_bootstrap && mod->K == mod->null_K)
			mod->max_logL_H0 = mod->max_logL;

		if (min_aic > mod->aic) {
			min_aic = mod->aic;
			mod->aic_K = mod->K;
		}
		if (min_bic > mod->bic) {
			min_bic = mod->bic;
			mod->bic_K = mod->K;
		}

		/* free parameters */
		free_model_data(mod, opt);

		/* choose next K */
		if (opt->n_bootstrap && mod->K == mod->null_K)
			mod->K = mod->alt_K;
		else if (!opt->n_bootstrap && mod->K < opt->max_K)
			mod->K++;
		else
			break;

	} while (1);

	/* if running bootstrap, record test statistic */
	if (opt->n_bootstrap) {
		double diff = mod->max_logL - mod->max_logL_H0;

		if (diff <= 0)	/* == OK, but prob. convergence error */
			return message(stderr, __FILE__, __func__, __LINE__,
				ERROR_MSG, INTERNAL_ERROR, "Null hypothesis "
				"likelihood exceeds alternative hypothesis "
				"likelihood.  Try increasing number of "
				"initializations (command-line option -n)\n");

		if (!bootstrap)
			mod->ts_obs = diff;
		else 
			mod->ts_bs = diff;
	}

	return err;
} /* estimate_model */


/**
 * Maximize the log likelihood.  Given a particular model (admixture/mixture
 * and _model::K), initialize the model _options::n_init initializations times,
 * maximize the likelihood, and record the result if a better solution is
 * found.  The best log likelihood is stored in _model::max_logL.  The mle
 * parameter estimates are stored in _model::mle_pKLM and _model::mle_etak (or
 * _model::mle_etaik), but only if the observed data is being fit to H0.  In
 * addition, if fitting the observed data, the best fitting solution is written
 * to file.
 *
 * @param opt options object
 * @param dat data object
 * @param mod model object
 * @param bootstrap is this a bootstrap run
 * @return error status
 */
int maximize_likelihood(options *opt, data *dat, model *mod, int bootstrap)
{
	int i;
	int err = NO_ERROR;


	/* resetting statistics recorded across multiple initializations */
        mod->first_max_logL = -Inf;
	mod->n_init = 0;
	mod->n_total_iter = 0;
	mod->n_maxll_times = 0;
	mod->n_maxll_init = -1;
	mod->n_targetll_times = 0;
	mod->n_max_iter = 0;
	mod->time_stop = 0;
	mod->ever_converged = 0;

	mod->start = clock();
	for (i=0; opt->target_revisit || opt->target_ll		/* targeting */
		|| opt->n_seconds || i<opt->n_init; i++) {	/* timing */
		mod->current_i = 0;
		mod->current_l = 0;
		mod->current_k = 0;
		mod->logL = 0.0;
		mod->converged = 0;
		mod->stopped = 0;
		mod->iter_stop = 0;

		/* initialize parameters */
		if ((err = initialize_model(opt, dat, mod)))
			return err;

		/* maximize likelihood */
		em(opt, dat, mod);
		/* !_model::converged b/c _model::iter_stop || _model::time_stop */

		if (mod->converged)
			mod->ever_converged = 1;

		/* add iter/init statistics if converged or timed out */
		if (mod->converged || (!mod->n_init && mod->time_stop)) {
			mod->n_total_iter += mod->n_iter;
			if (mod->n_max_iter < mod->n_iter)
				mod->n_max_iter = mod->n_iter;
			mod->n_init++;
		}

		/* solution already seen (up to convergence precision) */
		if (mod->converged && converged(opt, mod, mod->first_max_logL)) {
			mod->n_maxll_times++;
		/* first occurrence of better solution */
		} else if (mod->converged && mod->logL > mod->first_max_logL) {
			mod->n_maxll_times = 1;
			mod->first_max_logL = mod->logL;
			mod->n_maxll_init = mod->n_init;
		}

		/* record any better solution than previously seen */
		if (mod->logL > mod->max_logL) {
			mod->max_logL = mod->logL;
			mod->aic = aic(mod);
			mod->bic = bic(dat, mod);

			/* save mles if not bootstrap run and this is H0 */
			if (!bootstrap && opt->n_bootstrap && mod->K == mod->null_K) {
#ifndef OLDWAY
				COPY_3JAGGED_ARRAY(mod->mle_pKLM, mod->vpklm[mod->pindex], dat->uniquealleles);
#else
				COPY_3JAGGED_ARRAY(mod->mle_pKLM, mod->pKLM, dat->uniquealleles);
#endif
				if (!opt->admixture || opt->eta_constrained)
#ifndef OLDWAY
					COPY_1ARRAY(mod->mle_etak, mod->vetak[mod->pindex], mod->K);
#else
					COPY_1ARRAY(mod->mle_etak, mod->etak, mod->K);
#endif
				else
#ifndef OLDWAY
					COPY_2ARRAY(mod->mle_etaik, mod->vetaik[mod->pindex], mod->K);
#else
					COPY_2ARRAY(mod->mle_etaik, mod->etaik, mod->K);
#endif
			}

			/* write results to file if not bootstrap run */
			if (!bootstrap && opt->write_files) {
				/* TODO [KSD]: overwriting potentially many 
				 * times for big data is bad */
				if (opt->admixture) {
					partition_admixture(dat, mod);
					write_file_detail(opt, dat, mod);
					popq_admix(opt, dat, mod);
					indivq_admix(opt, dat, mod);
				} else {
					partition_mixture(dat, mod);
					write_file_detail(opt, dat, mod);
					popq_mix(opt, dat, mod);
					indivq_mix(opt, dat, mod);
				}
			
			}

			if (opt->afile) {
				if (!opt->write_files) {
					if (opt->admixture) 
						partition_admixture(dat, mod);
					else
						partition_mixture(dat, mod);
				}
				mod->arand = adj_rand(dat->I, opt->pK, mod->K,
					opt->partition_from_file, dat->I_K,
					ADJUSTED_RAND_INDEX);
			}
		}

//print_param(opt, dat, mod, mod->tindex);

		/* output information if sufficiently verbose and not -w */
		if (!bootstrap && opt->verbosity > QUIET && opt->write_files)
			fprintf(stdout, "K = %d, initialization = %d: %f "
				"(%s) in %3d iterations, %02d:%02d:%02d (%f; %d)\n", 
				mod->K, i, mod->logL,
				mod->converged?"converged":"not converged",
				mod->n_iter, 
				(int)(mod->seconds_run/3600),
				(int)((((int)mod->seconds_run)%3600)/60),
				(((int)mod->seconds_run)%60), mod->max_logL,
				mod->n_maxll_times);
		
		/* global mle if K=1; no need for multiple initializations */
		if (mod->K == 1)
			break;

		/* stop if used too much time */
		if (mod->time_stop)
			break;

		/* stop if sufficient repeat of same best solution */
		if (opt->target_revisit
			&& mod->n_maxll_times >= opt->target_revisit)
			break;

		/* stop if reach (sufficient times) target log likelihood */
		if (opt->target_ll && (mod->logL > opt->desired_ll
			|| converged(opt, mod, opt->desired_ll))) {
			if (!mod->n_targetll_times)
				mod->n_targetll_init = mod->n_init;
			mod->n_targetll_times++;
			if (!opt->target_revisit)
				break;
			else if (opt->target_revisit <= mod->n_targetll_times)
				break;
		}

	}

	return err;
} /* End of maximize_likelihood(). */


/**
 * Parametric bootstrap.  For __options::n_bootstrap times, simulate bootstrap
 * dataset using parameter estimates stored in _model::mle_etak or
 * _model::mle_etaik and _model::mle_pKLM.  Fit the H0 and HA models using
 * exactly the same procedure used to fit the observed data to the same models.
 * Record the number of times the test statistic, stored in _model::ts_bs,
 * exceeds the observed test statistic, stored in _model::ts_obs.  Store the
 * resulting estimated p-value in _model::pvalue.  While running, this function
 * doubles the amount of memory required to store the data, but it clears the
 * memory when done.
 *
 * @param opt options object
 * @param dat data object
 * @param mod model object
 * @return error status
 */
int run_bootstrap(options *opt, data *dat, model *mod)
{
	int i;
	int ntime = 0;
	int err = NO_ERROR;

	for (i=0; i < opt->n_bootstrap; i++) {
		fprintf(stdout, "Bootstrap dataset %d (of %d):", i+1,
			opt->n_bootstrap);

		/* generate bootstrap dataset under H0 */
		if ((err = parametric_bootstrap(opt, dat, mod)))
			return err;
/* temporary : if you want to generate some simulation data to play with
write_data(opt, dat, 1);
*/

		/* fit models H0 and HA */
		if ((err = estimate_model(opt, dat, mod, 1)))
			return err;

		/* mod->max_logL is maximum log likelihood under HA */
		if (mod->ts_bs >= mod->ts_obs)
			ntime++;
		fprintf(stdout, " test statistics bs=%f obs=%f (%f)\n",
			mod->ts_bs, mod->ts_obs, (double) ntime/(i+1));
	}

	mod->pvalue = ntime / opt->n_bootstrap;

	cleanup_parametric_bootstrap(dat);

	return err;
} /* run_bootstrap */

/**
 * Print state of model after maximize_likelihood().
 *
 * @param opt options object pointer
 * @param dat data object pointer
 * @param mod model object pointer
 * @param diff output of difftime()
 */
void print_model_state(options *opt, data *dat, model *mod, int diff, int newline)
{
	if (opt->compact) {
		fprintf(stdout, "%s %s %s %d %u %e %e %e %e %f %f %f ",
			opt->filename,
			opt->accel_abbreviation,
			opt->admixture?"admix":"mix", mod->K, opt->seed,
			opt->eta_lower_bound, opt->p_lower_bound,
			opt->abs_error, opt->rel_error,
			mod->max_logL, aic(mod), bic(dat, mod));
		if (opt->afile)
			fprintf(stdout, "%f ", mod->arand);
		else
			fprintf(stdout, "ND ");
		fprintf(stdout, "%s %02d:%02d:%02d %d %d %d %d",
			mod->ever_converged?"converged":"not",
			(int)(diff/3600),
			(int)((diff%3600)/60), (diff%60),
			mod->n_total_iter, 
			mod->n_init, mod->n_maxll_init,
			mod->n_maxll_times);
		if (opt->target_ll)
			fprintf(stdout, " %f %d %d", opt->desired_ll,
				mod->n_targetll_init,
				mod->n_targetll_times);
		if (mod->time_stop)
			fprintf(stdout, " time");
		if (newline)
			fprintf(stdout, "\n");
	} else {
		fprintf(stdout, "Dataset: %s\n", opt->filename);
		fprintf(stdout, "Method/Model: %s, %s, K=%d\n",
			opt->accel_abbreviation,
			opt->admixture?"admix":"mix", mod->K);
		fprintf(stdout, "Convergence: ae=%e, re=%e\n",
			opt->abs_error, opt->rel_error);
		fprintf(stdout, "Bounds: e=%e, p=%e\n",
			opt->eta_lower_bound, opt->p_lower_bound);
		fprintf(stdout, "Total number of iterations: %d\n",
			mod->n_total_iter);
		fprintf(stdout, "Total time: %02d:%02d:%02d\n",
			(int)(diff/3600), (int)((diff%3600)/60),
			(diff%60));
		fprintf(stdout, "Iteration of max log likelihood: %d "
			"of %d\n", mod->n_maxll_init, mod->n_init);
		fprintf(stdout, "Number of times reach max log "
			"likelihood: %d\n", mod->n_maxll_times);
		fprintf(stdout, "Maximum log likelihood: %f\n",
			mod->max_logL);
		fprintf(stdout, "AIC: %f\n", aic(mod));
		fprintf(stdout, "BIC: %f\n", bic(dat, mod));
		fprintf(stdout, "Converged: %s\n",
			mod->ever_converged?"yes":"no");
		if (opt->target_ll && mod->n_targetll_times) {
			fprintf(stdout, "Iteration of target log "
				"likelihood (%f): %d\n",
				opt->desired_ll,
				mod->n_targetll_init);
			fprintf(stdout, "Number of times reach target "
				"log likelihood (%f): %d\n",
				opt->desired_ll,
				mod->n_targetll_times);
		} else if (opt->target_ll && !opt->target_revisit)
			fprintf(stdout, "WARNING: Did not reach target log likelihood (%f).\n", opt->desired_ll);
		if (opt->target_revisit && opt->target_ll && mod->n_targetll_times < opt->target_revisit)
			fprintf(stdout, "WARNING: Did not reach target log likelihood (%f) %d times\n", opt->desired_ll, opt->target_revisit);
		else if (opt->target_revisit && !opt->target_ll && mod->n_maxll_times < opt->target_revisit)
			fprintf(stdout, "WARNING: Did not reach max. log likelihood %d times\n", opt->target_revisit);
		if (mod->time_stop)
			fprintf(stdout, "WARNING: Fitting stopped because ran out of time\n");
	}
} /* print_model_state */


/**
 * Synchronize options, data, and model.  This function checks to make sure
 * that the user has made self-consistent choices.  The user cannot fit K
 * subpopulations when there are fewer than K observations.  The current
 * bootstrap hypothesis testing procedure reads -k <k> from the command line
 * and tests H0: K = <k> - 1 vs. HA: K = <k>.  This choice is checked and set
 * up here.
 *
 * @param opt options object
 * @param dat data object
 * @param mod model object
 * @return error status
 */
int synchronize(options *opt, data *dat, model *mod)
{
	int err = NO_ERROR;

	/* [KSD TODO: there is really not much thinking in this...] */
	opt->lower_bound = MIN(opt->lower_bound, 1.0 / dat->I / dat->ploidy - 0.5 / dat->I / dat->ploidy);
	opt->eta_lower_bound = opt->lower_bound;
	opt->p_lower_bound = opt->lower_bound;

	/* no backtracking for quasi-Newton acceleration scheme */
	if (opt->accel_scheme >= QN) {
		opt->adjust_step = 0;
		opt->q = opt->accel_scheme - SQS3;
		if (opt->accel_abbreviation)
			free(opt->accel_abbreviation);
		if (opt->accel_name)
			free(opt->accel_name);
		MAKE_1ARRAY(opt->accel_abbreviation, 4 + (int) log10(opt->q));
		MAKE_1ARRAY(opt->accel_name, strlen(accel_method_names[QN]) + 7 + (int) log10(opt->q));
		sprintf(opt->accel_abbreviation, "Q%d", opt->q);
		sprintf(opt->accel_name, "%s (q=%d)", accel_method_names[QN], opt->q);
		if (mod->Ainv)
			FREE_1ARRAY(mod->Ainv);
		if (mod->cutu)
			FREE_1ARRAY(mod->cutu);
		MAKE_1ARRAY(mod->Ainv, opt->q*opt->q);
		MAKE_1ARRAY(mod->cutu, opt->q);
#ifdef LAPACK
		if (opt->q > 3) {
			MAKE_1ARRAY(mod->ipiv, opt->q + 1);
			MAKE_1ARRAY(mod->work, opt->q*opt->q);
		} else
			MAKE_1ARRAY(mod->A, opt->q * opt->q);
#else
		if (mod->A)
			FREE_1ARRAY(mod->A);
		if (opt->q > 3)
			return message(stderr, __FILE__, __func__, __LINE__,
				ERROR_MSG, INVALID_USER_SETUP, "Cannot use "
				"acceleration methods greater than 6 (QN3) "
				"without linking to lapack.\n");
		MAKE_1ARRAY(mod->A, opt->q * opt->q);
#endif
	} else {
		MAKE_1ARRAY(opt->accel_abbreviation, strlen(accel_method_abbreviations[opt->accel_scheme]) + 1);
		MAKE_1ARRAY(opt->accel_name, strlen(accel_method_names[opt->accel_scheme]) + 1);
		strcpy(opt->accel_abbreviation, accel_method_abbreviations[opt->accel_scheme]);
		strcpy(opt->accel_name, accel_method_names[opt->accel_scheme]);
	}

	if (dat->I < opt->max_K)
		return message(stderr, __FILE__, __func__, __LINE__, ERROR_MSG,
			INVALID_USER_SETUP, "Maximum number of clusters (%d) "
			"(set with command-line argument -k) cannot exceed "
			"the number of individuals (%d)\n", opt->max_K, dat->I);
	if (opt->n_bootstrap && opt->max_K <= 1)
		return message(stderr, __FILE__, __func__, __LINE__, ERROR_MSG,
			INVALID_USER_SETUP, "When bootstrapping, maximum K (%d) "
			"(set with command-line argument -k) must exceed 1.",
			opt->max_K);
	if (opt->n_bootstrap) {
		mod->null_K = opt->max_K - 1;
		mod->alt_K = opt->max_K;
	}

	/* if not timed or target seeking, then run once */
	if (!opt->target_ll && !opt->target_revisit && !opt->n_seconds
		&& !opt->n_init)
		opt->n_init = 1;

	if (opt->min_K > opt->max_K)
		return message(stderr, __FILE__, __func__, __LINE__, ERROR_MSG,
			INVALID_USER_SETUP, "Minimum K (%d) must not exceed "
			"maximum K (%d).", opt->min_K, opt->max_K);
	
	if (opt->afile)
		err = read_afile(opt, dat);

	return err;
} /* synchronize */


/**
 * Allocate options object and initialize.
 *
 * @param opt options object pointer reference
 * @return error status
 */
int make_options(options **opt)
{

	*opt = malloc(sizeof **opt);
	if (*opt == NULL)
		return message(stderr, __FILE__, __func__, __LINE__, ERROR_MSG,
			MEMORY_ALLOCATION, "options object");
	
	(*opt)->filename = NULL;
	(*opt)->R_format = 0;
	(*opt)->alleles_are_indices = 0;
	(*opt)->seed = 1234567;
	(*opt)->n_init = 50;

	/* convergence criterion set to match Lange's definition */
	(*opt)->max_iter = 0;	/*2000;*/
	(*opt)->rel_error = 0;	/*1e-6;*/
	(*opt)->abs_error = 1e-4;

	/* default tests K=6 */
	(*opt)->min_K = 6;
	(*opt)->max_K = 6;

	/* default initialization method: */
	(*opt)->initialization_method = RANDOM_CENTERS;	// TESTING;
	(*opt)->initialization_procedure = NOTHING;/*RAND_EM;*/
	(*opt)->n_rand_em_init = 50;
	(*opt)->lower_bound = 1e-8;
	(*opt)->eta_lower_bound = 1e-8;
	(*opt)->p_lower_bound = 1e-8;
	(*opt)->path = "./";
	(*opt)->admixture = 0;
	(*opt)->eta_constrained = 0;
	(*opt)->n_bootstrap = 0;
	(*opt)->block_relax = 0;
	(*opt)->accel_scheme = 0;
	(*opt)->accel_name = NULL;
	(*opt)->accel_abbreviation = NULL;
	(*opt)->q = 1;
	(*opt)->n_init_iter = 0;
	(*opt)->n_seconds = 0;
	(*opt)->adjust_step = 0;//1000;	/* Varadhan2008: \infty */
	(*opt)->verbosity = MINIMAL;
	(*opt)->compact = 1;
	(*opt)->qfile = NULL;
	(*opt)->pfile = NULL;
	(*opt)->afile = NULL;
	(*opt)->target_ll = 0;
	(*opt)->desired_ll = 0;
	(*opt)->target_revisit = 0;
	(*opt)->n_repeat = 1;
	(*opt)->repeat_seconds = 0;
	(*opt)->max_repeat_seconds = 0;
	(*opt)->write_files = 1;
	(*opt)->partition_from_file = NULL;
	(*opt)->parallel = 0;

	return NO_ERROR;
} /* make_options */


/**
 * Free options object.
 *
 * @param opt options object
 * @return void
 */
void free_options(options *opt)
{
	if (opt->partition_from_file)
		FREE_VECTOR(opt->partition_from_file);
	if (opt->accel_name)
		free(opt->accel_name);
	if (opt->accel_abbreviation)
		free(opt->accel_abbreviation);
	if (opt)
		free(opt);
	opt = NULL;
} /* free_options */


/**
 * Allocate data object and initialize.
 *
 * @param dat data object pointer reference
 * @return error status
 */
int make_data(data **dat)
{

	*dat = malloc(sizeof **dat);
	if (*dat == NULL)
		return message(stderr, __FILE__, __func__, __LINE__, ERROR_MSG,
			MEMORY_ALLOCATION, "data object");

	(*dat)->I = 0;
	(*dat)->L = 0;
	(*dat)->M = 0;
	(*dat)->ploidy = 2;
	(*dat)->IL = NULL;
	(*dat)->ila = NULL;
	(*dat)->uniquealleles = NULL;
	(*dat)->L_alleles = NULL;
	(*dat)->ILM = NULL;
	(*dat)->bs_ILM = NULL;
	(*dat)->pops = NULL;
	(*dat)->i_p = NULL;
	(*dat)->max_M = 0;
	(*dat)->missing_data = 0;
	(*dat)->I_K = NULL;
	(*dat)->IL_K = NULL;

	return NO_ERROR;
} /* make_data */


/**
 * Free data object.
 *
 * @param dat options object
 * @return void
 */
void free_data(data *dat)
{
	int i;
	if (dat->I_K)
		FREE_VECTOR(dat->I_K);
	if (dat->idv) {
		for (i = 0; i < dat->I; i++)
			if (dat->idv[i].name)
				free(dat->idv[i].name);
		FREE_VECTOR(dat->idv);
	}
	if (dat->IL_K)
		FREE_2ARRAY(dat->IL_K);
	if (dat->IL)
		FREE_2ARRAY(dat->IL);
	if (dat->uniquealleles)
		FREE_1ARRAY(dat->uniquealleles);
	if (dat->L_alleles)
		FREE_2ARRAY(dat->L_alleles);
	if (dat->ILM)
		FREE_3ARRAY(dat->ILM);
	if (dat->bs_ILM)
		FREE_3ARRAY(dat->bs_ILM);
	if (dat->pops) {
		for (i=0; i<dat->numpops; i++)
			free(dat->pops[i]);
		free(dat->pops);
		dat->pops = NULL;
	}
	if (dat->i_p)
		FREE_1ARRAY(dat->i_p);

	if (dat)
		free(dat);
	dat = NULL;
} /* free_data */


/**
 * Allocate model object and initialize.
 *
 * @param mod model data object pointer reference
 * @return error status
 */
int make_model(model **mod)
{
	*mod = malloc(sizeof **mod);
	if (*mod == NULL)
		return message(stderr, __FILE__, __func__, __LINE__, ERROR_MSG,
			MEMORY_ALLOCATION, "model object");
	
	(*mod)->K = 1;
#ifndef OLDWAY
	(*mod)->vpklm = NULL;
	(*mod)->vetaik = NULL;
	(*mod)->vetak = NULL;
	(*mod)->pindex = 0;	/* default behavior: unchanged */
	(*mod)->findex = 0;	/* default behavior: unchanged */
	(*mod)->tindex = 0;	/* default behavior: unchanged */
#else
	(*mod)->pKLM = NULL;
	(*mod)->etak = NULL;
	(*mod)->etaik = NULL;
#endif
	(*mod)->mle_pKLM = NULL;
	(*mod)->mle_etak = NULL;
	(*mod)->mle_etaik = NULL;
	(*mod)->vik = NULL;
	(*mod)->diklm = NULL;

#ifndef OLDWAY
	(*mod)->u_etaik = NULL;
	(*mod)->u_pklm = NULL;
	(*mod)->u_etak = NULL;
	(*mod)->v_etaik = NULL;
	(*mod)->v_pklm = NULL;
	(*mod)->v_etak = NULL;
	(*mod)->A = NULL;
	(*mod)->delta_index = 0;
	(*mod)->uindex = 0;
	(*mod)->vindex = 0;
#else
	(*mod)->U = NULL;
	(*mod)->V = NULL;
	(*mod)->init_etaik = NULL;
	(*mod)->iter1_etaik = NULL;
	(*mod)->iter2_etaik = NULL;
	(*mod)->init_pKLM=NULL;
	(*mod)->iter1_pKLM = NULL;
	(*mod)->iter2_pKLM = NULL;
	(*mod)->init_etak = NULL;
	(*mod)->iter1_etak = NULL;
	(*mod)->iter2_etak = NULL;
#endif
	(*mod)->A = NULL;
	(*mod)->Ainv = NULL;
	(*mod)->cutu = NULL;
#ifdef LAPACK
	(*mod)->ipiv = NULL;
	(*mod)->work = NULL;
#endif

	(*mod)->current_i = 0;
	(*mod)->current_l = 0;
	(*mod)->current_k = 0;
	(*mod)->etaupdate = 0;
	(*mod)->converged = 0;
	(*mod)->max_logL = -Inf;
	(*mod)->first_max_logL = -Inf;
	(*mod)->n_iter = 0;
	(*mod)->n_max_iter = 0;
	(*mod)->n_total_iter = 0;
	(*mod)->accel_abort = 1;
	(*mod)->accel_step = 0;

	(*mod)->n_maxll_times = 1;
	(*mod)->n_targetll_times = 0;
	(*mod)->n_targetll_init = 0;

	(*mod)->time_stop = 0;
	(*mod)->iter_stop = 0;

	return NO_ERROR;
} /* make_model */


/**
 * Allocate model.  Allocation failure leads to immediate return with error
 * code.
 *
 * [TODO] Nicer error handling.
 *
 * @param opt options object
 * @param mod model object
 * @param dat data object
 * @return error status
 */
int allocate_model_for_k(options *opt, model *mod, data *dat)
{

#ifndef OLDWAY
	/* need to keep x^n, while computing F(x^n) - x^n, F(F(X^n)) - F(x^n),
	   ultimately saving F(F(X^n)) for possible fallback */
	MAKE_4JAGGED_ARRAY(mod->vpklm, 3, mod->K, dat->L, dat->uniquealleles);
#else
	MAKE_3JAGGED_ARRAY(mod->pKLM, mod->K, dat->L, dat->uniquealleles);
#endif
	if (opt->n_bootstrap && !mod->mle_pKLM)
		MAKE_3JAGGED_ARRAY(mod->mle_pKLM, mod->K, dat->L, dat->uniquealleles);

	MAKE_1ARRAY(mod->count_K, mod->K);
	
	if (opt->admixture) {
		MAKE_4JAGGED_ARRAY(mod->diklm, dat->I, mod->K, dat->L,
			dat->uniquealleles);
		if (opt->eta_constrained) {
#ifndef OLDWAY
			MAKE_2ARRAY(mod->vetak, 3, mod->K);
#else
			MAKE_1ARRAY(mod->etak, mod->K);
#endif
			if (opt->n_bootstrap && !mod->mle_etak)
				MAKE_1ARRAY(mod->mle_etak, mod->K);
		} else {
#ifndef OLDWAY
			MAKE_3ARRAY(mod->vetaik, 3, dat->I, mod->K);
#else
			MAKE_2ARRAY(mod->etaik, dat->I, mod->K);
#endif
			if (opt->n_bootstrap && !mod->mle_etaik)
				MAKE_2ARRAY(mod->mle_etaik, dat->I, mod->K);
		}
				
	} else {
		MAKE_2ARRAY(mod->vik, dat->I, mod->K);
#ifndef OLDWAY
		MAKE_2ARRAY(mod->vetak, 3, mod->K);
#else
		MAKE_1ARRAY(mod->etak, mod->K);
#endif
		if (opt->n_bootstrap && !mod->mle_etak)
			MAKE_1ARRAY(mod->mle_etak, mod->K);
	}
	if (opt->accel_scheme) {
#ifndef OLDWAY
		MAKE_4JAGGED_ARRAY(mod->u_pklm, opt->q, mod->K, dat->L,
			dat->uniquealleles);
		MAKE_4JAGGED_ARRAY(mod->v_pklm, opt->q, mod->K, dat->L,
			dat->uniquealleles);
#else
		MAKE_3JAGGED_ARRAY(mod->init_pKLM, mod->K, dat->L,
			dat->uniquealleles);
		MAKE_3JAGGED_ARRAY(mod->iter1_pKLM, mod->K, dat->L,
			dat->uniquealleles);
		MAKE_3JAGGED_ARRAY(mod->iter2_pKLM, mod->K, dat->L,
			dat->uniquealleles);
		MAKE_1ARRAY(mod->U, dat->max_M);
		MAKE_1ARRAY(mod->V, dat->max_M);
#endif
		if (opt->admixture && !opt->eta_constrained) {
#ifndef OLDWAY
			MAKE_3ARRAY(mod->u_etaik, opt->q, dat->I, mod->K);
			MAKE_3ARRAY(mod->v_etaik, opt->q, dat->I, mod->K);
#else
			MAKE_2ARRAY(mod->init_etaik, dat->I, mod->K);
			MAKE_2ARRAY(mod->iter1_etaik, dat->I, mod->K);
			MAKE_2ARRAY(mod->iter2_etaik, dat->I, mod->K);
#endif
		} else {
#ifndef OLDWAY
			MAKE_2ARRAY(mod->u_etak, opt->q, mod->K);
			MAKE_2ARRAY(mod->v_etak, opt->q, mod->K);
#else
			MAKE_1ARRAY(mod->init_etak, mod->K);
			MAKE_1ARRAY(mod->iter1_etak, mod->K);
			MAKE_1ARRAY(mod->iter2_etak, mod->K);
#endif
		}
	}

	/* count number of parameters */
	if (!opt->admixture || opt->eta_constrained)
		mod->no_parameters = (mod->K - 1);
	else
		mod->no_parameters = dat->I * (mod->K - 1);

	for (int l = 0; l < dat->L; l++) {
		int m_start = dat->L_alleles && dat->L_alleles[l][0] == MISSING;
		mod->no_parameters += (dat->uniquealleles[l] - 1 - m_start) * mod->K;
	}

	return NO_ERROR;
} /* allocate_model_for_k */


/**
 * Free model object.  
 *
 * @param mod model object
 * @return void
 */
void free_model_data(model *mod, options *opt)
{
#ifndef OLDWAY
	FREE_4ARRAY(mod->vpklm);
#else
	FREE_3ARRAY(mod->pKLM);
#endif
	
	FREE_1ARRAY(mod->count_K);

	if (opt->admixture)
		FREE_4ARRAY(mod->diklm);
	if (!opt->admixture)
		FREE_2ARRAY(mod->vik);

#ifndef OLDWAY
	FREE_4ARRAY(mod->u_pklm);
	FREE_4ARRAY(mod->v_pklm);
#else
	FREE_3ARRAY(mod->init_pKLM);
	FREE_3ARRAY(mod->iter1_pKLM);
	FREE_3ARRAY(mod->iter2_pKLM);
	FREE_1ARRAY(mod->U);
	FREE_1ARRAY(mod->V);
#endif

	if (opt->admixture && !opt->eta_constrained) {
#ifndef OLDWAY
		FREE_3ARRAY(mod->u_etaik);
		FREE_3ARRAY(mod->v_etaik);
		FREE_3ARRAY(mod->vetaik);
#else
		FREE_2ARRAY(mod->init_etaik);
		FREE_2ARRAY(mod->iter1_etaik);
		FREE_2ARRAY(mod->iter2_etaik);
		FREE_2ARRAY(mod->etaik);
#endif
	} else {
#ifndef OLDWAY
		FREE_2ARRAY(mod->u_etak);
		FREE_2ARRAY(mod->v_etak);
		FREE_2ARRAY(mod->vetak);
#else
		FREE_1ARRAY(mod->init_etak);
		FREE_1ARRAY(mod->iter1_etak);
		FREE_1ARRAY(mod->iter2_etak);
		FREE_1ARRAY(mod->etak);
#endif
	}
} /* free_model_data */

/**
 * Free model MLEs. MLEs are collected over multiple initializations, so the
 * memory to store them are not released until the user requests it.
 *
 * @param mod model object
 */
void free_model_mles(model *mod)
{
	if (mod->mle_pKLM)
		FREE_3ARRAY(mod->mle_pKLM);
	if (mod->mle_etak)
		FREE_1ARRAY(mod->mle_etak);
	if (mod->mle_etaik)
		FREE_2ARRAY(mod->mle_etaik);
} /* free_model_mles */

/**
 * Free model object.
 *
 * @param mod model object
 * @return void
 */
void free_model(model *mod, options *opt)
{
	if (mod) {
		free_model_data(mod,opt);
		free_model_mles(mod);
		if (mod->Ainv)
			FREE_1ARRAY(mod->Ainv);
		if (mod->A)
			FREE_1ARRAY(mod->A);
		if (mod->cutu)
			FREE_1ARRAY(mod->cutu);
#ifdef LAPACK
		if (mod->work)
			FREE_1ARRAY(mod->work);
		if (mod->ipiv)
			FREE_1ARRAY(mod->ipiv);
#endif
		free(mod);
	}
	mod = NULL;
} /* free_model */

/**
 * Parse command line.
 *
 * The return value contains strings that should not be freed until program
 * exit.  See print_usage() for command-line options and usage.
 *
 * @param opt options object
 * @param dat data object
 * @param argc number of command-line arguments
 * @param argv command-line arguments
 * @return error status
 */
int parse_options(options *opt, data *dat, int argc, const char **argv)
{
	int i, j;
	int err = NO_ERROR;
	char a;

	for (i=1; i<argc; i++) {
		if (strlen(argv[i]) < 2)
			usage_error(argv, i, (void *)opt);

		/* skip to argument name */
		j = 1;
		a = argv[i][j];
		while (a == '-' && ++j < (int) strlen(argv[i]))
			a = argv[i][j];

		switch (a) {
			case 'a':
				opt->admixture = 1;
				break;
			case 'A':
				opt->afile = argv[++i];
				break;
			case 'b':
				opt->n_bootstrap = read_int(argc, argv, ++i,
					(void *)opt);
				if (opt->n_bootstrap < 0 || errno
					|| opt->block_relax)
					goto CMDLINE_ERROR;
				break;
			case 'c':
				opt->eta_constrained = 1;
				break;
			case 'd':
				opt->path = argv[++i];
				break;
			case 'e':
				opt->rel_error = read_double(argc, argv, ++i,
					(void *)opt);
				if (opt->rel_error < 0 || errno)
					goto CMDLINE_ERROR;
				break;
			case 'E':
				opt->abs_error = read_double(argc, argv, ++i,
					(void *)opt);
				if (opt->abs_error < 0 || errno)
					goto CMDLINE_ERROR;
				break;
			case 'f':
				opt->filename = argv[++i];
				break;
			case 'g':
				opt->adjust_step = read_int(argc, argv, ++i,
					(void *)opt);
				if (opt->adjust_step < 0 || errno)
					goto CMDLINE_ERROR;
				break;
			case 'h':
				fprint_usage(stdout, argv[0], (void *)opt);
				return CUSTOM_ERROR;
			case 'i':
				opt->n_init_iter = read_int(argc, argv, ++i,
					(void *)opt);
				if (opt->n_init_iter < 0 || errno)
					goto CMDLINE_ERROR;
				break;
			case 'I':
				opt->alleles_are_indices = 1;
				break;
			case '1':
				opt->min_K = read_int(argc, argv, ++i,
					(void *)opt);
				if (opt->min_K < 1 || errno)
					goto CMDLINE_ERROR;
				break;
			case '2':
				opt->max_K = read_int(argc, argv, ++i,
					(void *)opt);
				if (opt->max_K < 1 || errno)
					goto CMDLINE_ERROR;
				break;
			case 'k':
				opt->max_K = read_int(argc, argv, ++i,
					(void *)opt);
				if (opt->max_K < 1 || errno)
					goto CMDLINE_ERROR;
				opt->min_K = opt->max_K;
				break;
			case 'm':
				opt->n_rand_em_init = read_int(argc, argv, ++i,
					(void *)opt);
				if (opt->n_rand_em_init == 0)
					opt->initialization_procedure = NOTHING;
				if (opt->n_rand_em_init < 0 || errno)
					goto CMDLINE_ERROR;
				break;
			case 'M':
				if (opt->parallel == 1)
					opt->n_repeat = 1;
					opt->verbosity = SILENT;
			case 'n':
				opt->n_init = read_int(argc, argv, ++i,
					(void *)opt);
				if (errno)
					goto CMDLINE_ERROR;
				break;
			case 'p':
				dat->ploidy = read_int(argc, argv, ++i,
					(void *)opt);
				if (dat->ploidy < 1 || errno)
					goto CMDLINE_ERROR;
				break;
			case 'P':
				opt->pfile = argv[++i];
				break;
			case 'Q':
				opt->qfile = argv[++i];
				break;
			case 'R':
				opt->R_format = 1;
				break;
			case 'r':
				opt->seed = read_uint(argc, argv, ++i,
					(void *)opt);
				srand(opt->seed);
				break;
			case 'x':
				opt->block_relax = 1;
				if (opt->n_bootstrap > 0 || errno)
					goto CMDLINE_ERROR;
				//if ( || errno)
				//	goto CMDLINE_ERROR;
				break;
			case 's':
				opt->accel_scheme = read_int(argc, argv, ++i,
					(void *)opt);
				if (opt->accel_scheme < 0 || errno)
					goto CMDLINE_ERROR;
				break;
			case 't':
				opt->n_seconds = 60*read_uint(argc, argv, ++i,
					(void *)opt);
				break;
			case 'T':
				opt->max_iter = read_int(argc, argv, ++i,
					(void *)opt);
				if (opt->max_iter < 0 || errno)
					goto CMDLINE_ERROR;
				break;
			case 'u':
				while (++i < argc && argv[i][0] != '-') {
					switch (argv[i][0]) {
						case 'l':
							opt->target_ll = 1;
							opt->desired_ll =
								read_double(
								argc, argv, ++i,
								(void *)opt);
							break;
						case 'n':
							opt->target_revisit = 
								read_int(argc,
								argv, ++i,
								(void *)opt);
							if (opt->target_revisit
								< 0 || errno)
								goto CMDLINE_ERROR;
							break;
						default:
							err = INVALID_CMD_OPTION;
							goto CMDLINE_ERROR;
					}
				}
				i--;
				break;
			case 'v':
				if (i + 1 == argc)
					opt->verbosity = VERBOSE;
				else {
					char *ret_ptr;
					opt->verbosity = strtol(argv[++i], &ret_ptr, 0);
					if (errno)
						usage_error(argv, i, (void *)opt);
					else if (ret_ptr == argv[i]) {
						opt->verbosity = VERBOSE;
						i--;
					}
				}
				break;
			case 'w':
				while (++i < argc && argv[i][0] != '-') {
					switch (argv[i][0]) {
						case 't':
							opt->repeat_seconds =
								60*read_int(
								argc, argv, ++i,
								(void *)opt);
							break;
						case 'm':
							opt->max_repeat_seconds =
								60*read_int(
								argc, argv, ++i,
								(void *)opt);
							break;
						case 'n':
							opt->n_repeat = 
								read_int(argc,
								argv, ++i,
								(void *)opt);
							if (opt->n_repeat
								<= 0 || errno)
								goto CMDLINE_ERROR;
							break;
						default:
							err = INVALID_CMD_OPTION;
							goto CMDLINE_ERROR;
					}
				}
				i--;
				opt->write_files = 0;
				break;
			default:
				err = INVALID_CMD_OPTION;
				goto CMDLINE_ERROR;
		}
	}

	if (opt->filename == NULL)
		return message(stderr, __FILE__, __func__, __LINE__, ERROR_MSG,
			INVALID_CMDLINE, "You must specify the data file "
			"with command line option '-f'.  Try '-h' for help.\n");

	return NO_ERROR;

CMDLINE_ERROR:
	if (err == NO_ERROR) {
		err = INVALID_CMD_ARGUMENT;
		i--;
	}
	usage_error(argv, i, (void *)opt);
	return err;
} /* parse_options */

/**
 * Print command-line usage.
 *
 * @param fp file stream to print on
 * @param invocation_name name of command
 * @return void
 */
void fprint_usage(FILE *fp, const char *invocation_name, void *obj)
{
	options *opt = obj ? (options *) obj : NULL;
	/* strip the "./" from the beginning of the program invocation */
	const char *prog_name = strlen(invocation_name) > 2
		&& invocation_name[0] == '.' ? &(invocation_name[2])
		: invocation_name;
	fprintf(fp, "\nNAME\n");
	fprintf(fp, "\t%s - Maximum likelihood clustering of discrete data\n",
		prog_name);
	fprintf(fp, "\nSYNOPSIS\n");
	fprintf(fp,
	"\t%s [-k <n> | -1 <n> -2 <n>] [-a -b <n> -c -C <n> -d <s> -e <d> -E <d> -f <d> -g <d> -h"
	"\n\t\t-i <n> -I -m <n> -n <n> -p <n> -R -s <n> -t <n> -T <d> -u <s> -v -w <s> -x] -f <s>\n"
	"\n\t\twhere <n> stands for integer, <s> for string, <d> for double",
		prog_name);
	fprintf(fp, "\nDESCRIPTION\n");
	fprintf(fp, 
	"\t%s clusters multivariate discrete data observed on a sample of\n"
	"\tindividuals using the EM algorithm.  It handles data missing at\n"
	"\trandom.  It assumes coordinates within an individual are independent.\n"
	"\tIt allows the admixture model, where each coordinate is independently\n"
	"\tdrawn from a cluster, or the mixture model, where each individual is\n"
	"\tdrawn from a cluster.\n", prog_name);
	fprintf(fp, "\nOPTIONS\n"
		/* ---------------------------------------------------------- */
		"\t-a\t"
		"Choose admixture model (default: %s).\n"
		"\t-b\t"
		"Bootstrap test of H0: K=<k>-1 vs. Ha: K=<k>, where <k> is\n"
		"\t\tgiven by -k option.  Specify number of bootstraps as\n"
		"\t\targument (default: %d).\n"
		"\t-c\t"
		"Constrain mixing proportions identical across individuals\n"
		"\t\t(only enforced with -a; default: %s).\n"
		"\t-C\t"
		"The maximum number of iterations to fit (default: %d).\n"
		"\t-d\t"
		"Directory where output files are written (default: %s).\n"
		"\t-e\t"
		"Allowable log likelihood relative error for convergence\n"
		"\t\t(default: %.1e).\n"
		"\t-E\t"
		"Allowable log likelihood absolute error for convergence\n"
		"\t\t(default: %.1e).\n"
		"\t-f\t"
		"Name of data file (STRUCTURE format).\n"
		"\t-g\t"
		"Adjust step size at most this many times (default: %d)\n"
		"\t-h\t"
		"This help.\n"
		"\t-i\t"
		"Initial iterations prior to acceleration (default: %d)\n"
		"\t-I\t"
		"Alleles are indices (no sorting, etc.) (default: %s)\n"
		"\t-k\t"
		"The number of clusters to fit (default: %d).\n"
		"\t-1\t"
		"The minimum number of clusters to fit (default: %d).\n"
		"\t-2\t"
		"The maximum number of clusters to fit (default: %d).\n",
		opt->admixture?"yes":"no", opt->n_bootstrap, 
		opt->eta_constrained?"yes":"no", opt->max_iter,
		opt->path, opt->rel_error, opt->abs_error,
		opt->adjust_step,
		opt->n_init_iter, opt->alleles_are_indices?"yes":"no",
		opt->max_K, opt->min_K, opt->max_K
	);
	fprintf(fp, 
		"\t-m\t"
		"The number of Rand EM initializations, 0 to avoid Rand EM\n"
		"\t\t(default: %d).\n"
		"\t-M\t"
		"Parallelization option\n"
		"\t\t(default: no parallelization). \n"
		"\t-n\t"
		"Number of initializations to run EM to convergence\n"
		"\t\t(default: %d).\n"
		"\t-p\t"
		"The ploidy (default: 2).\n"
		"\t-r\t"
		"Random number (default: %u).\n"
		"\t-R\t"
		"Data file in R format (default: %s).\n"
		"\t-s\t"
		"The acceleration scheme (default: %s).\n"
		"\t\tAcceleration options:\n"
		"\t\t\t0 (default) - no acceleration\n"
		"\t\t\t1 - SQUAREM version 1\n"
		"\t\t\t2 - SQUAREM version 2\n"
		"\t\t\t3 - SQUAREM version 3\n"
		"\t\t\t4 - Quasi Newton version 1 (1 secant condition)\n"
		"\t\t\t5 - Quasi Newton version 2 (2 secant conditions)\n"
		"\t\t\t6 - Quasi Newton version 3 (3 secant conditions)\n",
		opt->n_rand_em_init, opt->n_init, opt->seed,
		opt->R_format?"yes":"no", accel_method_names[MIN(opt->accel_scheme, QN)]
	);
	if (opt->accel_scheme>=QN)
		fprintf(fp, "q=%d", opt->accel_scheme-SQS3);
	fprintf(fp,
		"\t-u\t"
		"Iterate until beat target\n"
		"\t\t-u n #: repeat until reach same max # times (default: %d)\n"
		"\t\t-u l #: repeat until reach max log likelihood # (default: %f)\n"
		"\t-t\t"
		"The time (in minutes) to maximize likelihood (default: %d).\n"
		"\t\tBe sure to check convergence if you set the above!\n"
		"\t-v\t"
		"Level of verbosity (default: %d)\n"
		"\t-w\t"
		"Write best fits to file (default: %s)\n"
		"\t\tAlso used to time repeated runs, as follows:\n"
		"\t\t-w n <n>: repeat at least <n> times (default: %d)\n"
		"\t\t-w t <n>: repeat at least <n> minutes (default: %d)\n"
		"\t\t-w m <n>: repeat at most <n> minutes (default: %d)\n"
		"\t-z\t"
		"Use block relaxation algorithm (default: %s).\n",
		opt->target_revisit, opt->target_ll?opt->desired_ll:0,
		opt->n_seconds, opt->verbosity, 
		opt->write_files ? "yes" : "no",
		opt->n_repeat, opt->repeat_seconds,
		opt->max_repeat_seconds, opt->block_relax?"yes":"no"
	);
} /* fprintf_usage */

/**
 * Compute (Adjusted) Rand index or E index.
 * 
 * @param n number of observations
 * @param k1 true number of clusters
 * @param k2 estimated number of clusters
 * @param cl1 true classification vector
 * @param cl2 estimated classification vector
 * @return Rand index, Adjusted Rand index, or E index
 */
double adj_rand(int n, int k1, int k2, int *cl1, int *cl2, int type)
{
	int i, j, **nmat;
	double *sumtr, *sumpr;
	double sumprsq, sumtrsq, sumsq, discordant, sumtrprsq;
	double term1, term2, term3;
	double nij2sum, nidot2sum, ndotj2sum;
	double index = 0;

	CMAKE_2ARRAY(nmat, k1, k2);
	MAKE_VECTOR(sumtr, k1);
	MAKE_VECTOR(sumpr, k2);

	for (i = 0; i < n; i++)
		nmat[cl1[i]][cl2[i]]++;
  
	sumtrsq=0.;
	for (i = 0; i < k1; i++) {
		sumtr[i] = 0.; 
		for (j = 0; j < k2; j++)
			sumtr[i] += nmat[i][j];
		sumtrsq += sumtr[i] * sumtr[i];
	}
  
	sumprsq = 0.;
	for (j = 0; j < k2; j++) {
		sumpr[j] = 0.;
		for (i = 0; i < k1; i++)
			sumpr[j] += nmat[i][j];
		sumprsq += sumpr[j] * sumpr[j];
	}

	if (type == E_INDEX) {
		sumtrprsq = 0.;
		for (i = 0;i < k1; i++)
			for (j = 0; j < k2; j++)
				sumtrprsq += sumtr[i] * sumtr[i] * sumpr[j] * sumpr[j];

		index = sumtrprsq / (n*(n-1) + (double) n*n/(n-1))
			- (sumprsq + sumtrsq)/(n-1);
		index *= 2.;
		index /= n*(n-1);
	}

	if (type == RAND_INDEX) {
		sumsq = 0.;
		for (i = 0; i < k1; i++)
			for (j = 0; j < k2; j++)
				sumsq += nmat[i][j] * nmat[i][j];

		discordant = 0.5*(sumtrsq + sumprsq) - sumsq;

		index = 1.0 - discordant / (n * (n-1.) / 2.);

	} else if (type == ADJUSTED_RAND_INDEX) {

		nidot2sum = 0.;
		for (i = 0; i < k1; i++)
			nidot2sum += sumtr[i] * (sumtr[i]-1) / 2;

		nij2sum = 0.;
		for (i = 0; i < k1; i++)
			for (j = 0; j < k2; j++)
				nij2sum += nmat[i][j] * (nmat[i][j]-1) / 2.0;

		ndotj2sum = 0.;
		for (i = 0; i < k2; i++)
			ndotj2sum += sumpr[i] * (sumpr[i]-1) / 2.0;
	

		term3 = nidot2sum * ndotj2sum / (n * (n-1.) / 2.);
		term1 = nij2sum - term3;
		term2 = (nidot2sum + ndotj2sum)/2 - term3;
		index = term1 / term2;
	}

	FREE_VECTOR(sumtr);
	FREE_VECTOR(sumpr);
	FREE_MATRIX(nmat);

	return index;
} /* adj_rand */
