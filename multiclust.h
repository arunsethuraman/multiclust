/**
 * @file multiclust.h
 * @author Arun Sethurman
 * @author Karin Dorman, kdorman@iastate.edu
 *
 * Header file for multiclust structs, functions and defines.
 */

#ifndef __H_MULTICLUST__
#define __H_MULTICLUST__


#ifdef __HAVE_R_
        #include <R.h>
	#include <Rmath.h>
	#undef printf
	#define printf Rprintf
#else
	#include <stdio.h>
	#include <stdlib.h>
	#include <math.h>
	#include <float.h>
#endif

#define e_step_admixture(A, B, C) e_step_admixture_orig(A, B, C)
#define m_step_admixture(A, B, C) m_step_admixture_orig(A, B, C)

/* quiet unused parameter warnings */
#ifdef UNUSED
#elif defined(__GNUC__)
# define UNUSED(x) UNUSED_ ## x __attribute__((unused))
#elif defined(__LCLINT__)
# define UNUSED(x) /*@unused@*/ x
#else
# define UNUSED(x) x
#endif

/* min and max macros */
#define MAX(a,b) ((a)>(b)?(a):(b))
#define MIN(a,b) ((a)<(b)?(a):(b))

/* memory allocation macros */
#include <string.h>
#include <time.h>
#include <errno.h>
#include "message.h"
#include "array.h"

/* new macros to redefine array.h error handling */
/* note: *_ERR versions are not entirely thread-safe */
#define MAKE_1ARRAY_RETURN(a,n) do {                                           \
	(a) = malloc((n) * sizeof *(a));                                       \
	if ((a) == NULL)                                                       \
		return message(stderr, __FILE__, __func__, __LINE__,           \
			ERROR_MSG, MEMORY_ALLOCATION, NULL);                   \
} while (0)
#define MAKE_1ARRAY_ERR(a,n) do {                                              \
	(a) = malloc((n) * sizeof *(a));                                       \
	if ((a) == NULL)                                                       \
		errno = message(stderr, __FILE__, __func__, __LINE__,          \
			ERROR_MSG, MEMORY_ALLOCATION, NULL);                   \
} while (0)
#define CMAKE_1ARRAY_RETURN(a,n) do {                                          \
	(a) = calloc((n), sizeof *(a));                                        \
	if ((a) == NULL)                                                       \
		return message(stderr, __FILE__, __func__, __LINE__,           \
			ERROR_MSG, MEMORY_ALLOCATION, NULL);                   \
} while (0)
#define CMAKE_1ARRAY_ERR(a,n) do {                                             \
	(a) = calloc((n), sizeof *(a));                                       \
	if ((a)==NULL)                                                         \
		errno = message(stderr, __FILE__, __func__, __LINE__,          \
			ERROR_MSG, MEMORY_ALLOCATION, NULL);                   \
} while (0)
#define COPY_1ARRAY_RETURN(a,b,n) do {                                         \
	(a) = memcpy((a), (b), (n) * sizeof *(a));                             \
	if ((a) == NULL)                                                       \
		return message(stderr, __FILE__, __func__, __LINE__,           \
			ERROR_MSG, MEMCPY_ERROR, NULL);                        \
} while (0)
#define COPY_1ARRAY_ERR(a,n) do {                                              \
	(a) = memcpy((a), (b), (n) * sizeof *(a));                             \
	if ((a) == NULL)                                                       \
		errno = message(stderr, __FILE__, __func__, __LINE__,          \
			ERROR_MSG, MEMORY_ALLOCATION, NULL);                   \
} while (0)
/* avoid gcc warnings by undef'ing */
#undef MAKE_1ARRAY
#undef CMAKE_1ARRAY
#undef COPY_1ARRAY
/* before reseting to rewritten macros */
#define MAKE_1ARRAY MAKE_1ARRAY_RETURN
#define CMAKE_1ARRAY CMAKE_1ARRAY_RETURN
#define COPY_1ARRAY COPY_1ARRAY_RETURN
//#define LAPACK
/**
 * Initialization methods (how initial parameters are chosen).
 */
enum {	RANDOM_PARTITION,	/*!< random partition of data */
	RANDOM_CENTERS,		/*!< random partition based on random centers */
	TESTING,		/*!< random testing initialization */
	NUM_INITIALIZATION_METHODS	/*!< no. of initialization methods */
};

/**
 * Initialization procedures (what is done with initializations).
 */
enum {	NOTHING,	/*!< no initialization procedure */
	RAND_EM,	/*!< so-called Rand EM */
	NUM_INITIALIZATION_PROCEDURES	/*!< no. of initialization procedures */
};

/**
 * Indices (to compare cluster solutions).
 */
enum {	E_INDEX,		/*!< E index */
	RAND_INDEX,		/*!< RAND index */
	ADJUSTED_RAND_INDEX,	/*!< adjusted RAND index */
	NUM_INDICES		/*!< no. indices */
};

/**
 * Acceleration methods.
 */
enum {	NONE,		/*!< no acceleration */
	SQS1,
	SQS2,
	SQS3,
	QN,
	NUM_ACCELERATION_METHODS
};

/**
 * Output formats.
 */
enum {	STRUCTURE,	/*!< our structure-type format */
	PED,		/*!< Plink's *.ped format */
};

#define MISSING         -9

extern const char *accel_method_abbreviations[NUM_ACCELERATION_METHODS];
extern const char *accel_method_names[NUM_ACCELERATION_METHODS];


typedef struct _options options;
typedef struct _data data;
typedef struct _model model;
typedef struct _indiv indiv;

/**
 * Run options.  Options, many chosen by user, that describe the data and
 * running conditions.
 */
struct _options {
	int initialization_method;	/*!< method of initialization */
	int initialization_procedure;	/*!< procedure of initialization */
	int eta_constrained;		/*!< constrained mixing proportions */
	int admixture;			/*!< admixture model */
	int n_init;			/*!< number of random initializations */
	int n_rand_em_init;		/*!< number of RandEM initializations */
	int max_iter;			/*!< max. number of iterations */
	int min_K;			/*!< min. K to try */
	int max_K;			/*!< max. K to try */
	int missing_value;		/*!< missing data value in input file */
	int imputation_method;		/*!< set optional imputation method for missing data */
	char const *imputed_outfile;	/*!< data file to write with results after imputation */
	int output_format;		/*!< data file output format */
	double lower_bound;		/*!< minimum probability */
	double rel_error;		/*!< relative log like. error for convergence */
	double abs_error;		/*!< absolute log like. error for convergence */
	const char *filename;		/*!< data file name */
	const char *path;		/*!< location of data files */
	const char *admix_pfile;	/*!< name of admixture P file */
	const char *admix_qfile;	/*!< name of admixture Q file */
	int R_format;			/*!< data file in R format */
	int interleaved;		/*!< chromosomes interleaved */
	int alleles_are_indices;	/*!< allele names are indices */
	int one_plus;			/*!< alleles are indices + 1 */
	int write_plus_one;		/*!< when write data, +1 to alleles */
	unsigned int seed;		/*!< random number seed */
	int n_bootstrap;		/*!< no. of bootstrap trials */
	int block_relax;		/*!< block relaxation algorithm */
	int accel_scheme;		/*!< acceleration scheme */
	int do_projection;		/*!< do simplex projection */
	char *accel_name;		/*!< name of chosen acceleration scheme */
	char *accel_abbreviation;	/*!< abbreviation of chosen acceleration scheme */
	int q;				/*!< number of secant approximations (QN) */
	int n_init_iter;		/*!< no. of initial iterations */
	unsigned int n_seconds;		/*!< no. of seconds to run */
	int target_ll;			/*!< stop when reach target log likelihood */
	int target_revisit;		/*!< number of times to reach (same) log likelihood */
	double desired_ll;		/*!< log likelihood to target */
	int adjust_step;		/*!< max. attempts to adjust step size according to Varadhan2008 */
	int verbosity;			/*!< amount of verbosity */
	int compact;			/*!< compact output */
	double eta_lower_bound;		/*!< lower bound on mixing proportions */
	double p_lower_bound;		/*!< lower bound on allele frequencies */
	int n_repeat;			/*!< min. number of times to repeat run */
	unsigned int repeat_seconds;	/*!< min. seconds to repeat run */
	unsigned int max_repeat_seconds;/*!< max. seconds to repeat run */
	int write_files;		/*!< write files with best solution */
	const char *qfile;		/*!< file with estimated eta */
	const char *pfile;		/*!< file with estimated p */
	const char *afile;		/*!< file with partition (true) assignments */
	int *partition_from_file;	/*!< partition from some file */
	int pK;				/*!< K in afile */
	int parallel;			/*!< determines parallelization of code */
	const char *outfile_name; 	/*!< allows user to create unique outfile name */
	int simulate;			/*!< simulate data */
	const char *simulate_outfile;	/*!< simulation data outfile */
}; /* options */

/**
 * The data.  Store information about the genetic data, including the
 * dimensions, the alleles, the unique alleles at each locus, the sufficient
 * statistics, room for bootstrapped data, individual names, and sampling
 * location names.
 */
struct _data {
	/* dimensions of data */
	int I;			/*!< number of sampled individuals */
	int L;			/*!< number of loci sampled */
	int M;			/*!< max. no. distinct alleles at any locus */
	int ploidy;		/*!< ploidy level */
	int missing_data;	/*!< indicate if there is any missing data */

	/* genetic data */
	int **IL;		/*!< original data */
	int ***ila;		/*!< integer transformed data */
	int *uniquealleles;	/*!< no. unique alleles at each locus */
	int **L_alleles;	/*!< list unique alleles at each locus */
	int ***ILM;		/*!< data as unique allele counts */
	int ***bs_ILM;		/*!< stores ILM when bootstrapping */

	/* individuals */
	indiv *idv;		/*!< list of information about individuals */
	int **IL_K;		/*!< partition of alleles into clusters */
	int *I_K;		/*!< partition of individuals into clusters */
	

	/* sampling locales */
	int numpops;		/*!< number of locales */
	char **pops;		/*!< name of physical locales */
	int *i_p;		/*!< sample size per locale */
	int max_M;		/*!< greatest number of alleles across all loci */
}; /* data */

/**
 * The model.  Stores information about the model, including its parameters,
 * maximum likelihood estimates of those parameters, variables for initializing
 * the EM (which are model-specific), and summary statistics obtained from a
 * fitted model.
 *
 */
struct _model {
	/* current hypothesis */
	int K;			/*!< number of subpopulations */
	int null_K, alt_K;	/*!< H0 and Ha K */
	/* parameters */
	int no_parameters;	/*!< number of parameters */
#ifndef OLDWAY
	double ****vpklm;
	double **vetak;
	double ***vetaik;
	int pindex;		/*!< index of previous iterate */
	int findex;		/*!< index of current from x */
	int tindex;		/*!< index of current to F(x) */
#else
	double ***pKLM;		/*!< allele probabilities */
	double *etak;		/*!< mixing proportions */
	double **etaik;		/*!< individual-specific mixing proportions */
#endif
	double ***mle_pKLM;
	double *mle_etak;	/*!< max. likelihood estimates of parameters */
	double **mle_etaik;
	double ****diklm;	/*!< posterior prob. allele from subpopn k */
	double **vik;		/*!< posterior prob. indiv. from subpopn k */

	/* acceleration */
#ifndef OLDWAY
	double ***u_etaik;	/*!< U matrices: */
	double ****u_pklm;
	double **u_etak;
	double ***v_etaik;	/*!< V matrices: */
	double ****v_pklm;
	double **v_etak;
	int delta_index;	/*!< index to place next deltas */
	int uindex;		/*!< index with u (penultimate) deltas */
	int vindex;		/*!< index with v (last) deltas */
#else
	double *init_etak;	/*!< saves all initialized etak values */
	double **init_etaik;	/*!< init_etaik saves all initialized etaik values */
	double ***init_pKLM;	/*!< saves all initialized pKLM values */
	double *iter1_etak;	/*!< saves all etak values from first iteration of EM */
	double *iter2_etak;	/*!< saves all etak values from second iteration of EM */
	double **iter1_etaik;	/*!< saves all iteration 1 etaik values */
	double **iter2_etaik;	/*!< saves all iteration 2 etaik values */
	double ***iter1_pKLM;	/*!< saves all iteration 1 pKLM values */
	double ***iter2_pKLM;	/*!< saves all iteration 2 pKLM values */
	double *V;		/*!< V vector = X_2-X_1 */
	double *U;		/*!< U vector = X_1-X_0*/
#endif

	/* QN acceleration */
	double *A;
	double *Ainv;
	double *cutu;
#ifdef LAPACK
	double *work;
	int *ipiv;
#endif

	/* state for current iteration */
	double logL;		/*!< current max. log likelihood */
	int current_i;		/*!< current_i is the current i to be updated in EM if block relaxation fails */
	int current_l;		/*!< current_l is the current l to be updated in EM if block relaxation fails */
	int current_k;		/*!< current_k is the current k to be updated in EM if block relaxation fails */
	int etaupdate;		/*!< indexes whether am currently updating etas or pklas 0=>pklas, 1=> etas */
	int accel_abort;	/*!< since we don't adjust step size: indicates step too big */
	int accel_step;		/*!< indicate if last step was accelerated */

	/* state for current initalization */
	int converged;		/*!< Indicator variable to flag if likelihood has converged */
	int stopped;		/*!< Indicator variable to flag if iterations must stop */
	int *count_K;		/*!< no. in each cluster in final assignment */
	int n_iter;		/*!< current number of iterations */

	/* bootstrap */
	double ts_obs;		/*!< observed test statistic */
	double ts_bs;		/*!< current bootstrap test statistic */
	double pvalue;		/*!< p-value for H0 via parametric bootstrap */

	/* maintain state across multiple initializations */
	int ever_converged;	/*!< has it converged on any initialization */
	double max_logL_H0;	/*!< max. log likelihood under H0 */
	double max_logL;	/*!< max. log likelihood (under HA) */
	double first_max_logL;	/*!< max. log likelihood (under HA) when first met */
	double aic;		/*!< AIC for max_logL solution */
	double bic;		/*!< BIC for max_logL solution */
	double arand;		/*!< adjusted RAND index for max_logL solution */
	int n_init;		/*!< number of initializations */
	int n_total_iter;	/*!< total number of iterations (over inits) */
	int n_maxll_init;	/*!< initialization producing max ll */
	int n_maxll_times;	/*!< no. initializations reaching max ll */
	int n_targetll_times;	/*!< no. inializations reaching target ll */
	int n_targetll_init;	/*!< first initialization reaching target ll */
	int time_stop;		/*!< estimation stopped b/c of time */
	int iter_stop;		/*!< estimation stopped b/c too many iterations */
	int n_max_iter;		/*!< maximum number of iterations */
	clock_t start;		/*!< time first initialization started */
	double seconds_run;	/*!< number of seconds run */

	/* state across multiple models */
	int aic_K;		/*!< estimated K by AIC */
	int bic_K;		/*!< estimated K by BIC */
}; /* model */

/**
 * Information about sampled individuals.  Stores the name and sampling locale
 * of each observed individual.
 */
struct _indiv {
	char *name;		/*!< name of individual */
	int locale;		/*!< locale of individual */
};

/* EM algorithm */
int initialize_model(options *opt, data *dat, model *mod);
void em(options *opt, data *dat, model *mod);
int em_step(options *opt, data *dat, model *mod);
double em_e_step(options *opt, data *dat, model *mod);
int em_2_steps(model *mod, data *dat, options *opt);
int stop(options *opt, model *mod, double loglik);
int converged(options *opt, model *mod, double loglik);
double log_likelihood(options *opt, data *dat, model *mod, int);
int loglikelihood_decrease(options *opt, double pll, double ll);
void m_step_admixture(options *opt, data *dat, model *mod);
double aic(model *mod);
double bic(data *dat, model *mod);

/* acceleration methods */
int accelerated_em_step(options *opt, data *dat, model *mod);
int simplex_project_pklm(model *mod, data *dat, options *opt, int, int);
int simplex_project_eta(model *mod, options *opt, int);

/* Block Relaxation by Zhou's Algorithm */
double block_update_params(options *opt, data *dat, model *mod);
void initialize_etaiks(model *mod);
void initialize_pklas(model *mod, data *dat);
void initialize_etaks(model *mod);
int set_params_em(model *mod, data *dat, options *opt);

/* bootstrap */
int parametric_bootstrap(options *opt, data *dat, model *mod);
void cleanup_parametric_bootstrap(data *dat);

/* input */
int read_file(options *opt, data *dat);
int read_admixture_qfile(data *dat, model *smod, char const *pfile);
int read_admixture_pfile(data *dat, model *smod, char const *qfile);
int make_ila(data *dat);
int read_pfile(options *opt, data *dat, model *mod);
int read_qfile(options *opt, data *dat, model *mod);
int read_afile(options *opt, data *dat);

/* output */
int write_file(char *outfile, options *opt, data *dat, model *mod);
int write_file_detail(options *opt, data *dat, model *mod);
int write_data(options *opt, data *dat, char const *outfile, int useILM);
void partition_admixture(data *dat, model *mod);
void partition_mixture(data *dat, model *mod);
int popq_admix(options *opt, data *dat, model *mod);
int popq_mix(options *opt, data *dat, model *mod);
int indivq_admix(options *opt, data *dat, model *mod);
int indivq_mix(options *opt, data *dat, model *mod);
void print_param(options *opt, data *dat, model *mod, int);
int print_param_diff(options *opt, data *dat, model *mod, int, int);
void print_udiff(options *opt, data *dat, model *mod, int);
void print_vdiff(options *opt, data *dat, model *mod, int);
void print_number_unique_alleles(data *dat);
void print_unique_alleles(data *dat);
void print_allele_counts(data *dat);
void print_sufficient_statistics(data *dat);

#endif
