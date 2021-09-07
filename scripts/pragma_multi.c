#include <omp.h>

if (opt->test_run){
	logL = 0.0;
	max_logL = -INFINITY;

	#pragma omp parallel for reduction(max:logL){
		for (i = 0; opt->target_revisit || opt->target_ll	/* targeting */
			|| opt->n_seconds || i < opt->n_init; i++) {	/* timing */
				
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

			/* record any better solution than previously seen */
			// needs to be run uniquely in multiprocessing
			if (mod->logL > mod->max_logL) {
				mod->max_logL = mod->logL;
			}
		}
	}
}