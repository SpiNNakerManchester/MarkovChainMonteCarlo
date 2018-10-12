#include <mcmc_model.h>
#include "lighthouse.h"

#include <spin1_api.h>

#ifndef use
#define use(x) do {} while ((x)!=(x))
#endif

uint32_t mcmc_model_get_params_n_bytes() {
    return sizeof(struct mcmc_params);
}

uint32_t mcmc_model_get_state_n_bytes() {
    return sizeof(struct mcmc_state);
}

/*
 likelihood of a single data point x given the two parameters alpha and beta
 that describe the position of the lighthouse = P( X | ALPHA, BETA )

 - alpha defines distance along the shoreline from reference zero point

 - beta defines distance from the shoreline

 see accompanying documents and Sivia book for more detail
 */
CALC_TYPE mcmc_model_likelihood(
        CALC_TYPE *data, uint32_t n_pts, mcmc_params_pointer_t params,
		mcmc_state_pointer_t state) {
	use(params);
	CALC_TYPE beta = state->beta;
	// Need to do the loop in here now
	CALC_TYPE sum = ZERO;
	for (unsigned int i=0; i < n_pts; i++) {
		CALC_TYPE value = PI*(SQR(state->beta)+SQR(data[i]-state->alpha));
#if TYPE_SELECT == 2
//	return LN(beta/value);
		CALC_TYPE value2 = LN(beta) - LN(value);
	    sum += value2; // note this only works with version 5.0+
	                   // of arm-none-eabi-gcc; if you have an earlier
	                   // version install the new one or comment out
	                   // and use the line above instead
#else
	    sum += LN(beta) - LN(value);
#endif
	}
	return sum;
}

/*
 prior probability for the two parameters alpha and beta = P( ALPHA, BETA )

 impossible that they are outside ranges, otherwise uniform
 */
CALC_TYPE mcmc_model_prior_prob(
        mcmc_params_pointer_t params, mcmc_state_pointer_t state) {
	// we made the likelihood into log-likelihood, so this needs to do the same?
    if (state->alpha < params->alpha_min || state->alpha > params->alpha_max ||
            state->beta < params->beta_min || state->beta > params->beta_max) {
    	// log(0) = -inf
        return BIGNEG; // ZERO;
    }
    // log(1) = 0
    return ZERO; // ONE;
}


/*
 generate jumps from MH transition distribution which is uncorrelated bivariate
 t with degrees_freedom degrees of freedom
 */
void mcmc_model_transition_jump(
        mcmc_params_pointer_t params, mcmc_state_pointer_t state,
        mcmc_state_pointer_t new_state, uint32_t timestep) {
	//use(timestep);
    new_state->alpha = state->alpha + (t_deviate() * params->alpha_jump_scale);
    new_state->beta = state->beta + (t_deviate() * params->beta_jump_scale);
}

/*
 exit function from MCMC: single-core executable, so simply spin1_exit
 */
void mcmc_exit_function() {
	spin1_exit(0);
}

/*
 set up address and key function: in this single-core exec. case, empty function
 */
void mcmc_get_address_and_key() {

}
