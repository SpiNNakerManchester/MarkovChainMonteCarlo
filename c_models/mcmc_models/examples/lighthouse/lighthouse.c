#include "../../mcmc_model.h"
#include "lighthouse.h"
//#include <debug.h>

// definition of Pi for use in likelihood
//#if CALC_TYPE == double
//CALC_TYPE PI = 3.141592653589793;
//#elif CALC_TYPE == float
//CALC_TYPE pi = 3.1415927
//#elif CALC_TYPE == accum
//CALC_TYPE pi = 3.141592k
//#endif

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
        CALC_TYPE x, mcmc_params_pointer_t params, mcmc_state_pointer_t state) {
//    log_info("mcmc_model_likelihood");
//	return state->beta / (PI * ( SQR( state->beta ) + SQR(x - state->alpha)));
	CALC_TYPE beta = state->beta;
	CALC_TYPE value = PI*(SQR(state->beta)+SQR(x-state->alpha));
//	CALC_TYPE logbeta = log_test(beta);
//	CALC_TYPE logvalue = log_test(value);
//	return logbeta - logvalue;
#if TYPE_SELECT == 2
	return LN(beta/value);
#else
	return LN(beta) - LN(value);
#endif
}

/*
 prior probability for the two parameters alpha and beta = P( ALPHA, BETA )

 impossible that they are outside ranges, otherwise uniform
 */
CALC_TYPE mcmc_model_prior_prob(
        mcmc_params_pointer_t params, mcmc_state_pointer_t state) {
    if (state->alpha < params->alpha_min || state->alpha > params->alpha_max ||
            state->beta < params->beta_min || state->beta > params->beta_max) {
        return ZERO;
    }
    return ONE;
}


/*
 generate jumps from MH transition distribution which is uncorrelated bivariate
 t with degrees_freedom degrees of freedom
 */
void mcmc_model_transition_jump(
        mcmc_params_pointer_t params, mcmc_state_pointer_t state,
        mcmc_state_pointer_t new_state) {
    new_state->alpha = state->alpha + (t_deviate() * params->alpha_jump_scale);
    new_state->beta = state->beta + (t_deviate() * params->beta_jump_scale);
}
