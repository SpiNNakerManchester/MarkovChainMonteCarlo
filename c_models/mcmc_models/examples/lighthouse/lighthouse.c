#include "../../mcmc_model.h"
#include "lighthouse.h"

// definition of Pi for use in likelihood
double pi = 3.141592653589793;

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
double mcmc_model_likelihood(
        double x, mcmc_params_pointer_t params, mcmc_state_pointer_t state) {
    return state->beta / (pi * ( SQR( state->beta ) + SQR(x - state->alpha)));
}

/*
 prior probability for the two parameters alpha and beta = P( ALPHA, BETA )

 impossible that they are outside ranges, otherwise uniform
 */
double mcmc_model_prior_prob(
        mcmc_params_pointer_t params, mcmc_state_pointer_t state) {
    if (state->alpha < params->alpha_min || state->alpha > params->alpha_max ||
            state->beta < params->beta_min || state->beta > params->beta_max) {
        return 0.0;
    }
    return 1.0;
}


/*
 generate jumps from MH transition distribution which is uncorrelated bivariate
 t with degrees_freedom degrees of freedom
 */
void mcmc_model_transition_jump(
        t_deviate_function_t t_deviate, mcmc_params_pointer_t params,
        mcmc_state_pointer_t state, mcmc_state_pointer_t new_state) {
    new_state->alpha = mcmc_model_next_transition_jump(
        state->alpha, t_deviate, params->alpha_jump_scale);
    new_state->beta = mcmc_model_next_transition_jump(
        state->beta, t_deviate, params->beta_jump_scale);
}
