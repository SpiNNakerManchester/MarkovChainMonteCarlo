#include <stdint.h>

typedef struct mcmc_params* mcmc_params_pointer_t;
typedef struct mcmc_state* mcmc_state_pointer_t;
typedef double (*t_deviate_function_t)();

// macro to define square( x ) = x^2
#define SQR( x ) (( x ) * ( x ))

uint32_t mcmc_model_get_params_n_bytes();

uint32_t mcmc_model_get_state_n_bytes();

double mcmc_model_likelihood(
    double x, mcmc_params_pointer_t params, mcmc_state_pointer_t state);

double mcmc_model_prior_prob(
    mcmc_params_pointer_t params, mcmc_state_pointer_t state);

void mcmc_model_transition_jump(
    t_deviate_function_t t_deviate, mcmc_params_pointer_t params,
    mcmc_state_pointer_t state, mcmc_state_pointer_t new_state);

static inline double mcmc_model_next_transition_jump(
        double value, t_deviate_function_t t_deviate, double jump_scale) {
    return value + (t_deviate() * jump_scale);
}
