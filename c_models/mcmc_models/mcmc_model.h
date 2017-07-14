#include <stdint.h>
#include <stdfix.h>
//#include <log.h>

//#define CALC_TYPE double
#define CALC_TYPE accum

typedef struct mcmc_params* mcmc_params_pointer_t;
typedef struct mcmc_state* mcmc_state_pointer_t;

extern CALC_TYPE t_deviate();

// macro to define square( x ) = x^2
#define SQR( x ) (( x ) * ( x ))

// macro for log (from s1615)
//#define LOG( x ) logk(x)

// !\brief Get the number of bytes in the params struct (usually just sizeof)
uint32_t mcmc_model_get_params_n_bytes();

// !\brief Get the number of bytes in the state struct (usually just sizeof)
uint32_t mcmc_model_get_state_n_bytes();

// !\brief Given a value x, calculate the likelihood of the value
CALC_TYPE mcmc_model_likelihood(
    CALC_TYPE x, mcmc_params_pointer_t params, mcmc_state_pointer_t state);

// !\brief Get the prior probability of a given state
CALC_TYPE mcmc_model_prior_prob(
    mcmc_params_pointer_t params, mcmc_state_pointer_t state);

// !\brief Jump to a new state from the current state
void mcmc_model_transition_jump(
    mcmc_params_pointer_t params, mcmc_state_pointer_t state,
    mcmc_state_pointer_t new_state);
