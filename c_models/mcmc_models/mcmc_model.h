#include <stdint.h>

#define TYPE_SELECT 2 // 0 - double, 1 - float, 2 - accum

// Concatenate to add suffix to type where required
#define CONCAT_HELPER(a,b) a ## b
#define CONCAT(a,b) CONCAT_HELPER(a,b)

#if TYPE_SELECT == 0

#define CALC_TYPE double
#define SUFFIX 00000000000

#include <math.h>

#define LN( x ) log(x)
#define EXP( x ) exp(x)

#define POW( x, p ) pow((x), (p))
#define ABS( x ) fabs(x)

#define PI 3.141592653589793

#elif TYPE_SELECT == 1

#define CALC_TYPE float
#define SUFFIX f

#include <math.h>

#define LN( x ) logf(x)
#define EXP( x ) expf(x)

#define POW( x, p ) powf((x), (p))
#define ABS( x ) fabs(x)

#define PI 3.141593f

#elif TYPE_SELECT == 2

#define CALC_TYPE accum
#define SUFFIX k

#include <stdfix.h>
#include <stdfix-exp.h>
#include <log.h>
#include <math.h>

#define LN( x ) logk(x)
#define EXP( x ) expk(x)

#define POW( x, p ) pow((x), (p))
#define ABS( x ) absk(x)

#define PI 3.141593k
#define DEVIATE_TOL 0.002k

#endif

// Define constants by adding defined suffix value at end
#define ONE CONCAT(1.000000, SUFFIX)
#define HALF CONCAT(0.500000, SUFFIX)
#define ZERO CONCAT(0.000000, SUFFIX)
#define TWO CONCAT(2.000000, SUFFIX)
#define THREE CONCAT(3.000000, SUFFIX)
#define FOUR CONCAT(4.000000, SUFFIX)

typedef struct mcmc_params* mcmc_params_pointer_t;
typedef struct mcmc_state* mcmc_state_pointer_t;

extern CALC_TYPE t_deviate();

// macro to define square( x ) = x^2
#define SQR( x ) (( x ) * ( x ))

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

