#include <stdint.h>

#define TYPE_SELECT 2 // 0 - double, 1 - float, 2 - accum

#if TYPE_SELECT == 0

#define CALC_TYPE double

#include <math.h>

#define LN( x ) log(x)
#define EXP( x ) exp(x)

#define POW( x, p ) pow((x), (p))
#define ABS( x ) fabs(x)

#define ONE 1.00000000000000000
#define HALF 0.50000000000000000
#define ZERO 0.00000000000000000
#define TWO 2.00000000000000000
#define THREE 3.00000000000000000
#define FOUR 4.00000000000000000
#define TWOFIVESIX 256.00000000000000000
#define PI 3.141592653589793

#elif TYPE_SELECT == 1

#define CALC_TYPE float

#include <math.h>

#define LN( x ) logf(x)
#define EXP( x ) expf(x)

#define POW( x, p ) powf((x), (p))
#define ABS( x ) fabs(x)

#define ONE 1.0000000f
#define HALF 0.5000000f
#define ZERO 0.0000000f
#define TWO 2.0000000f
#define THREE 3.0000000f
#define FOUR 4.0000000f
#define TWOFIVESIX 256.0000000f
#define PI 3.141593f

#elif TYPE_SELECT == 2

#define CALC_TYPE accum

#include <stdfix.h>
#include <stdfix-exp.h>
#include <log.h>

#define LN( x ) logk(x)
#define EXP( x ) expk(x)

#define POW( x, p ) pow((x), (p))
#define ABS( x ) absk(x)

#define ONE 1.000000k
#define HALF 0.500000k
#define ZERO 0.000000k
#define TWO 2.000000k
#define THREE 3.000000k
#define FOUR 4.000000k
#define TWOFIVESIX 256.000000k
#define PI 3.141593k

#endif


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

//CALC_TYPE log_test(CALC_TYPE value)
//{
//	return LN(value);
//}

