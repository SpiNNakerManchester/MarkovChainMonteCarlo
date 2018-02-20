#include "../../mcmc_model.h"

#define PPOLYORDER 9
#define QPOLYORDER 9

struct mcmc_params {

    // scaling of transition distribution for MH jumps
    CALC_TYPE p_jump_scale[PPOLYORDER];
    CALC_TYPE q_jump_scale[QPOLYORDER];
    CALC_TYPE mu_jump_scale;
    CALC_TYPE sigma_jump_scale;

};

struct mcmc_state {
    // parameters has size order_p (p poly) + order_q (q poly) + 2 (mu and sigma)
    CALC_TYPE parameters[PPOLYORDER+QPOLYORDER+2];
};
