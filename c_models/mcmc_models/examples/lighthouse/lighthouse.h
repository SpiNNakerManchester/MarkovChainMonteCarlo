struct mcmc_params {

    // scaling of t transition distribution for MH jumps
    CALC_TYPE alpha_jump_scale;
    CALC_TYPE beta_jump_scale;

    // Alpha range
    CALC_TYPE alpha_min;
    CALC_TYPE alpha_max;

    // Beta range
    CALC_TYPE beta_min;
    CALC_TYPE beta_max;
};

struct mcmc_state {
    CALC_TYPE alpha;
    CALC_TYPE beta;
};
