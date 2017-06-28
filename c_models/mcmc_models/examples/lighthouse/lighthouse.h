struct mcmc_params {

    // scaling of t transition distribution for MH jumps
    double alpha_jump_scale;
    double beta_jump_scale;

    // Alpha range
    double alpha_min;
    double alpha_max;

    // Beta range
    double beta_min;
    double beta_max;
};

struct mcmc_state {
    double alpha;
    double beta;
};
