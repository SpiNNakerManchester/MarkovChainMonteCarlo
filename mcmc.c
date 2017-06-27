/*

 simple_MCMC.c

 A 2D inference implementation using basic Metropolis-Hastings MCMC.

 This is known as the 'lighthouse problem', where light flashes are observed
 along a shoreline containing sensors (with infinite spatial precision!) and
 one wants to infer the position and distance of the source lighthouse from
 this set of data.

 The light flashes are uniformly distributed around the circular rotation of
 the light source.

 Various supporting documents which give more detail have been attached with
 this C code file.

 Non-trivial inference problem because of the distribution induced by the
 problem specification.

 Created by Michael Hopkins on 28/4/16.

 */

#include <limits.h>
#include <math.h>
#include <spin1_api.h>
#include <debug.h>
#include <data_specification.h>
#include <simulation.h>
#include <recording.h>

// The type of the seed
typedef uint32_t uniform_seed[5];

enum regions {
    RECORDING = 0,
    PARAMETERS,
    BUFFER_STATE_REGION,
    RECORDED_DATA
};

struct parameters {

    // Generic MCMC
    uint32_t burn_in;
    uint32_t thinning;
    uint32_t n_samples;
    uint32_t n_data_points;
    double degrees_of_freedom;
    uniform_seed seed;

    // SpiNNaker Data Loading
    uint32_t data_window_size;
    uint32_t sequence_mask;
    uint32_t acknowledge_key;
    uint32_t data_tag;
    uint32_t timer;

    // Problem specific
    double alpha_jump_scale;
    double beta_jump_scale;
    double alpha_min;
    double alpha_max;
    double beta_min;
    double beta_max;
};

// definition of Pi for use in likelihood
double pi = 3.141592653589793;

// 1
double ONE = 1.00000000000000;

// macro to define square( x ) = x^2
#define SQR( x ) (( x ) * ( x ))

// setup variables for uniform PRNG
double uint_max_scale = 1.0 / UINT_MAX;

// The parameters
uint32_t burn_in;
uint32_t thinning;
uint32_t n_samples;
uint32_t n_data_points;
uint32_t data_window_size;
uint32_t sequence_mask;
uint32_t acknowledge_key;
uint32_t data_tag;
double alpha_jump_scale;
double beta_jump_scale;
double alpha_min;
double alpha_max;
double beta_min;
double beta_max;
double degrees_of_freedom;
uniform_seed seed;

// Pointer to receive the data with
uint32_t *data_receive_ptr;

// The data to process
double *data;

// The next sequence number expected
uint32_t next_sequence = 0;

// The last sequence number seen
uint32_t last_sequence = 0xFFFFFFFF;

// The timer value
uint32_t timer;

// Do the DMA likelihood or not
uint dma_likelihood = 0;

// The number of points in the buffers
#define N_BUFFER_POINTS 128

// The size of the buffers in bytes
#define DMA_BUFFER_SIZE (N_BUFFER_POINTS << 3)

// The local data buffers
double *dma_buffers[2];

// The current data buffer being read
uint32_t dma_read_buffer = 0;

// Flag to indicate when likelihood data transfer has been done
uint likelihood_done = 0;

struct double_uint {
    uint first_word;
    uint second_word;
};

union double_to_ints {
    double double_value;
    struct double_uint int_values;
};

void print_value(double d_value, char *buffer) {
    union double_to_ints converter;
    converter.double_value = d_value;
    io_printf(
        buffer, "0x%08x%08x",
        converter.int_values.second_word, converter.int_values.first_word);
}

// returns a high-quality Uniform[0,1] random variate -
// Marsaglia KISS32 algorithm
double uniform(uniform_seed seed) {

    int t;

    seed[1] ^= (seed[1] << 5);
    seed[1] ^= (seed[1] >> 7);
    seed[1] ^= (seed[1] << 22);

    t = seed[2] + seed[3] + seed[4];
    seed[2] = seed[3];
    seed[4] = t < 0;
    seed[3] = t & 2147483647;
    seed[0] += 1411392427;

    return (double)
        ((unsigned int) seed[0] + seed[1] + seed[3]) * uint_max_scale;
}

// Returns a standard t-distributed deviate - from Ripley
double t_deviate(const double df, uniform_seed seed) {
    double u, u1, x, v;

    again: u = uniform(seed);
    u1 = uniform(seed);

    if (u < 0.5) {
        x = ONE / (4.0 * u - ONE);
        v = (ONE / SQR(x)) * u1;
    } else {
        x = 4.0 * u - 3.0;
        v = u1;
    }

    if (v < (ONE - 0.5 * fabs(x)))
        return x;

    if (v >= pow((ONE + SQR( x ) / df), -(df + ONE) / 2.0))
        goto again;

    return x;
}

/*
 likelihood of a single data point x given the two parameters alpha and beta
 that describe the position of the lighthouse = P( X | ALPHA, BETA )

 - alpha defines distance along the shoreline from reference zero point

 - beta defines distance from the shoreline

 see accompanying documents and Sivia book for more detail
 */
// Problem Specific
double likelihood(double x, double alpha, double beta) {
    return beta / (pi * ( SQR( beta ) + SQR(x - alpha)));
}

/*
 prior probability for the two parameters alpha and beta = P( ALPHA, BETA )

 impossible that they are outside ranges, otherwise uniform
 */
// Problem Specific
double prior_prob(
        double alpha, double beta, double alpha_min, double alpha_max,
        double beta_min, double beta_max) {
    if (alpha < alpha_min || alpha > alpha_max ||
            beta < beta_min || beta > beta_max) {
        return 0.0;
    } else {
        return 1.0;
    }
}

/*
 Metropolis-Hastings MCMC update sample test

 if posterior probability at new point > posterior probability at old point,
     keep new point
 else if ( posterior probability at new point /
           posterior probability at old point ) > ~Uniform[0, 1],
     keep new point
 else
     keep old point

 */
bool MH_MCMC_keep_new_point(double old_pt_posterior_prob,
        double new_pt_posterior_prob, uniform_seed seed) {
    if (new_pt_posterior_prob > old_pt_posterior_prob)
        return true;
    else if ((new_pt_posterior_prob / old_pt_posterior_prob) > uniform(seed))
        return true;
    else
        return false;
}


/*
 generate jumps from MH transition distribution which is uncorrelated bivariate
 t with degrees_freedom degrees of freedom
 */
// Problem Specific (t_deviate could be passed in?)
void transition_jump(double alpha, double beta, double *new_alpha,
        double *new_beta, double degrees_freedom, uniform_seed seed,
        double alpha_jump_scale, double beta_jump_scale) {
    *new_alpha = alpha + t_deviate(degrees_freedom, seed) * alpha_jump_scale;
    *new_beta = beta + t_deviate(degrees_freedom, seed) * beta_jump_scale;
}

void do_transfer(double *dataptr, uint bytes) {
    likelihood_done = 0;
    dma_read_buffer = (dma_read_buffer + 1) & 1;
    spin1_dma_transfer(
        0, dataptr, dma_buffers[dma_read_buffer], DMA_READ, bytes);
}

/*
 computes the likelihood for any given parameter combination over the whole
 data set

 **** this is one of the massive parallelisation opportunities i.e. distribute
 **** one likelihood calculation per core and all returned to a collating core
 **** for product (or sum if using log-likelihoods)

 */
// Not Problem Specific in general - only because it uses alpha and beta in
// this case
double full_data_set_likelihood(double alpha, double beta) {
    double l = ONE;
    if (!dma_likelihood) {

        // distribute these data points across cores?
        for (unsigned int i = 0; i < n_data_points; i++) {
            l *= likelihood(data[i], alpha, beta);
        }
        return l;
    }

    // Keep a count of the points to be processed
    uint points_to_process = n_data_points;
    uint bytes_to_get = n_data_points * 8;
    uint bytes = DMA_BUFFER_SIZE;
    double *dataptr = data;

    // Continue until all the points have been processed
    while (points_to_process > 0) {

        // Wait for the transfer to complete
        while (!likelihood_done) {
            spin1_wfi();
        }
        uint32_t dma_process_buffer = dma_read_buffer;

        // Work out how many points have been processed
        uint points = bytes >> 3;

        // Start the next transfer
        bytes_to_get -= bytes;
        bytes = bytes_to_get;
        if (bytes > DMA_BUFFER_SIZE) {
            bytes = DMA_BUFFER_SIZE;
        }
        if (bytes > 0) {
            dataptr = &(dataptr[points]);
            do_transfer(dataptr, bytes);
        } else {

            // If the data is all processed, to the first transfer for the next
            // calculation
            do_transfer(data, DMA_BUFFER_SIZE);
        }

        // Process the points in the buffer
        for (uint i = 0; i < points; i++) {
            l *= likelihood(dma_buffers[dma_process_buffer][i], alpha, beta);
        }

        points_to_process -= points;
    }
    return l;
}

void dma_callback(uint unused0, uint unused1) {
    use(unused0);
    use(unused1);
    likelihood_done = 1;
}

void run(uint unused0, uint unused1) {
    use(unused0);
    use(unused1);

    // Parameterise these with a struct
    double current_alpha = 0.0;
    double current_beta = 1.0;
    double new_alpha;
    double new_beta;

    // Generic
    double current_posterior;
    double new_posterior;
    unsigned int i;
    unsigned int sample_count = 0;
    unsigned int accepted = 0;
    unsigned int likelihood_calls = 0;

    // Try to copy data in to DTCM
    double *data_ptr = (double *) sark_tag_ptr(data_tag, sark_app_id());
    data = (double *) spin1_malloc(n_data_points * sizeof(double));
    if (data != NULL) {
        spin1_memcpy(data, data_ptr, n_data_points * sizeof(double));
        dma_likelihood = 0;
    } else {
        uint space = sark_heap_max(sark.heap, 0);
        log_warning(
            "Could not allocate data of size %d to DTCM (%d bytes free)"
            "- using DMAs", n_data_points * sizeof(double), space);
        data = data_ptr;

        // Allocate the buffers
        dma_buffers[0] = (double *) spin1_malloc(DMA_BUFFER_SIZE);
        dma_buffers[1] = (double *) spin1_malloc(DMA_BUFFER_SIZE);

        // Set up the callback handler
        spin1_callback_on(DMA_TRANSFER_DONE, dma_callback, 0);
        dma_likelihood = 1;

        // Do the first transfer
        do_transfer(data, DMA_BUFFER_SIZE);
    }

    // first posterior calculated at initialisation values of alpha and beta
    current_posterior = full_data_set_likelihood(current_alpha, current_beta) *
        prior_prob(
            current_alpha, current_beta, alpha_min, alpha_max,
            beta_min, beta_max);

    // update likelihood function counter for diagnostics
    likelihood_calls++;

    // burn-in phase
    for (i = 0; i < burn_in; i++) {

        // make a jump around parameter space using bivariate t distribution
        // with 3 degrees of freedom
        transition_jump(
            current_alpha, current_beta, &new_alpha, &new_beta,
            degrees_of_freedom, seed, alpha_jump_scale, beta_jump_scale);

        // update likelihood function counter for diagnostics
        likelihood_calls++;

        // calculate joint probability at this point
        new_posterior = full_data_set_likelihood(new_alpha, new_beta) *
            prior_prob(
                new_alpha, new_beta, alpha_min, alpha_max,
                beta_min, beta_max);

        // if accepted, update current state, otherwise leave it as is
        if (MH_MCMC_keep_new_point(current_posterior, new_posterior, seed)) {

            current_posterior = new_posterior;
            current_alpha = new_alpha;
            current_beta = new_beta;

            // update acceptance count
            accepted++;
        };

    }

    log_info("Burn-in accepted %d of %d", accepted, likelihood_calls);

    // reset diagnostic statistics
    accepted = likelihood_calls = 0;

    uint samples_to_go = thinning - 1;

    // now burnt in, start to collect inference sample output
    do {

        // **** this is the other massive parallelisation opportunities
        // i.e. multiple parallel chains ****

        // make a jump around parameter space
        transition_jump(
            current_alpha, current_beta, &new_alpha, &new_beta,
            degrees_of_freedom, seed, alpha_jump_scale, beta_jump_scale);

        likelihood_calls++; // update likelihood function counter for diagnostics

         // calculate joint probability at this point
        new_posterior = full_data_set_likelihood(new_alpha, new_beta) *
            prior_prob(
                new_alpha, new_beta, alpha_min, alpha_max,
                beta_min, beta_max);

        // if accepted, update current state, otherwise leave it as is
        if (MH_MCMC_keep_new_point(current_posterior, new_posterior, seed)) {
            current_posterior = new_posterior;
            current_alpha = new_alpha;
            current_beta = new_beta;

            // update acceptance count
            accepted++;
        };

        // output every THINNING samples
        if (samples_to_go == 0) {
            recording_record(0, &current_alpha, 8);
            recording_record(0, &current_beta, 8);
            sample_count++;
            samples_to_go = thinning;
        }
        samples_to_go--;

    // until you have enough posterior samples
    } while (sample_count < n_samples);

    recording_finalise();

    log_info("Sampling accepted %d of %d", accepted, likelihood_calls);
    spin1_exit(0);
}

void multicast_callback(uint key, uint payload) {
    uint sequence = key & sequence_mask;
    if (sequence == next_sequence) {
        data_receive_ptr[0] = payload;
        data_receive_ptr++;
        last_sequence = sequence;
        next_sequence = (sequence + 1) & sequence_mask;
    }
}

void timer_callback(uint time, uint unused) {
    spin1_delay_us(timer >> 1);
    use(time);
    use(unused);
    spin1_send_mc_packet(acknowledge_key, last_sequence, WITH_PAYLOAD);
}

void empty_multicast_callback(uint key, uint payload) {
    use(key);
    use(payload);
}

void trigger_run(uint unused0, uint unused1) {
    use(unused0);
    use(unused1);
    spin1_callback_off(TIMER_TICK);
    spin1_schedule_callback(run, 0, 0, 2);
}


/*
 main program which implements MCMC and outputs sample points from the
 posterior distribution using Bayes' theorem

 P( ALPHA, BETA | X ) ~= P( ALPHA, BETA ) * P( X | ALPHA, BETA )
 */
void c_main() {
    char buffer[1024];

    // Read the data specification header
    address_t data_address = data_specification_get_data_address();
    if (!data_specification_read_header(data_address)) {
        rt_error(RTE_SWERR);
    }

    // Setup recording
    address_t recording_address = data_specification_get_region(
        RECORDING, data_address);
    uint8_t region_ids[1] = {RECORDED_DATA};
    uint32_t recording_flags = 0;
    if (!recording_initialize(
            1, region_ids, recording_address, BUFFER_STATE_REGION,
            &recording_flags)) {
        rt_error(RTE_SWERR);
    }

    // Read the parameters
    address_t parameters_address = data_specification_get_region(
        PARAMETERS, data_address);
    struct parameters *params = (struct parameters *) parameters_address;

    // MCMC setup definitions
    // no of MCMC transitions to reach apparent equilibrium before generating
    // inference samples
    burn_in = params->burn_in;
    log_info("Burn in = %d", burn_in);

    // subsequent MCMC samples are correlated, so thin the chain to avoid this
    thinning = params->thinning;
    log_info("Thinning = %d", thinning);

    // number of posterior samples required
    n_samples = params->n_samples;
    log_info("N Samples = %d", n_samples);

    // number of data points
    n_data_points = params->n_data_points;
    log_info("N Data Points = %d", n_data_points);

    // data receive window size - 0 if not a receiver
    data_window_size = params->data_window_size;
    log_info("Data Window Size = %d", data_window_size);

    // Sequence mask for data reception - 0 if not a receiver
    sequence_mask = params->sequence_mask;
    last_sequence = sequence_mask;
    log_info("Sequence mask = 0x%08x", sequence_mask);

    // Acknowledge key for data reception - 0 if not a receiver
    acknowledge_key = params->acknowledge_key;
    log_info("Acknowledge key = 0x%08x", acknowledge_key);

    // Tag which is allocated for the data
    data_tag = params->data_tag;
    log_info("Data tag = %d", data_tag);

    // Timer for data acknowledgement - 0 if not receiver
    timer = params->timer;
    log_info("Timer = %d", timer);

    // scaling of t transition distribution for MH jumps in alpha direction
    alpha_jump_scale = params->alpha_jump_scale;
    print_value(alpha_jump_scale, buffer);
    log_info("Alpha jump scale = %s", buffer);

    // scaling of t transition distribution for MH jumps in beta direction
    beta_jump_scale = params->beta_jump_scale;
    print_value(beta_jump_scale, buffer);
    log_info("Beta jump scale = %s", buffer);

    // Alpha range
    alpha_min = params->alpha_min;
    print_value(alpha_min, buffer);
    log_info("Alpha min = %s", buffer);
    alpha_max = params->alpha_max;
    print_value(alpha_max, buffer);
    log_info("Alpha max = %s", buffer);

    // Beta range
    beta_min = params->beta_min;
    print_value(beta_min, buffer);
    log_info("Beta min = %s", buffer);
    beta_max = params->beta_max;
    print_value(beta_max, buffer);
    log_info("Beta max = %s", buffer);

    // The number of degrees of freedom to jump around with
    degrees_of_freedom = params->degrees_of_freedom;
    print_value(degrees_of_freedom, buffer);
    log_info("Degrees of freedom = %s", buffer);

    // The random seed
    sark_word_cpy(seed, params->seed, sizeof(uniform_seed));

    // Allocate the data receive space if this is the nominated receiver
    if (data_window_size > 0) {
        data_receive_ptr = (uint32_t *) sark_xalloc(
            sv->sdram_heap, n_data_points * sizeof(double),
            data_tag, ALLOC_LOCK);
        if (data_receive_ptr == NULL) {
            log_error(
                "Could not allocate data array of size %d in SDRAM",
                n_data_points * sizeof(double));
            rt_error(RTE_SWERR);
        }

        // Register for the data
        spin1_callback_on(MCPL_PACKET_RECEIVED, multicast_callback, -1);

        // Set up the timer for acknowledgement
        spin1_set_timer_tick(timer);
        spin1_callback_on(TIMER_TICK, timer_callback, 0);
    } else {
        spin1_callback_on(MCPL_PACKET_RECEIVED, empty_multicast_callback, -1);
    }

    // Register for the start message, even if not the receiver
    spin1_callback_on(MC_PACKET_RECEIVED, trigger_run, -1);

    spin1_start(SYNC_WAIT);
}

