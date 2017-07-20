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

#define TYPE_SELECT 2 // 0 - double, 1 - float, 2 - accum (as set in mcmc_model.h)

#include <limits.h>
//#include <math.h>
#include <spin1_api.h>
#include <debug.h>
#include <data_specification.h>
#include <simulation.h>
#include <recording.h>
#include "mcmc_model.h"
#include "examples/lighthouse/lighthouse.h"

// Define spin1_wfi
extern void spin1_wfi();

// The type of the seed
typedef uint32_t uniform_seed[5];

enum regions {
    RECORDING,
    PARAMETERS,
    MODEL_PARAMETERS,
    MODEL_STATE
};

struct parameters {

    // no of MCMC transitions to reach apparent equilibrium before generating
    // inference samples
    uint32_t burn_in;

    // subsequent MCMC samples are correlated, so thin the chain to avoid this
    uint32_t thinning;

    // number of posterior samples required
    uint32_t n_samples;

    // number of data points
    uint32_t n_data_points;

    // data receive window size - 0 if not a receiver
    uint32_t data_window_size;

    // Sequence mask for data reception - 0 if not a receiver
    uint32_t sequence_mask;

    // Acknowledge key for data reception - 0 if not a receiver
    uint32_t acknowledge_key;

    // Tag which is allocated for the data
    uint32_t data_tag;

    // Timer for data acknowledgement - 0 if not receiver
    uint32_t timer;

    // The random seed
    uniform_seed seed;

    // The number of degrees of freedom to jump around with
    CALC_TYPE degrees_of_freedom;
};

// setup variables for uniform PRNG
#if TYPE_SELECT == 0
CALC_TYPE uint_max_scale = 1.0 / UINT_MAX;
#elif TYPE_SELECT == 1
CALC_TYPE uint_max_scale = 1.0f / UINT_MAX;
#endif
// Don't need max_scale for fixed-point

// The general parameters
struct parameters parameters;

// The model-specific parameters
mcmc_params_pointer_t params;

// The model-specific state
mcmc_state_pointer_t state;

// Pointer to receive the data with
uint32_t *data_receive_ptr;

// The data to process
CALC_TYPE *data;

// The next sequence number expected
uint32_t next_sequence = 0;

// The last sequence number seen
uint32_t last_sequence = 0xFFFFFFFF;

// Do the DMA likelihood or not
uint dma_likelihood = 0;

// The number of points in the buffers
#define N_BUFFER_POINTS 128

// The size of the buffers in bytes
#define DMA_BUFFER_SIZE (N_BUFFER_POINTS << 3)

// The local data buffers
CALC_TYPE *dma_buffers[2];

// The current data buffer being read
uint32_t dma_read_buffer = 0;

// Flag to indicate when likelihood data transfer has been done
uint likelihood_done = 0;

//#ifdef USE_FP
// No print_value function required for fixed-point
//#else
#if TYPE_SELECT != 2
struct double_uint {
    uint first_word;
    uint second_word;
};

union double_to_ints {
    CALC_TYPE double_value;
    struct double_uint int_values;
};

void print_value(CALC_TYPE d_value, char *buffer) {
    union double_to_ints converter;
    converter.double_value = d_value;
    io_printf(
        buffer, "0x%08x%08x",
        converter.int_values.second_word, converter.int_values.first_word);
}
#endif
//#endif

// returns a high-quality Uniform[0,1] random variate -
// Marsaglia KISS32 algorithm
CALC_TYPE uniform(uniform_seed seed) {

    int t;

    seed[1] ^= (seed[1] << 5);
    seed[1] ^= (seed[1] >> 7);
    seed[1] ^= (seed[1] << 22);

    t = seed[2] + seed[3] + seed[4];
    seed[2] = seed[3];
    seed[4] = t < 0;
    seed[3] = t & 2147483647;
    seed[0] += 1411392427;

//    log_info("uniform, 0.000015k = %.10k", 0.000015k);

    uint int_value = seed[0] + seed[1] + seed[3];

#if TYPE_SELECT == 2
    return (CALC_TYPE) ulrbits(int_value);
#else
    return (CALC_TYPE) int_value * uint_max_scale; // float?
#endif
}

// Returns a standard t-distributed deviate - from Ripley
CALC_TYPE t_deviate() {
//	char buffer[1024];
    CALC_TYPE x;
    CALC_TYPE v;
    CALC_TYPE df = parameters.degrees_of_freedom;
#if TYPE_SELECT == 2
    CALC_TYPE rhs;
#endif

    do {
        CALC_TYPE u = uniform(parameters.seed);
        CALC_TYPE u1 = uniform(parameters.seed);

//        log_info("t_deviate, do 0.25*2048 = %k", 0.25k*2048.0k);
        if (u < HALF) {
#if TYPE_SELECT == 2
        	if (ABS(FOUR * u - ONE) < 0.002k) {
        		log_info("t_deviate(), abs(FOUR*u-ONE) = %k, u = %k",
        				ABS(FOUR * u - ONE), u);
        		x = 1000.0k; // or some large number...
        	} else {
#endif
        	x = ONE / (FOUR * u - ONE);
#if TYPE_SELECT == 2
        	}
#endif
        	v = (ONE / SQR(x)) * u1;
        } else {
        	x = FOUR * u - THREE;
        	v = u1;
        }
//        print_value(x, buffer);
//        log_info("t_deviate(), x = %s", buffer);

        if (v < (ONE - HALF * ABS(x))) {
//            print_value(x, buffer);
//            log_info("t_deviate(), x = %s", buffer);
            return x;
        }

#if TYPE_SELECT == 2
        if (df == THREE) {
        	rhs = ONE / SQR( ONE + SQR( x ) / THREE );
        } /*else {
        	rhs = POW((ONE + SQR( x ) / df), -(df + ONE) / TWO);
        }*/
    } while (v >= rhs);
#else
    } while (v >= POW((ONE + SQR( x ) / df), -(df + ONE) / TWO));
#endif

//    print_value(x, buffer);
//    log_info("t_deviate(), x = %s", buffer);
    return x;
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
bool MH_MCMC_keep_new_point(CALC_TYPE old_pt_posterior_prob,
        CALC_TYPE new_pt_posterior_prob, uniform_seed seed) {
	// Now we are taking log-likelihood, new_pt > old_pt when new_pt is zero
	// - to follow the same earlier logic with likelihood we need to call
	// the uniform(seed) function but always return false
	if (new_pt_posterior_prob == ZERO) {
		CALC_TYPE dummy = uniform(seed);
		return false;
	}

	// If old_pt is zero, then return true!
	if (old_pt_posterior_prob == ZERO)
		return true;

	// Now do actual tests if required
	// we have taken logs, but if
	if (new_pt_posterior_prob > old_pt_posterior_prob)
        return true;
	// log-domain, so can turn division into subtraction and use EXP
    else if (EXP(new_pt_posterior_prob-old_pt_posterior_prob) > uniform(seed))
//    else if ((new_pt_posterior_prob/old_pt_posterior_prob) > uniform(seed))
        return true;
    else
        return false;
}

void do_transfer(CALC_TYPE *dataptr, uint bytes) {
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
CALC_TYPE full_data_set_likelihood(mcmc_state_pointer_t state_to_use) {
//	char buffer[1024];
//    CALC_TYPE l = ONE;
    CALC_TYPE l = ZERO;
    if (!dma_likelihood) {
        // distribute these data points across cores?
        for (unsigned int i = 0; i < parameters.n_data_points; i++) {
//            l *= mcmc_model_likelihood(data[i], params, state_to_use);
            l += mcmc_model_likelihood(data[i], params, state_to_use);
        }
        return l;
    }

    // Keep a count of the points to be processed
    uint points_to_process = parameters.n_data_points;
    uint bytes_to_get = parameters.n_data_points * 8;
    uint bytes = DMA_BUFFER_SIZE;
    CALC_TYPE *dataptr = data;

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
//            l *= mcmc_model_likelihood(
        	l += mcmc_model_likelihood(
                    dma_buffers[dma_process_buffer][i], params, state_to_use);
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
//	char buffer[1024];
    use(unused0);
    use(unused1);

//    log_info("UINT_MAX is %u", UINT_MAX);

    // Create a new state pointer
    uint32_t state_n_bytes = mcmc_model_get_state_n_bytes();
    mcmc_state_pointer_t new_state = (mcmc_state_pointer_t) spin1_malloc(
        state_n_bytes);
    if (new_state == NULL) {
        log_error("Could not allocate space for new state");
        rt_error(RTE_SWERR);
    }

//    CALC_TYPE log_likelihood;
//    CALC_TYPE likelihood;
    CALC_TYPE current_posterior;
    CALC_TYPE new_posterior;
    unsigned int sample_count = 0;
    unsigned int accepted = 0;
    unsigned int likelihood_calls = 0;

    // Try to copy data in to DTCM
    CALC_TYPE *data_ptr = (CALC_TYPE *) sark_tag_ptr(
        parameters.data_tag, sark_app_id());
    data = (CALC_TYPE *) spin1_malloc(
    		parameters.n_data_points*sizeof(CALC_TYPE));
    if (data != NULL) {
        spin1_memcpy(data, data_ptr,
        		parameters.n_data_points*sizeof(CALC_TYPE));
        dma_likelihood = 0;
    } else {
        uint space = sark_heap_max(sark.heap, 0);
        log_warning(
            "Could not allocate data of size %d to DTCM (%d bytes free)"
            "- using DMAs", parameters.n_data_points*sizeof(CALC_TYPE), space);
        data = data_ptr;

        // Allocate the buffers
        dma_buffers[0] = (CALC_TYPE *) spin1_malloc(DMA_BUFFER_SIZE);
        dma_buffers[1] = (CALC_TYPE *) spin1_malloc(DMA_BUFFER_SIZE);

        // Set up the callback handler
        spin1_callback_on(DMA_TRANSFER_DONE, dma_callback, 0);
        dma_likelihood = 1;

        // Do the first transfer
        do_transfer(data, DMA_BUFFER_SIZE);
    }

    // exponential to get likelihood back from log likelihood
//    log_likelihood = full_data_set_likelihood(state);
//    likelihood = expk(log_likelihood);
//
//    log_info("log_likelihood = %k, likelihood = %k",
//    		log_likelihood, likelihood);

    // first posterior calculated at initial state values
    current_posterior = full_data_set_likelihood(state) *
    		mcmc_model_prior_prob(params, state);
//            full_data_set_likelihood(state) * mcmc_model_prior_prob(params, state);


    // update likelihood function counter for diagnostics
    likelihood_calls++;

    uint samples_to_go = parameters.thinning;
    bool burn_in = true;

    do {

        // make a jump around parameter space using bivariate t distribution
        // with 3 degrees of freedom
        mcmc_model_transition_jump(params, state, new_state);

//        log_info("state->beta = %k, new_state->beta = %k",
//        		state->beta, new_state->beta);
//        print_value(new_state->alpha, buffer);
//        log_info("new_state->alpha = %s", buffer);
//        print_value(new_state->beta, buffer);
//        log_info("new_state->beta = %s", buffer);

        // update likelihood function counter for diagnostics
        likelihood_calls++;

        // get likelihood from log likelihood
//        likelihood = expk(full_data_set_likelihood(new_state));

        // calculate joint probability at this point
        new_posterior = full_data_set_likelihood(new_state) *
				mcmc_model_prior_prob(params, new_state);

//        if (!burn_in) {
//#if TYPE_SELECT == 2
//        	log_info("new_posterior = %k, current_posterior = %k",
//        			new_posterior, current_posterior);
//#endif
//        }

//        log_info("seed is %d %d %d %d %d", parameters.seed[0], parameters.seed[1],
//        		parameters.seed[2], parameters.seed[3], parameters.seed[4]);
        // if accepted, update current state, otherwise leave it as is
        if (MH_MCMC_keep_new_point(
                current_posterior, new_posterior, parameters.seed)) {

            current_posterior = new_posterior;
            spin1_memcpy(state, new_state, state_n_bytes);

//            print_value(EXP(2.00000), buffer);
//            log_info("point accepted %d", accepted);

            // update acceptance count
            accepted++;
        };

        if (burn_in) {
            if (likelihood_calls == parameters.burn_in) {
                log_info(
                    "Burn-in accepted %d of %d", accepted, likelihood_calls);

                // reset diagnostic statistics
                accepted = 0;
                likelihood_calls = 0;
                burn_in = false;
            }
        } else {

//        	log_info("starting sampling, sample_count = %d, to_go = %d",
//        			sample_count, samples_to_go);
            // output every THINNING samples
            if (samples_to_go == 0) {
//                log_info("recording: state->beta = %k, state->alpha = %k",
//                		state->beta, state->alpha);
                recording_record(0, state, state_n_bytes);
                sample_count++;
                samples_to_go = parameters.thinning;
            }
            samples_to_go--;
        }

    } while (sample_count < parameters.n_samples);

    recording_finalise();

    log_info("Sampling accepted %d of %d", accepted, likelihood_calls);
    spin1_exit(0);
}

void multicast_callback(uint key, uint payload) {
    uint sequence = key & parameters.sequence_mask;
    if (sequence == next_sequence) {
        data_receive_ptr[0] = payload;
        data_receive_ptr++;
        last_sequence = sequence;
        next_sequence = (sequence + 1) & parameters.sequence_mask;
    }
}

void timer_callback(uint time, uint unused) {
    spin1_delay_us(parameters.timer >> 1);
    use(time);
    use(unused);
    spin1_send_mc_packet(
        parameters.acknowledge_key, last_sequence, WITH_PAYLOAD);
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
#if TYPE_SELECT != 2
    char buffer[1024];
#endif

    // Read the data specification header
    address_t data_address = data_specification_get_data_address();
    if (!data_specification_read_header(data_address)) {
        rt_error(RTE_SWERR);
    }

    // Setup recording
    address_t recording_address = data_specification_get_region(
        RECORDING, data_address);
    uint32_t recording_flags = 0;
    if (!recording_initialize(recording_address, &recording_flags)) {
        rt_error(RTE_SWERR);
    }

    // Read the parameters
    address_t parameters_address = data_specification_get_region(
        PARAMETERS, data_address);
    struct parameters *sdram_params = (struct parameters *) parameters_address;
    spin1_memcpy(&parameters, sdram_params, sizeof(struct parameters));

    // Set the last sequence to the one "before" 0
    last_sequence = parameters.sequence_mask;

    // Read the model parameters
    address_t model_params_address = data_specification_get_region(
        MODEL_PARAMETERS, data_address);
    uint32_t params_n_bytes = mcmc_model_get_params_n_bytes();
    params = (mcmc_params_pointer_t) spin1_malloc(params_n_bytes);
    if (params == NULL) {
        log_error("Could not allocate model parameters");
        rt_error(RTE_SWERR);
    }
    spin1_memcpy(params, model_params_address, params_n_bytes);

    // Read the model state
    address_t model_state_address = data_specification_get_region(
        MODEL_STATE, data_address);
    uint32_t state_n_bytes = mcmc_model_get_state_n_bytes();
    state = (mcmc_state_pointer_t) spin1_malloc(state_n_bytes);
    if (state == NULL) {
        log_error("Could not allocate model state");
        rt_error(RTE_SWERR);
    }
    spin1_memcpy(state, model_state_address, state_n_bytes);

    log_info("Burn in = %d", parameters.burn_in);
//    log_info("Thinning = %d", parameters.thinning);
//    log_info("N Samples = %d", parameters.n_samples);
//    log_info("N Data Points = %d", parameters.n_data_points);
//    log_info("Data Window Size = %d", parameters.data_window_size);
//    log_info("Sequence mask = 0x%08x", parameters.sequence_mask);
//    log_info("Acknowledge key = 0x%08x", parameters.acknowledge_key);
//    log_info("Data tag = %d", parameters.data_tag);
//    log_info("Timer = %d", parameters.timer);
#if TYPE_SELECT == 2
    log_info("Degrees of freedom = %k", parameters.degrees_of_freedom);
#else
    print_value(parameters.degrees_of_freedom, buffer);
    log_info("Degrees of freedom = %s", buffer);
#endif

    // Allocate the data receive space if this is the nominated receiver
    if (parameters.data_window_size > 0) {
        data_receive_ptr = (uint32_t *) sark_xalloc(
            sv->sdram_heap, parameters.n_data_points * sizeof(CALC_TYPE),
            parameters.data_tag, ALLOC_LOCK);
        if (data_receive_ptr == NULL) {
            log_error(
                "Could not allocate data array of size %d in SDRAM",
                parameters.n_data_points * sizeof(CALC_TYPE));
            rt_error(RTE_SWERR);
        }

        // Register for the data
        spin1_callback_on(MCPL_PACKET_RECEIVED, multicast_callback, -1);

        // Set up the timer for acknowledgement
        spin1_set_timer_tick(parameters.timer);
        spin1_callback_on(TIMER_TICK, timer_callback, 0);
    } else {
        spin1_callback_on(MCPL_PACKET_RECEIVED, empty_multicast_callback, -1);
    }

    // Register for the start message, even if not the receiver
    spin1_callback_on(MC_PACKET_RECEIVED, trigger_run, -1);

    spin1_start(SYNC_WAIT);
}

