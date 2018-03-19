#include "mcmc_cholesky.h"
#include "../mcmc_models/mcmc_model.h"
//#include "../mcmc_models/examples/arma/arma.h"

#include <spin1_api.h>
#include <stdint.h>
#include <stdbool.h>
#include <debug.h>
#include <data_specification.h>


// Value of the pointer to the location in SDRAM to get parameters
uint32_t *parameter_rec_ptr;  // will need something like this...

uint32_t *t_jump_ptr;

// Functions required to print floats as hex values (uncomment for debug)
//struct double_uint {
//    uint first_word;
//    uint second_word;
//};
//
//union double_to_ints {
//    CALC_TYPE double_value;
//    struct double_uint int_values;
//};
//
//void print_value_rf(CALC_TYPE d_value, char *buffer) {
//    union double_to_ints converter;
//    converter.double_value = d_value;
//    io_printf(
//        buffer, "0x%08x%08x",
//        converter.int_values.second_word, converter.int_values.first_word);
//}

enum regions {
	// add a recording region if required
	PARAMETERS
};

struct cholesky_parameters {

    // Acknowledge key for parameter location in SDRAM
    uint32_t acknowledge_key;

};

// The general parameters
struct cholesky_parameters cholesky_parameters;

// Acknowledge key global variable
uint32_t ack_key;

// Global variable for samples?
DataMat samples;

// Global variable for number of samples so far
uint32_t n_samples_read;

void trigger_run(uint key, uint payload);

// Helper functions for Cholesky decomposition go here
void zero_upper_triang(Mat A, uint32_t size) {
	uint32_t i, j;

	FOR( i, size )
		for( j = i+1; j < size; j++ )
			A[i][j] = 0.0f;
}

void cholesky(Mat A, const uint32_t size, bool zero_upper) {
	uint32_t i, j;
	int32_t k;
	float sum;

	log_info("Cholesky: A[0][0] %k", (accum) A[0][0]);
	log_info("Cholesky: A[10][10] %k", (accum) A[10][10]);
	log_info("Cholesky: A[10][15] %k", (accum) A[10][15]);

	log_info("Cholesky: size is %d", size);

	FOR( i, size ) {
		//log_info("i=%d", i);
		for( j = i; j < size; j++ ) {
			//log_info("j=%d", j);
			//sum = A[i][j];
			for( sum = A[i][j], k = i-1; k >= 0; k-- ) {
				//log_info("k=%d, sum=%k", k, (accum) sum);
				sum -= A[i][k] * A[j][k];
			}

			//log_info("i=%d j=%d, sum is %k", i, j, (accum) sum);
			if ( i == j ) {
				if ( sum <= 0.0f ) { // the other possibility here is to fail with an error
					log_info("Warning: possible non-pds matrix in cholesky()\n");
					A[i][i] = 1.0e-30f;  // i.e. very small positive value perhaps float epsilon * 10
				}
				else
					A[i][i] = sqrtf( sum );
			}
			else
				A[j][i] = sum / A[i][i];
		}
	}

	if( zero_upper )
		zero_upper_triang( A, size ); 	// zero upper triangle if necessary - for our application it is
}

//void vec_times_mat_scaled(const Vec restrict vec, const Mat restrict mat,
//		Vec res, const float scale, const uint32_t size) {
void vec_times_mat_scaled(const Vec vec, const Mat mat, Vec res,
		const float scale, const uint32_t size) {
	float sum;
	uint32_t i, j;

	FOR( i, size ) {
		sum = 0.0f;
		FOR( j, size )
			sum += vec[j] * mat[j][i];  // need to check if this correct way around! i.e. it could be mat[i][j]

		res[i] = sum * scale;
	}
}

void mean_covar_of_mat_n(const DataMat data, Vec mean, Mat cov,
		const uint32_t n, const uint32_t d) {
	uint32_t i, j, k;
	float 	sum, xi, xj, covar;

	log_info("Covariance matrix");

	FOR( i, d ) {									// calculate means
		for ( sum = 0.0f, j = 0; j < n; j++ )
			sum += data[j][i];

		mean[i] = sum / (float)n;
	}

	FOR( i, d ) {
		for ( j = i; j < d; j++ ) {
			covar = 0.0f;

			if( i == j )
				FOR( k, n ) {			// calculate variances on diagonal
					xi = data[k][i] - mean[i];
					covar += xi * xi;
				}
			else
				FOR( k, n ) {					// covariances off diagonal
					xi = data[k][i] - mean[i];
					xj = data[k][j] - mean[j];
					covar += xi * xj;
				}

			cov[i][j] = covar / ( (float)n - 1.0f ); // n must be > 1

			if( i != j )
				cov[j][i] = cov[i][j];

			log_info("i=%d j=%d, cov=%k", i, j, (accum) cov[i][j]);

		}
	}

}



// Sensible to get this size if it's directly available to us
uint32_t mcmc_model_get_params_n_bytes() {
    return sizeof(struct mcmc_params);
}
uint32_t mcmc_model_get_state_n_bytes() {
    return sizeof(struct mcmc_state);
}

void t_jump_callback(uint key, uint payload) {
	use(key);
	t_jump_ptr[0] = payload;
	log_info("Cholesky: t_jump_callback payload %d", payload);
	//spin1_callback_off(TIMER_TICK);
}

//void call_run(uint key, uint payload) {
//	use(key);
//	parameter_rec_ptr[0] = payload;
//	log_info("Cholesky: call_run, payload %d", payload);
//	run(0,0);
//}

// Function to collect parameter history and run Cholesky algorithm
void run(uint unused0, uint unused1) {
    use(unused0);
    use(unused1);

    // for debug writing values
 //   char buffer[1024];
    bool zero_upper = true;

    uint32_t params_n_bytes, state_n_bytes, n;
    uint32_t p, q, i;
    uint status;
    Vec mean, rot_scaled_t_jump, t_jump;
    Mat cov;

    // p and q defined in header
    p = PPOLYORDER;
    q = QPOLYORDER;

    n = p+q+2;

    CALC_TYPE state_parameters[n];


    // Get the size of the state parameters
	//uint32_t state_n_bytes = mcmc_model_get_state_n_bytes();

	// Get the parameters from SDRAM (from the correct place!)
	// the problem here is that these haven't been stored each time... ! :-)
	// But they could be: DTCM per core is 64K, so we can store approx
	//					  65536 / (8 * 20) ~= 410 samples, minus a bit for
	//			          other stored information. I think this is how to do it
	//spin1_memcpy(state_parameters, parameter_rec_ptr[0], state_n_bytes);

    log_info("Cholesky: n_samples_read is %d %d", n_samples_read, NCOVSAMPLES);

    log_info("Cholesky: pointer values %d %d", parameter_rec_ptr[0],
    		t_jump_ptr[0]);

    // Build up the sample
	state_n_bytes = mcmc_model_get_state_n_bytes();
	spin1_memcpy(state_parameters, parameter_rec_ptr[0], state_n_bytes);

	for (i=0; i<n; i++) {
		samples[n_samples_read][i] = state_parameters[i];
	}

	// send ptr back to the main vertex?
	while (!spin1_send_mc_packet(ack_key, parameter_rec_ptr[0],
			WITH_PAYLOAD)) {
		spin1_delay_us(1);
	}

	// We have now read in this sample
	n_samples_read++;

	// Do the work here every NCOVSAMPLES timesteps
	if (n_samples_read == NCOVSAMPLES) {

//		status = spin1_int_disable();
		spin1_callback_off(MCPL_PACKET_RECEIVED);

		// register to get the correct address for the jump scale parameters
	    spin1_callback_on(MCPL_PACKET_RECEIVED, t_jump_callback, -1);

	    log_info("Cholesky: pointer values %d %d", parameter_rec_ptr[0],
	    		t_jump_ptr[0]);

	    // We have reached the point where we want to do the calculation
		log_info("Cholesky: doing the calculation %d", NCOVSAMPLES);
		log_info("Cholesky: samples[0][0] %k", (accum) samples[0][0]);
		log_info("Cholesky: samples[1][0] %k", (accum) samples[1][0]);
		log_info("Cholesky: samples[2][0] %k", (accum) samples[2][0]);
		log_info("Cholesky: samples[5][0] %k", (accum) samples[5][0]);
		log_info("Cholesky: samples[15][0] %k", (accum) samples[15][0]);
		log_info("Cholesky: samples[95][0] %k", (accum) samples[95][0]);
		log_info("Cholesky: samples[99][0] %k", (accum) samples[99][0]);

		// Get the mean and covariance of the samples matrix
		mean_covar_of_mat_n(samples, mean, cov, n_samples_read, n);

		log_info("Cholesky: decompostion time, cov[0][0] %k",
				(accum) cov[0][0]);
		log_info("Cholesky: decompostion time, cov[1][1] %k",
				(accum) cov[1][1]);
		log_info("Cholesky: decompostion time, cov[10][10] %k",
				(accum) cov[10][10]);

		// Do the Cholesky decomposition
		cholesky(cov, n, zero_upper);

		log_info("Cholesky: register for message %d", NCOVSAMPLES);

	    log_info("Cholesky: pointer values %d %d", parameter_rec_ptr[0],
	    		t_jump_ptr[0]);

		// get the "old" t_jump from memory
		params_n_bytes = mcmc_model_get_params_n_bytes();
		spin1_memcpy(t_jump, t_jump_ptr[0], params_n_bytes);

		// Should probably check here what's in t_jump
		log_info("t_jump[0] %k", (accum) t_jump[0]);

		// generate a new t_jump vector
		vec_times_mat_scaled(t_jump, cov, rot_scaled_t_jump, 0.25, n);

		log_info("Cholesky: vec_times_mat_scaled, rot_scaled_t_jump[0]=%k",
				(accum) rot_scaled_t_jump[0]);

		// Copy this to relevant location
		spin1_memcpy(t_jump_ptr[0], rot_scaled_t_jump, params_n_bytes);

		// send ptr back to the main vertex?
		while (!spin1_send_mc_packet(ack_key, t_jump_ptr[0],
				WITH_PAYLOAD)) {
			spin1_delay_us(1);
		}

		log_info("Cholesky: after ptr sent back to ARMA");

		// Reset number of samples back to original value
		n_samples_read = 0;

		// turn this callback off somehow...
		spin1_callback_off(MCPL_PACKET_RECEIVED);

		spin1_callback_on(MCPL_PACKET_RECEIVED, trigger_run, -1);

//		spin1_mode_restore(status);

		// or something like this anyway...
	}

	// End of required functions

}

void trigger_run(uint key, uint payload) {
	use(key);
	// Get the pointer value to the location in SDRAM
	parameter_rec_ptr[0] = payload;
	log_info("Cholesky: trigger_run payload %d", payload);
	// Get ready to run the cholesky algorithm
	spin1_callback_off(TIMER_TICK);
    spin1_schedule_callback(run, 0, 0, 2);
}

void end_callback(uint unused0, uint unused1) {
	use(unused0);
	use(unused1);
	// End message has arrived from other vertex, so exit
	log_info("Cholesky: exit");
	spin1_exit(0);

}

void c_main() {
	// Get the acknowledge key from rf_parameters
	address_t data_address = data_specification_get_data_address();
	address_t cholesky_parameters_address = data_specification_get_region(
	        PARAMETERS, data_address);  // is this PARAMETERS or RECORDING?
										// Do we need both here?
	struct cholesky_parameters *cholesky_sdram_params =
			(struct cholesky_parameters *) cholesky_parameters_address;
	spin1_memcpy(&cholesky_parameters, cholesky_sdram_params,
			sizeof(struct cholesky_parameters));

	ack_key = cholesky_parameters.acknowledge_key;

	log_info("ack_key = 0x%08x", ack_key);

	// Set the global variable
	n_samples_read = 0;

	// register for the start message ?
    spin1_callback_on(MCPL_PACKET_RECEIVED, trigger_run, -1);

    // register for the shutdown message
    spin1_callback_on(MC_PACKET_RECEIVED, end_callback, -1);

	// start in sync_wait
    spin1_start(SYNC_WAIT);
}
