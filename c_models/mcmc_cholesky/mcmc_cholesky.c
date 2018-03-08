#include "mcmc_cholesky.h"
#include "../mcmc_models/mcmc_model.h"

#include <spin1_api.h>
#include <stdint.h>
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
	RECORDING, // I think this is probably necessary as it is where the history
	           // data that we need to use will be stored... ?
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

// Helper functions for Cholesky decomposition go here
void zero_upper_triang(Mat A, uint32_t size) {
	uint32_t i, j;

	FOR( i, size )
		for( j = i+1; j < size; j++ )
			A[i][j] = 0.0f;
}

void cholesky(Mat A, const uint32_t size, bool zero_upper) {
	uint32_t i, j, k;
	float sum;

	FOR( i, size ) {
		for( j = i; j < size; j++ ) {
			for( sum = A[i][j], k = i-1; k >= 0; k-- )
				sum -= A[i][k] * A[j][k];

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

		}
	}

}



// Sensible to get this size if it's directly available to us
uint32_t mcmc_model_get_state_n_bytes() {
    return sizeof(struct mcmc_state);
}

// Function to collect parameter history and run Cholesky algorithm
void run(uint unused0, uint unused1) {
    use(unused0);
    use(unused1);

    // for debug writing values
    char buffer[1024];

    uint8_t p, q, n, i;
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

	// Do the work here
	if (n_samples_read == NCOVSAMPLES) {
		// We have reached the point where we want to do the calculation

		// Get the mean and covariance of the samples matrix
		mean_covar_of_mat_n(samples, mean, cov, n_samples_read, n);

		// Do the Cholesky decomposition
		cholesky(cov, n, true);

		// get the "old" t_jump from memory, somehow: remember it's the
		// same size as the model state parameters
		uint32_t state_n_bytes = mcmc_model_get_state_n_bytes();
		spin1_memcpy(t_jump, t_jump_ptr[0], state_n_bytes);

		// Should probably check here what's in t_jump

		// generate a new t_jump vector
		vec_times_mat_scaled(t_jump, cov, rot_scaled_t_jump, 0.25, n);

		// Copy this to relevant location
		spin1_memcpy(t_jump_ptr[0], rot_scaled_t_jump, state_n_bytes);

		// send ptr back to the main vertex?

		// Reset number of samples back to zero
		n_samples_read = 0;

		// or something like this anyway...
	}
	else {
		// Build up the sample
		uint32_t state_n_bytes = mcmc_model_get_state_n_bytes();
		spin1_memcpy(state_parameters, parameter_rec_ptr[0], state_n_bytes);

		for (i=0; i<n; i++) {
			samples[n_samples_read][i] = state_parameters[i];
		}

		// We have now read in this sample
		n_samples_read++;
	}


	// At this point send the calculated Cholesky matrix back to the
	// model transition jump function in arma.c ?
	// Or is it better to do the matrix * parameters calc here and send
	// that result back?  ITCM may decide this...
	// wait here until packet is acknowledged/sent...
	float returnval = ONE;
	while (!spin1_send_mc_packet(ack_key, returnval, WITH_PAYLOAD)) {
		spin1_delay_us(1);
	}

	// End of required functions

}

void trigger_run(uint key, uint payload) {
	use(key);
	// Get the pointer value to the location in SDRAM
	parameter_rec_ptr[0] = payload;
	// Get ready to run the root_finder algorithm
	spin1_callback_off(TIMER_TICK);
    spin1_schedule_callback(run, 0, 0, 2);
}

void end_callback(uint unused0, uint unused1) {
	use(unused0);
	use(unused1);
	// End message has arrived from other vertex, so exit
//	log_info("Root finder: exit");
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

    // register for message from arma about location of t_jump vector ?

    // register for the shutdown message
    spin1_callback_on(MC_PACKET_RECEIVED, end_callback, -1);

	// start in sync_wait
    spin1_start(SYNC_WAIT);
}
