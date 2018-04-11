#include "mcmc_cholesky.h"
#include "../mcmc_models/mcmc_model.h"

#include <spin1_api.h>
#include <stdint.h>
#include <stdbool.h>
#include <debug.h>
#include <data_specification.h>


// Value of the pointer to the location in SDRAM to get parameters
uint32_t *parameter_rec_ptr;  // will need something like this...

uint32_t *t_variate_ptr;

uint32_t t_var_global;

// Functions required to print floats as hex values (uncomment for debug)
struct double_uint {
    uint first_word;
    uint second_word;
};

union double_to_ints {
    CALC_TYPE double_value;
    struct double_uint int_values;
};

void print_value_ch(CALC_TYPE d_value, char *buffer) {
    union double_to_ints converter;
    converter.double_value = d_value;
    io_printf(
        buffer, "0x%08x%08x",
        converter.int_values.second_word, converter.int_values.first_word);
}

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

// Global variable for samples? Needs a rethink for large NCOVSAMPLES > 800
//DataMat samples;

float** sample_data;

//float *dtcm_params;

// Global variable for number of samples so far
uint32_t n_samples_read;

// Global variable to use dtcm or sdram
uint32_t sdram_cholesky;

// Global variable for which calculation to use
bool do_calculation_using_cov_matrix;

void trigger_run(uint key, uint payload);

// Helper functions for Cholesky decomposition go here
void zero_upper_triang(Mat A, uint32_t size) {
	uint32_t i, j;

	FOR( i, size )
		for( j = i+1; j < size; j++ )
			A[i][j] = LA_ZERO;
}

void cholesky(Mat A, const uint32_t size, bool zero_upper) {
	uint32_t i, j;
	int32_t k;
	LA_TYPE sum;

	FOR( i, size ) {
		for( j = i; j < size; j++ ) {
			for( sum = A[i][j], k = i-1; k >= 0; k-- ) {
				sum -= A[i][k] * A[j][k];
			}

			if ( i == j ) {
				// 1) try everything with doubles (reduce small positive value accordingly)

				// 2) output more diagnostics (A matrix, i, etc.?)

				// 3) fail by exiting at this point and spitting out information
				//    (with some method to tell the gatherer not to include this sample)
				if ( sum <= LA_ZERO ) { // the other possibility here is to fail with an error
					log_info("Warning: possible non-pds matrix in cholesky()\n");
					A[i][i] = LA_SMALL;  // i.e. very small positive value perhaps float epsilon * 10
				}
				else
					A[i][i] = LA_SQRT( sum );
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
		const LA_TYPE scale, const uint32_t size) {
	LA_TYPE sum;
	uint32_t i, j;

	FOR( i, size ) {
		sum = LA_ZERO;
		FOR( j, size )
			sum += vec[j] * mat[j][i];  // need to check if this correct way around! i.e. it could be mat[i][j]

		res[i] = sum * scale;
	}
}

void mean_covar_of_mat_n(float **data, Vec mean, Mat cov,
		const uint32_t n, const uint32_t d) {
	uint32_t i, j, k;
	LA_TYPE sum, xi, xj, covar;

//	char buffer[1024];

//	print_value_ch(data[0][0], buffer);
//	log_info("mean_covar: data[0][0] = %k, hex is %s", (accum) data[0][0],
//			buffer);
//	print_value_ch(data[10][10], buffer);
//	log_info("mean_covar: data[10][10] = %k, hex is %s", (accum) data[10][10],
//			buffer);

	FOR( i, d ) {									// calculate means
		for ( sum = LA_ZERO, j = 0; j < n; j++ )
			sum += data[j][i];

		mean[i] = sum / (LA_TYPE)n;

	}

	FOR( i, d ) {
		for ( j = i; j < d; j++ ) {
			covar = LA_ZERO;

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

			cov[i][j] = covar / ( (LA_TYPE)n - LA_ONE ); // n must be > 1

			if( i != j )
				cov[j][i] = cov[i][j];

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

void t_variate_callback(uint key, uint payload) {
	use(key);
	t_variate_ptr[0] = payload;
	t_var_global = payload;
}

// Function to collect parameter history and run Cholesky algorithm
void run(uint unused0, uint unused1) {
    use(unused0);
    use(unused1);

    char buffer[1024];
    bool zero_upper = true;

    uint32_t params_n_bytes, state_n_bytes, n;
    uint32_t p, q, i, ii;
    Vec mean, rot_scaled_t_variate, t_variate;
    Mat cov;

    // p and q defined in header
    p = PPOLYORDER;
    q = QPOLYORDER;

    n = p+q+2;

    CALC_TYPE state_parameters[n];
    CALC_TYPE t_variate_data[n];
    CALC_TYPE rot_t_variate_data[n];

    // Build up the sample
	state_n_bytes = mcmc_model_get_state_n_bytes();
	spin1_memcpy(state_parameters, parameter_rec_ptr[0], state_n_bytes);


//	// samples is going to be what size now?
//	sample_data = (CALC_TYPE *) spin1_malloc(1*state_n_bytes);
//	log_info("sample data(1) is %d", sample_data);
//
//	sample_data = (CALC_TYPE *) spin1_malloc(100*state_n_bytes);
//	log_info("sample data(100) is %d", sample_data);
//
//	sample_data = (CALC_TYPE *) spin1_malloc(1000*state_n_bytes);
//	log_info("sample data(1000) is %d", sample_data);
//
//	sample_data = (CALC_TYPE *) spin1_malloc(5000*state_n_bytes);
//	log_info("sample data(5000) is %d", sample_data);


	// Store this sample data in DTCM if possible, but if the size of this
	// vector is going to end up bigger than DTCM can handle
	// (i.e. NCOVSAMPLES is (approx) greater than 65536 / (20 * 4) = 819),
	// which MH advises it should be (the Matlab code used NCOVSAMPLES=5000),
	// then we need to use DMAs/SDRAM instead.
	// allocate here? sample_data[n_samples_read] = (CALC_TYPE *)
	//sample_data[n_samples_read] = (CALC_TYPE *) spin1_malloc(state_n_bytes);
//	for (i=0; i<n; i++) {
//		dtcm_params[i] = state_parameters[i];
//	}
//	spin1_dma_transfer(DMA_WRITE, sample_data[n_samples_read],
//			dtcm_params, DMA_WRITE, state_n_bytes);
	for (i=0; i<n; i++) {
		sample_data[n_samples_read][i] = state_parameters[i];
	}

//	log_info("Cholesky, n_samples read: %d", n_samples_read);

	// send ptr back to the main vertex?
	while (!spin1_send_mc_packet(ack_key, parameter_rec_ptr[0],
			WITH_PAYLOAD)) {
		spin1_delay_us(1);
	}

//	print_value_ch(sample_data[n_samples_read][1], buffer);
//	log_info("Cholesky: sample_data[%d][1] = %k, hex is %s", n_samples_read,
//			(accum) sample_data[n_samples_read][1], buffer);

	// We have now read in this sample
	n_samples_read++;

	// Set global var
	t_var_global = 2;

	// Get ready to read the t_variate vector next
	spin1_callback_off(MCPL_PACKET_RECEIVED);

	// register to get the correct address for the parameters
	spin1_callback_on(MCPL_PACKET_RECEIVED, t_variate_callback, -1);

	// wait until this arrives
	while (t_var_global == 2) {
		spin1_wfi();
	}

	// Do the work here every NCOVSAMPLES timesteps
	if (n_samples_read == NCOVSAMPLES) {

		log_info("Performing Cholesky decomposition");

		// Get the mean and covariance of the samples matrix
		mean_covar_of_mat_n(sample_data, mean, cov, n_samples_read, n);

		// Do the Cholesky decomposition
		cholesky(cov, n, zero_upper);

		// get the t_variate from memory
		params_n_bytes = mcmc_model_get_params_n_bytes();
		// Data that arrives here is CALC_TYPE (i.e. float)
		spin1_memcpy(t_variate_data, t_variate_ptr[0], params_n_bytes);
		// Convert to LA_TYPE
		for (i=0; i<n; i++) {
			t_variate[i] = (LA_TYPE) t_variate_data[i];
		}

		// generate a new t_variate vector
		vec_times_mat_scaled(t_variate, cov, rot_scaled_t_variate, 0.25f, n);

		// Convert back to CALC_TYPE again
		for (i=0; i<n; i++) {
			rot_t_variate_data[i] = (CALC_TYPE) rot_scaled_t_variate[i];
		}

		// Copy this to relevant location
		spin1_memcpy(t_variate_ptr[0], rot_t_variate_data, params_n_bytes);

		// send ptr back to the main vertex?
		while (!spin1_send_mc_packet(ack_key, t_variate_ptr[0],
				WITH_PAYLOAD)) {
			spin1_delay_us(1);
		}

		// Reset number of samples back to original value
		n_samples_read = 0;

		// Tell the code that it now needs to do the calculation
		// of the new jump scale vector using the covariance matrix
		do_calculation_using_cov_matrix = true;

	}
	else {
		// Not the Nth timestep, so do the calculation and return values

		// Get the t_variate from memory
		params_n_bytes = mcmc_model_get_params_n_bytes();
		spin1_memcpy(t_variate_data, t_variate_ptr[0], params_n_bytes);

		// If we are over N timesteps then this uses the covariance matrix
		if (do_calculation_using_cov_matrix) {
			// Convert to LA_TYPE
			for (i=0; i<n; i++) {
				t_variate[i] = (LA_TYPE) t_variate_data[i];
			}

			// generate a new t_variate vector using covariance
			vec_times_mat_scaled(t_variate, cov,
					rot_scaled_t_variate, LA_QUARTER, n);

			// Convert back to CALC_TYPE again
			for (i=0; i<n; i++) {
				rot_t_variate_data[i] = (CALC_TYPE) rot_scaled_t_variate[i];
			}
		}
		else {
			// We are under N timesteps overall so this is a vector*vector
			// (element-wise) calculation instead... but it can't be done here
			// so send it back for the main ARMA vertex to deal with

			for (ii=0; ii<n; ii++) {
				rot_t_variate_data[ii] = t_variate_data[ii];
			}
		}

		// Copy this to relevant location
		spin1_memcpy(t_variate_ptr[0], rot_t_variate_data, params_n_bytes);

		// send ptr back to the main vertex?
		while (!spin1_send_mc_packet(ack_key, t_variate_ptr[0],
				WITH_PAYLOAD)) {
			spin1_delay_us(1);
		}
	}

	// Turn the t_variate callback off
	spin1_callback_off(MCPL_PACKET_RECEIVED);

	// Turn the run trigger callback back on again
	spin1_callback_on(MCPL_PACKET_RECEIVED, trigger_run, -1);
}

void trigger_run(uint key, uint payload) {
	use(key);
	// Get the pointer value to the location in SDRAM
	parameter_rec_ptr[0] = payload;
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
	        PARAMETERS, data_address);
	struct cholesky_parameters *cholesky_sdram_params =
			(struct cholesky_parameters *) cholesky_parameters_address;
	spin1_memcpy(&cholesky_parameters, cholesky_sdram_params,
			sizeof(struct cholesky_parameters));

	ack_key = cholesky_parameters.acknowledge_key;

	log_info("ack_key = 0x%08x", ack_key);

	// sample_data is going to be what size now?
	uint32_t ii;
	uint32_t state_n_bytes = mcmc_model_get_state_n_bytes();
//	sample_data = (CALC_TYPE **) spin1_malloc(NCOVSAMPLES*state_n_bytes);
//
//	if (sample_data!=NULL) {
//		// we can do this using dtcm
//		log_info("Doing calculation using DTCM: %d", sample_data);
//		//sample_data = (CALC_TYPE **)
//		sdram_cholesky = 0;
//	}
//	else {
		// we can't do it using dtcm
	uint32_t n_samples = NCOVSAMPLES;
	sample_data = (CALC_TYPE **) sark_xalloc(sv->sdram_heap,
			n_samples*sizeof(CALC_TYPE*), 0, ALLOC_LOCK);
	log_info("calculation using SDRAM: %d", sample_data);
	//uint32_t param_size = PPOLYORDER+QPOLYORDER+2;
	sample_data[0] = (CALC_TYPE *) sark_xalloc(sv->sdram_heap,
			n_samples*state_n_bytes, 0, ALLOC_LOCK);
	log_info("sample_data[0]: %d", sample_data[0]);
	//dtcm_params = (CALC_TYPE *) sark_alloc(param_size, sizeof(CALC_TYPE));
	// Set pointers for the remaining rows
	for (ii=1; ii<n_samples; ii++) {
		sample_data[ii] = sample_data[ii-1]+(PPOLYORDER+QPOLYORDER+2);
	}
	sdram_cholesky = 1;
//	}

	// Set the global variables
	n_samples_read = 0;
	do_calculation_using_cov_matrix = false;

	// register for the start message ?
    spin1_callback_on(MCPL_PACKET_RECEIVED, trigger_run, -1);

    // register for the shutdown message
    spin1_callback_on(MC_PACKET_RECEIVED, end_callback, -1);

	// start in sync_wait
    spin1_start(SYNC_WAIT);
}
