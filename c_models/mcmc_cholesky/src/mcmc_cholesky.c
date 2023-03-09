/*
 * Copyright (c) 2016 The University of Manchester
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */
#include "mcmc_cholesky.h"
#include <mcmc_model.h>

#include <spin1_api.h>
#include <stdint.h>
#include <stdbool.h>
#include <debug.h>
#include <data_specification.h>

// Value of the locations in SDRAM to get parameters / t_variate data
uint32_t parameter_rec_add;
uint32_t t_variate_add;

// Value to be received when waiting for t_variate data
uint32_t t_var_global;

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
//void print_value_ch(LA_TYPE d_value, char *buffer) {
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

// Define spin1_wfi
extern void spin1_wfi();

// Acknowledge key global variable
uint32_t ack_key;

// Global variable for samples? Needs a rethink for large NCOVSAMPLES > 800
CALC_TYPE** sample_data;

// Global variable for covariance matrix (which gets reused each step after
// the first burn-in period, so it needs to be global!)
Mat covariance;

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

	for (i=0; i<size; i++) {
		for (j=i+1; j<size; j++) {
			A[i][j] = LA_ZERO;
		}
	}
}

void cholesky(Mat A, const uint32_t size, bool zero_upper) {
//	char buffer[1024];
	uint32_t i, j;
	int32_t k;
	LA_TYPE sum;

	for (i=0; i<size; i++) {
		for (j=i; j<size; j++) {
			for (sum=A[j][i], k=i-1; k>=0; k--) {
				sum -= A[j][k] * A[i][k];
			}

			if (i == j) {
				if (sum <= LA_ZERO) {
					// do we fail by exiting at this point and spitting out information?
					// (with some method to tell the gatherer not to include this sample)
					//log_info("Warning: possible non-pds matrix in cholesky()\n");
					io_printf(IO_BUF, "Warning: possible non-pds matrix in cholesky()\n");
					A[i][i] = LA_SMALL;  // i.e. very small positive value
				}
				else {
					A[i][i] = LA_SQRT(sum);
				}
			}
			else {
				A[j][i] = sum / A[i][i];
			}
		}
	}

	if (zero_upper) {
		zero_upper_triang(A, size); 	// zero upper triangle if necessary - for our application it is
	}
}

void vec_times_mat_scaled(const Vec vec, const Mat mat, Vec res,
		const LA_TYPE scale, const uint32_t size) {
	LA_TYPE sum;
	uint32_t i, j;

	for (i=0; i<size; i++) {
		sum = LA_ZERO;
		for (j=0; j<size; j++) {
			sum += vec[j] * mat[i][j];
		}
		// need to check if this correct way around! i.e. it could be mat[j][i]
		// (I think this is the correct way round: it's a bit confusing compared
		//  directly to MATLAB as the cholesky command there gives an upper-
		//  triangular matrix whereas this (correctly) uses a lower-triangular)

		res[i] = sum * scale;
	}
}

void mean_covar_of_mat_n(const uint32_t n, const uint32_t d) {
	uint32_t i, j, k;
	LA_TYPE sum, xi, xj, covar;
	Vec mean;

	for (i=0; i<d; i++) {									// calculate means
		for (sum=LA_ZERO, j=0; j<n; j++) {
			sum += sample_data[j][i];
		}

		mean[i] = sum / (LA_TYPE)n;

	}

	for (i=0; i<d; i++) {
		for (j=i; j<d; j++) {
			covar = LA_ZERO;

			if (i == j) {
				for (k=0; k<n; k++) {			// calculate variances on diagonal
					xi = sample_data[k][i] - mean[i];
					covar += xi * xi;
				}
			}
			else {
				for (k=0; k<n; k++) {					// covariances off diagonal
					xi = sample_data[k][i] - mean[i];
					xj = sample_data[k][j] - mean[j];
					covar += xi * xj;
				}
			}

			covariance[i][j] = covar / ((LA_TYPE)n-LA_ONE); // n must be > 1

			if (i != j) {
				covariance[j][i] = covariance[i][j];
			}
		}
	}
}

// Sensible to get this size if it's directly available to us
uint32_t mcmc_model_get_params_n_bytes(void) {
    return sizeof(struct mcmc_params);
}
uint32_t mcmc_model_get_state_n_bytes(void) {
    return sizeof(struct mcmc_state);
}

void t_variate_callback(uint key, uint payload) {
	use(key);
	t_variate_add = payload;
	t_var_global = payload;
}

// Function to collect parameter history and run Cholesky algorithm
void run(uint unused0, uint unused1) {
    use(unused0);
    use(unused1);

//    char buffer[1024];
    bool zero_upper = true;

    uint32_t params_n_bytes, state_n_bytes, n;
    uint32_t i, ii;
    Vec rot_scaled_t_variate, t_variate;

    // vector/matrix dimension defined in header
    n = MATDIM;

    CALC_TYPE state_parameters[n];
    CALC_TYPE t_variate_data[n];
    CALC_TYPE rot_t_variate_data[n];

    // Build up the sample
	state_n_bytes = mcmc_model_get_state_n_bytes();
	uint32_t *param_ptr = (uint32_t *) parameter_rec_add;
	spin1_memcpy(state_parameters, param_ptr, state_n_bytes);

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

	// The allocation is done in main(), I've left the useful comment for now
	for (i=0; i<n; i++) {
		sample_data[n_samples_read][i] = state_parameters[i];
	}

	// send ptr back to the main vertex
	while (!spin1_send_mc_packet(ack_key, parameter_rec_add,
			WITH_PAYLOAD)) {
		spin1_delay_us(1);
	}

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
		// Print out so we know it's got here
		io_printf(IO_BUF, "Performing Cholesky decomposition \n");

		// Get the mean and covariance of the samples matrix
		mean_covar_of_mat_n(n_samples_read, n);
		// Note: there's no need to "return" the mean here, as it isn't used

		// Do the Cholesky decomposition
		cholesky(covariance, n, zero_upper);

		// get the t_variate from memory
		params_n_bytes = mcmc_model_get_params_n_bytes();

		// Data that arrives here is CALC_TYPE (i.e. float)
		uint32_t *t_var_ptr = (uint32_t *) t_variate_add;
		spin1_memcpy(t_variate_data, t_var_ptr, params_n_bytes);

		// Convert to LA_TYPE
		for (i=0; i<n; i++) {
			t_variate[i] = (LA_TYPE) t_variate_data[i];
		}

		// generate a new t_variate vector
		vec_times_mat_scaled(t_variate, covariance,
				rot_scaled_t_variate, LA_SCALE, n);

		// Convert back to CALC_TYPE again
		for (i=0; i<n; i++) {
			rot_t_variate_data[i] = (CALC_TYPE) rot_scaled_t_variate[i];
		}

		// Copy this to relevant location
		spin1_memcpy(t_var_ptr, rot_t_variate_data, params_n_bytes);

		// send ptr back to the main vertex
		while (!spin1_send_mc_packet(ack_key, t_variate_add,
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
		uint32_t *t_var_ptr = (uint32_t *) t_variate_add;
		spin1_memcpy(t_variate_data, t_var_ptr, params_n_bytes);

		// If we are over N timesteps then this uses the covariance matrix
		if (do_calculation_using_cov_matrix) {
			// Convert to LA_TYPE
			for (i=0; i<n; i++) {
				t_variate[i] = (LA_TYPE) t_variate_data[i];
			}

			// generate a new t_variate vector using cholesky covariance
			vec_times_mat_scaled(t_variate, covariance,
					rot_scaled_t_variate, LA_SCALE, n);

			// Convert back to CALC_TYPE again
			for (i=0; i<n; i++) {
				rot_t_variate_data[i] = (CALC_TYPE) rot_scaled_t_variate[i];
			}
		}
		else {
			// We are under N timesteps overall so this is a vector*vector
			// (element-wise) calculation instead... but it can't be done here
			// so send the vector back for the main ARMA vertex to deal with
			for (ii=0; ii<n; ii++) {
				rot_t_variate_data[ii] = t_variate_data[ii];
			}
		}

		// Copy this to relevant location
		spin1_memcpy(t_var_ptr, rot_t_variate_data, params_n_bytes);

		// send ptr back to the main vertex
		while (!spin1_send_mc_packet(ack_key, t_variate_add,
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
	parameter_rec_add = payload;
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

void c_main(void) {

//	Testing Cholesky: note, you will need to alter the cholesky function to
//	take double[test_size][test_size] to make this one work instead of Mat
//
//  (Using examples from https://rosettacode.org/wiki/Cholesky_decomposition)
//	uint32_t i, j;
//	char buffer[1024];
//	bool zero_upper = true;
//
////	log_info("Cholesky test matrix (25 15 -5, 15 18 0, -5 0 11)");
//	log_info("Cholesky test matrix (18 22 54 42, 22 70 86 62, 54 86 174 134,",
//			" 42 62 134 106");
//
//	uint32_t test_size = 4; // 3;
//
//    double test_matrix[test_size][test_size];
//    test_matrix[0][0]=18.0;  // this one works, tested ADG 04/05/18
//    test_matrix[1][0]=22.0;
//    test_matrix[2][0]=54.0;
//    test_matrix[3][0]=42.0;
//    test_matrix[0][1]=22.0;
//    test_matrix[1][1]=70.0;
//    test_matrix[2][1]=86.0;
//    test_matrix[3][1]=62.0;
//    test_matrix[0][2]=54.0;
//    test_matrix[1][2]=86.0;
//    test_matrix[2][2]=174.0;
//    test_matrix[3][2]=134.0;
//    test_matrix[0][3]=42.0;
//    test_matrix[1][3]=62.0;
//    test_matrix[2][3]=134.0;
//    test_matrix[3][3]=106.0;
//
////    test_matrix[0][0]=25.0;  // this one works, tested ADG 02/05/18
////    test_matrix[1][0]=15.0;
////    test_matrix[2][0]=-5.0;
////    test_matrix[0][1]=15.0;
////    test_matrix[1][1]=18.0;
////    test_matrix[2][1]=0.0;
////    test_matrix[0][2]=-5.0;
////    test_matrix[1][2]=0.0;
////    test_matrix[2][2]=11.0;
//
//    // What's in test_matrix now?
//	for (i=0; i<test_size; i++) {
//		for (j=0; j<test_size; j++) {
//			print_value_ch(test_matrix[j][i], buffer);
//			log_info("pre test_matrix[%u][%u] = %s ", j, i, buffer);
//		}
//	}
//
//    // Cholesky on this
//    cholesky(test_matrix, test_size, zero_upper);
//
//    // What's in test_matrix now?
//	for (i=0; i<test_size; i++) {
//		for (j=0; j<test_size; j++) {
//			print_value_ch(test_matrix[j][i], buffer);
//			log_info("test_matrix[%u][%u] = %s ", j, i, buffer);
//		}
//	}

	// Get the acknowledge key from rf_parameters
	data_specification_metadata_t *data_address = data_specification_get_data_address();
	address_t cholesky_parameters_address = data_specification_get_region(
	        PARAMETERS, data_address);
	struct cholesky_parameters *cholesky_sdram_params =
			(struct cholesky_parameters *) cholesky_parameters_address;
	spin1_memcpy(&cholesky_parameters, cholesky_sdram_params,
			sizeof(struct cholesky_parameters));

	ack_key = cholesky_parameters.acknowledge_key;

	log_info("ack_key = 0x%08x", ack_key);

	// sample_data is going to be what size?
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
		sample_data[ii] = sample_data[ii-1]+(MATDIM);
	}
	sdram_cholesky = 1;
//	}

	// Set the global variables
	n_samples_read = 0;
	do_calculation_using_cov_matrix = false;

	// register for the start message
    spin1_callback_on(MCPL_PACKET_RECEIVED, trigger_run, -1);

    // register for the shutdown message
    spin1_callback_on(MC_PACKET_RECEIVED, end_callback, -1);

	// start in sync_wait
    spin1_start(SYNC_WAIT);
}
