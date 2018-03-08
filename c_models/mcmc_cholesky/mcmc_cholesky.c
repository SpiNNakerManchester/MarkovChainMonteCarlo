#include "mcmc_cholesky.h"

#include <spin1_api.h>
#include <stdint.h>
#include <debug.h>
#include <data_specification.h>


// Value of the pointer to the location in SDRAM to get parameters
uint32_t *parameter_rec_ptr;  // will need something like this...

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
struct cholesky_parameters rf_parameters;

// Acknowledge key global variable
uint32_t ack_key;

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
					A[i][i] = 1.0e-30;  // i.e. very small positive value perhaps float epsilon * 10
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

void vec_times_mat_scaled(const Vec restrict vec, const Mat restrict mat,
		Vec res, const float scale, const uint32_t size) {
	float sum;
	uint32_t i, j;

	FOR( i, size ) {
		sum = 0.0f;
		FOR( j, size )
			sum += vec[j] * mat[j][i];  // need to check if this correct way around! i.e. it could be mat[i][j]

		res[i] = sum * scale;
	}
}

void mean_covar_of_mat_n(const Mat data, Vec mean, Mat cov,
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

    uint8_t i, p, q;
    CALC_TYPE state_parameters[PPOLYORDER+QPOLYORDER+2];

    // p and q defined in header
    p = PPOLYORDER;
    q = QPOLYORDER;

    // Get the size of the state parameters
	uint32_t state_n_bytes = mcmc_model_get_state_n_bytes();

	// Get the parameters from SDRAM (from the correct place!)
	spin1_memcpy(state_parameters, parameter_rec_ptr[0], state_n_bytes);

	// Do the work here

	// At this point send the calculated Cholesky matrix back to the
	// model transition jump function in arma.c ?
	// Or is it better to do the matrix * parameters calc here and send
	// that result back?  ITCM may decide this...
	// wait here until packet is acknowledged/sent...
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

	ack_key = rf_parameters.acknowledge_key;

	log_info("ack_key = 0x%08x", ack_key);

	// register for the start message
    spin1_callback_on(MCPL_PACKET_RECEIVED, trigger_run, -1);

    // register for the shutdown message
    spin1_callback_on(MC_PACKET_RECEIVED, end_callback, -1);

	// start in sync_wait
    spin1_start(SYNC_WAIT);
}

//#define ORDER 18

/*
	This is a test program - set up some polynomial coefficients and find their complex roots
	seems to works up to high orders but needs testing on some real ARMA coefficients
*/
//int test_main( void )
//int main( void )
//{
//	int  i, m = ORDER;
//	bool pos = true;  // alternator for coeff signs
//	complex float a[ORDER+1], roots[ORDER+1];
//
//	for( i = 0; i <= ORDER; i++ ) {  // set up some polynomial coefficients with alternating signs
//		if( pos ) {
//			a[i] = ( (float)(ORDER - i) + 0.5f ) + 0.0f * I;
//			pos = false;
//			}
//		else {
//			a[i] = -( (float)(ORDER - i) + 0.3f ) + 0.0f * I;
//			pos = true;
//			}
//		printf("\n real coeff of order:%3d  = %9.3f", i,  crealf(a[i]) );
//		}
//
//	printf("\n\n");
//
//	zroots( a, m, roots, true );  // this function does the actual work
//
//	for( i = 1; i <= ORDER; i++ )	// output complex roots
//		printf("\n root:%3d  real = %12.8f   imag = %12.8f", i,  crealf(roots[i]),  cimagf(roots[i]) );
//
//	printf("\n\n");
//
//	return 0;
//}
