/*
 * Copyright (c) 2016-2021 The University of Manchester
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
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
// compile with gcc -std=c99 -O3 -fcx-limited-range root_finder.c -l m -o rf

#include "mcmc_root_finder.h"

#include <mcmc_model.h>
#include <arma.h>

#include <spin1_api.h>
#include <stdint.h>
#include <debug.h>
#include <data_specification.h>

#define ROOT_FAIL 1 // define it this way here to send "boolean" back

// Value of the location in SDRAM to get parameters
uint32_t parameter_rec_add;

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
	// add a recording region for debug to check the roots calculation if required
	PARAMETERS
};

struct rf_parameters {

    // Acknowledge key for parameter location in SDRAM
    uint32_t acknowledge_key;

};

// The general parameters
struct rf_parameters rf_parameters;

// Define spin1_wfi
extern void spin1_wfi();

// Acknowledge key global variable
uint32_t ack_key;

// Helper functions for finding the root

// multiply complex by scalar
complex float RCmul( float x, complex float a )
{
	return ( x * crealf( a )) + ( x * cimagf( a )) * I;
}

// edit these constants if necessary for fixed point
// (Note, these could plausibly be passed in from the Python frontend)
// (Note also, it's probably a good idea to replace all "float" with "CALC_TYPE")
#define EPSS 1.0e-7
#define MR 8
#define MT 10
#define MAXIT (MT*MR)
/*
	Given degree m and the m+1 complex coefficients a[0..m] of the polynomial a[i]*x^i
	and a complex value x, this function improves x by Laguerre's method until it
	converges - within the achievable roundoff limit - to a root of the given
	polynomial. Number of iterations taken is returned as its
*/
void laguerre_poly_root( complex float a[], int m, complex float *x, int *its )
{
	int 				iter, j;
	float 			abx, abp, abm, err;
	complex float	dx, x1, b, d, f, g, h, sq, gp, gm, g2;
	static float 	frac[MR+1] = { 0.0, 0.5, 0.25, 0.75, 0.13, 0.38, 0.62, 0.88, 1.0 };

	for ( iter = 1; iter <= MAXIT; iter++ ) {  // loop over iterations up to maximum
		*its = iter;
		b = a[m];
		err = cabsf( b );
		d = f = 0.0f + 0.0f * I;
		//d = f = ZERO * ZERO * I;
		abx = cabsf( *x );
		for ( j = m-1; j >= 0; j-- ) {	// efficient computation of polynomial and first two derivatives
			f = ( *x * f ) + d;
			d = ( *x * d ) + b;
			b = ( *x * b ) + a[j];
			err = cabsf( b ) + abx * err;
			}

		err *= EPSS;	// estimate of roundoff error in evaluating polynomial

		if ( cabsf( b ) <= err ) {
			return;	// we are on the root
		}

		g = d / b;								// the generic case so use Laguerre's formula
		g2 = g * g;
		//h = g2 - RCmul( TWO, f / b );
		h = g2 - RCmul( 2.0f, f / b );
		sq = csqrtf( RCmul( (float) (m-1), RCmul((float) m, h ) - g2 ));
		gp = g + sq;
		gm = g - sq;
		abp = cabsf( gp );
		abm = cabsf( gm );

		if ( abp < abm ) {
			gp = gm;
		}

		dx = FMAX( abp, abm ) > 0.0f ? ( ((float) m ) + 0.0f * I ) / gp :
		RCmul( 1.0f + abx, cosf((float)iter) + sinf((float)iter) * I );
//		dx = FMAX( abp, abm ) > ZERO ? ( ((float) m ) + ZERO * I ) / gp :
//		RCmul( ONE + abx, cosf((float)iter) + sinf((float)iter) * I );

		x1 = *x - dx;

		if ( crealf(*x) == crealf(x1) && cimagf(*x) == cimagf(x1) ) {
			return; // converged
		}

		if ( iter % MT ) {
			*x = x1;
		}
		else {
			// occasionally take fractional step to break rare limit cycle
			*x = *x - RCmul( frac[ iter / MT ], dx );
		}
	}

	// Note: Not sure this is necessary, but could be a log_info to see if it ever occurs
	// log_info("Too many iterations in laguerre_poly_root()");
	// very unusual and only for complex roots
	// printf("Too many iterations in laguerre_poly_root()"); exit(1);

	return;
}

/* from NR in C errata
*** 40,42 ****
  		dx=((FMAX(abp,abm) > 0.0 ? Cdiv(Complex((float) m,0.0),gp)
! 			: RCmul(exp(log(1+abx)),Complex(cos((float)iter),sin((float)iter)))));
  		x1=Csub(*x,dx);
--- 40,42 ----
  		dx=((FMAX(abp,abm) > 0.0 ? Cdiv(Complex((float) m,0.0),gp)
! 			: RCmul(1+abx,Complex(cos((float)iter),sin((float)iter)))));
  		x1=Csub(*x,dx);

*/

// This is the actual function that does the work

#define EPS 2.0e-6
#define MAXM 100
/*
	Given the degree m and m+1 complex coefficients a[0..m] of the polynomial a[i]*x^i
	this function successively calls laguerre_poly_root() and finds all m complex
	roots in roots[1..m]. The bool polish should be input as true if polishing is
	desired (i.e. almost always)
*/
void zroots( complex float a[], int m, complex float roots[], bool polish )
{
	int 				i, its, j, jj;
	complex float 	x, b, c, ad[MAXM];

	for ( j = 0; j <= m; j++ ) {
		ad[j] = a[j]; // copy coeffs for successive deflation
	}

	for ( j = m; j >= 1; j-- ) {	// loop over each root to be found

		//x = ZERO + ZERO * I;  		// start at zero to favour convergence to smallest root
		x = 0.0f + 0.0f * I;  		// start at zero to favour convergence to smallest root

		laguerre_poly_root( ad, j, &x, &its );

		if ( fabs(cimagf(x)) <= 2.0f * EPS * fabs(crealf(x))) {
			x = crealf(x) + 0.0f * I; // set imaginary part to zero
		}
//		if ( fabs(cimagf(x)) <= TWO * EPS * fabs(crealf(x))) {
//			x = crealf(x) + ZERO * I; // set imaginary part to zero
//		}

		roots[j] = x;

		b = ad[j]; // forward deflation

		for ( jj = j-1; jj>=0; jj-- ) {
			c = ad[jj];
			ad[jj] = b;
			b = (x * b) + c;
		}
	}

	if ( polish ) {
		for ( j = 1; j <= m; j++ ) {
			// polish roots using undeflated coeffs
			laguerre_poly_root( a, m, &roots[j], &its );
		}
	}

	for ( j = 2; j <= m; j++ ) {		// sort roots by their real parts

		x = roots[j];

		for ( i = j-1; i >= 1; i-- ) {
			if ( crealf(roots[i]) <= crealf(x) ) {
				break;
			}
			roots[i+1] = roots[i];
		}

		roots[i+1] = x;
	}
}

// Sensible to get this size if it's directly available to us
uint32_t mcmc_model_get_state_n_bytes(void) {
    return sizeof(struct mcmc_state);
}

// Function to collect parameters and run root_finder algorithm
void run(uint unused0, uint unused1) {
    use(unused0);
    use(unused1);

    // for debug writing values
//    char buffer[1024];

    uint32_t i, p, q;

    // p and q defined in header
    p = PPOLYORDER;
    q = QPOLYORDER;

    // allocate some DTCM for this array
	CALC_TYPE *state_parameters = spin1_malloc((p+q+2)*sizeof(CALC_TYPE));
//	CALC_TYPE state_parameters[p+q+2];

//    log_info("ROOTFINDER: running root finder");

    // Get the size of the state parameters
	uint32_t state_n_bytes = mcmc_model_get_state_n_bytes();

	// Get the parameters from SDRAM
	uint32_t *param_ptr = (uint32_t *) parameter_rec_add;
	spin1_memcpy(state_parameters, param_ptr, state_n_bytes);

	// Set up data structures for the coefficients of a polynomial
	// characteristic equation for AR and MA model respectively.
	// C99 compile flag required for complex type and variable length arrays
	// Trying to do this using memory allocation rather than on the fly
	CALC_TYPE *AR_eq = spin1_malloc((p+1)*sizeof(CALC_TYPE));
	CALC_TYPE *MA_eq = spin1_malloc((q+1)*sizeof(CALC_TYPE));

	//log_info("ROOTFINDER: sizeof(complex float)=%u", sizeof(complex float));
	complex float *AR_param = spin1_malloc((p+1)*sizeof(complex float));
	complex float *MA_param = spin1_malloc((q+1)*sizeof(complex float));
	complex float *AR_rt = spin1_malloc((p+1)*sizeof(complex float));
	complex float *MA_rt = spin1_malloc((q+1)*sizeof(complex float));
	//complex float AR_param[p+1], MA_param[q+1], AR_rt[p+1], MA_rt[q+1];

	// Create an array with the coefficients for the characteristic AR equation
	// The characteristic equation is -a_p*x^p-a_(p-1)*x^(p-1)-...+1=0;

	// zroots needs them in order from lowest power to highest
	AR_eq[0] = ONE;
	for (i=0; i < p; i++) {
		AR_eq[i+1] = -state_parameters[i];
	}

	// Create an array with the coefficients for the characteristic MA equation
	// The characteristic equation is -b_q*x^q-b_(q-1)*x^(q-1)-...+1=0;
	MA_eq[0] = ONE;
	for (i=0; i < q; i++) {
		MA_eq[i+1] = -state_parameters[p+i];
	}

	// read parameters into complex vectors so that we can calculate roots - from 0?
	for (i=0; i <= p; i++) {
		// only real part relevant - so imaginary part = 0
		AR_param[i] = (float)AR_eq[i] + ZERO * I;
	}

	for (i=0; i <= q; i++ ) {
		// only real part relevant - so imaginary part = 0
		MA_param[i] = (float)MA_eq[i] + ZERO * I;
	}

	// These calls find the complex roots for each case
	zroots( AR_param, p, AR_rt, true);
	zroots( MA_param, q, MA_rt, true);

	// Set the return value
	uint32_t returnval = 0;

	// test for root magnitude <= 1 and if so return a fail result
	// zroots returns values from array index 1 upwards
	for (i=1; i <= p; i++) {
		if ( cabsf(AR_rt[i]) <= ONE ) { returnval = ROOT_FAIL; }
	}

	for (i=1; i <= q; i++) {
		if ( cabsf(MA_rt[i]) <= ONE ) { returnval = ROOT_FAIL; }
	}

	// At this point send returnval to normal ARMA vertex:
	// wait here until packet is acknowledged/sent...
	while (!spin1_send_mc_packet(ack_key, returnval, WITH_PAYLOAD)) {
		spin1_delay_us(1);
	}

	// free up memory
	sark_free(state_parameters);
	sark_free(AR_eq);
	sark_free(MA_eq);
	sark_free(AR_param);
	sark_free(MA_param);
	sark_free(AR_rt);
	sark_free(MA_rt);

	// End of required functions

}

void trigger_run(uint key, uint payload) {
	use(key);
	// Get the value of the location in SDRAM
	parameter_rec_add = payload;
	// Get ready to run the root_finder algorithm
	spin1_callback_off(TIMER_TICK);
    spin1_schedule_callback(run, 0, 0, 2);
}

void end_callback(uint unused0, uint unused1) {
	use(unused0);
	use(unused1);
	// End message has arrived from other vertex, so exit
	log_info("Root finder: exit");
	spin1_exit(0);
}

void c_main(void) {
	// Get the acknowledge key from rf_parameters
	data_specification_metadata_t *data_address = data_specification_get_data_address();
	address_t rf_parameters_address = data_specification_get_region(
	        PARAMETERS, data_address);
	struct rf_parameters *rf_sdram_params =
			(struct rf_parameters *) rf_parameters_address;
	spin1_memcpy(&rf_parameters, rf_sdram_params,
			sizeof(struct rf_parameters));

	ack_key = rf_parameters.acknowledge_key;

	log_info("ROOTFINDER ack_key = 0x%08x", ack_key);

	// register for the start message
    spin1_callback_on(MCPL_PACKET_RECEIVED, trigger_run, -1);

    // register for the shutdown message
    spin1_callback_on(MC_PACKET_RECEIVED, end_callback, -1);

	// start in sync_wait
    spin1_start(SYNC_WAIT);
}

/*
	This is a test program - set up some polynomial coefficients and find their complex roots
	seems to works up to high orders but needs testing on some real ARMA coefficients
*/
//#define ORDER 18
//
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
