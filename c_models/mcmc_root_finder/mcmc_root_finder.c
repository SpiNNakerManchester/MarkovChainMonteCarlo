
// compile with gcc -std=c99 -O3 -fcx-limited-range root_finder.c -l m -o rf

#include "mcmc_root_finder.h"

#include "../mcmc_models/mcmc_model.h"
#include "../mcmc_models/examples/arma/arma.h"

#include <spin1_api.h>
#include <stdint.h>
//#include <stdbool.h>
#include <debug.h>
#include <data_specification.h>

//#define CALC_TYPE float
//#define ROOT_FAIL -100.000000000000
#define ROOT_FAIL 1 // define it this way here to send "boolean" back

//#ifndef use
//#define use(x) do {} while ((x)!=(x))
//#endif

// Timeout between sending and receiving results in number of timer ticks
//#define TIMEOUT 3

// The model-specific state
//mcmc_state_pointer_t state;

uint32_t *parameter_rec_ptr;

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
	// add a recording region for debug to check the roots calculation

	PARAMETERS
};

struct more_parameters {

    // Acknowledge key for something or other
    uint32_t acknowledge_key;

};

// The general parameters
struct more_parameters more_parameters;

uint32_t ack_key;

// Helper functions for finding the root

// multiply complex by scalar
complex float RCmul( float x, complex float a )
{
	return ( x * crealf( a )) + ( x * cimagf( a )) * I;
}

// edit these constants if necessary for fixed point
#define EPSS 1.0e-7
#define MR 8
#define MT 10
#define MAXIT (MT*MR)
/*
	Given degree m and the m+1 complex coefficients a[0..m] of the polynomial a[i]*x^i and a complex value x,
	this function improves x by Laguerre's method until it converges - within the achievable roundoff limit -
	to a root of the given polynomial. Number of iterations taken is returned as its
*/
void laguerre_poly_root( complex float a[], int m, complex float *x, int *its )
{
	int 				iter,j;
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

		if ( cabsf( b ) <= err ) return;	// we are on the root

		g = d / b;								// the generic case so use Laguerre's formula
		g2 = g * g;
		//h = g2 - RCmul( TWO, f / b );
		h = g2 - RCmul( 2.0f, f / b );
		sq = csqrtf( RCmul( (float) (m-1), RCmul((float) m, h ) - g2 ));
		gp = g + sq;
		gm = g - sq;
		abp = cabsf( gp );
		abm = cabsf( gm );

		if ( abp < abm ) gp = gm;

		dx = FMAX( abp, abm ) > 0.0f ? ( ((float) m ) + 0.0f * I ) / gp :
		RCmul( 1.0f + abx, cosf((float)iter) + sinf((float)iter) * I );
//		dx = FMAX( abp, abm ) > ZERO ? ( ((float) m ) + ZERO * I ) / gp :
//		RCmul( ONE + abx, cosf((float)iter) + sinf((float)iter) * I );

		x1 = *x - dx;

		if ( crealf(*x) == crealf(x1) && cimagf(*x) == cimagf(x1) ) return; // converged

		if ( iter % MT )
			*x = x1;
		else
			*x = *x - RCmul( frac[ iter / MT ], dx ); // occasionally take fractional step to break rare limit cycle
		}

//	printf("Too many iterations in laguerre_poly_root()"); exit(1); // very unusual and only for complex roots

	return;
}

//#undef EPSS
//#undef MR
//#undef MT
//#undef MAXIT

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
	Given the degree m and m+1 complex coefficients a[0..m] of the polynomial a[i]*x^i this function
	successively calls laguerre_poly_root() and finds all m complex roots in roots[1..m]. The bool
	polish should be input as true if polishing is desired (i.e. almost always)
*/
void zroots( complex float a[], int m, complex float roots[], bool polish )
{
	int 				i, its, j, jj;
	complex float 	x, b, c, ad[MAXM];

	for ( j = 0; j <= m; j++ ) ad[j] = a[j]; // copy coeffs for successive deflation

	for ( j = m; j >= 1; j-- ) {	// loop over each root to be found

		//x = ZERO + ZERO * I;  		// start at zero to favour convergence to smallest root
		x = 0.0f + 0.0f * I;  		// start at zero to favour convergence to smallest root

		laguerre_poly_root( ad, j, &x, &its );

		if ( fabs(cimagf(x)) <= 2.0f * EPS * fabs(crealf(x)))
			x = crealf(x) + 0.0f * I; // set imaginary part to zero
//		if ( fabs(cimagf(x)) <= TWO * EPS * fabs(crealf(x)))
//			x = crealf(x) + ZERO * I; // set imaginary part to zero

		roots[j] = x;

		b = ad[j]; // forward deflation

		for ( jj = j-1; jj>=0; jj-- ) {
			c = ad[jj];
			ad[jj] = b;
			b = (x * b) + c;
			}
		}

	if ( polish )
		for ( j = 1; j <= m; j++ )  	// polish roots using undeflated coeffs
			laguerre_poly_root( a, m, &roots[j], &its );

	for ( j = 2; j <= m; j++ ) {		// sort roots by their real parts

		x = roots[j];

		for ( i = j-1; i >= 1; i-- ) {
			if ( crealf(roots[i]) <= crealf(x) ) break;
			roots[i+1] = roots[i];
			}
		roots[i+1] = x;
		}
}

//#undef EPS
//#undef MAXM

uint32_t mcmc_model_get_state_n_bytes() {
    return sizeof(struct mcmc_state);
}

void run(uint unused0, uint unused1) {
    use(unused0);
    use(unused1);

    // for debug writing values
    char buffer[1024];

    uint8_t i, p, q;
    CALC_TYPE state_parameters[PPOLYORDER+QPOLYORDER+2];

    p = PPOLYORDER;
    q = QPOLYORDER;

	// This is probably the point where the executables need to be separated;
	// So if sigma > ZERO then we need to write the parameters to SDRAM?
	// Then the other executable reads these parameters back from SDRAM,
	// sets the data structure up and tests for the root magnitudes, sending
	// back ZERO or ROOT_FAIL... ?
	// I'm assuming here that somehow I can write to SDRAM and then the other
	// executable can read it... there needs to be some form of control
	// whereby this is possible...

//    log_info("ROOT FINDER: in run() function, p=%d", p);

    // Set up the parameters and the state
	uint32_t state_n_bytes = mcmc_model_get_state_n_bytes();

//	log_info("ROOT FINDER: state_n_bytes: %d", state_n_bytes);

	// parameter_rec_ptr should have arrived at this point
//	log_info("ROOT FINDER: parameter_rec_ptr[0]: %d", parameter_rec_ptr[0]);

	// get parameters from sdram
	spin1_memcpy(state_parameters, parameter_rec_ptr[0], state_n_bytes);

//	spin1_callback_off(MCPL_PACKET_RECEIVED); // turn it off?

//	spin1_memcpy(state, model_state_address, state_n_bytes);

//	state_parameters = state->parameters;

//	print_value(state_parameters[1], buffer);
//	log_info("ROOT FINDER: state_parameters[1] = %s", buffer);
//	print_value(state_parameters[10], buffer);
//	log_info("ROOT FINDER: state_parameters[10] = %s", buffer);
//	print_value(state_parameters[18], buffer);
//	log_info("ROOT FINDER: state_parameters[18] = %s", buffer);

//	log_info("ROOT_FINDER: state_parameters[0]: %d", state_parameters[0]);





	// set up data structures for the coefficients of a polynomial
	// characteristic equation for AR and MA model respectively.
	// C99 compile flag required for complex type and variable length arrays
	CALC_TYPE AR_eq[p+1], MA_eq[q+1];
	complex float AR_param[p+1], MA_param[q+1], AR_rt[p+1], MA_rt[q+1];

	// This command reverses the sequence of the first p parameters,
	// get their negative values and adds an element 1 at the end
	// to form the coeffients of AR characteristic equation.
	// The characteristic equation is -a_p*x^p-a_(p-1)*x^(p-1)-...+1=0;

	// we don't need to do this reversal here as the c function assumes the
	// original order is correct!  gah!

	//AR_eq[p] = ONE;  // REAL_CONST( 1.0 );  // ONE;
	AR_eq[0] = ONE;
	for(i=0; i < p; i++) {
		//AR_eq[i] = -state_parameters[p-i-1];
		AR_eq[i+1] = -state_parameters[i];
//		print_value_rf(AR_eq[i], buffer);
//		log_info("check parameter value (AR) = %s", buffer);
//		log_info("check parameter value %d (AR) = %k", i, (accum) AR_eq[i]);
	}

	// This command reverses the sequence from the p+1 to q elements
	// of the parameters, get their negative values and add an element 1
	// at the end to form the coeffients of MA characteristic equation.
	// The characteristic equation is -b_q*x^q-b_(q-1)*x^(q-1)-...+1=0;

	// we don't need to do this reversal here as the c function assumes the
	// original order is correct!  gah!

//	MA_eq[q] = ONE;  // REAL_CONST( 1.0 );  // ONE:
	MA_eq[0] = ONE;
	for(i=0; i < q; i++) {
//		MA_eq[i] = -state_parameters[p+(q-i-1)];
		MA_eq[i+1] = -state_parameters[p+i];
//		print_value_rf(MA_eq[i], buffer);
//		log_info("check parameter value (MA) = %s", buffer);
//		log_info("check parameter value %d (MA) = %k", i, (accum) MA_eq[i]);
	}

	// read parameters into complex vectors so that we can calculate roots - from 0?
	for(i=0; i <= p; i++) {
		// only real part relevant - so imaginary part = 0
		AR_param[i] = (float)AR_eq[i] + ZERO * I;
	}

	for( i = 0; i <= q; i++ ) {
		// only real part relevant - so imaginary part = 0
		MA_param[i] = (float)MA_eq[i] + ZERO * I;
	}

//    log_info("ROOT FINDER: call zroots for AR");

	// this function finds the complex roots in each case
	zroots( AR_param, p, AR_rt, true);

//	log_info("ROOT FINDER: call zroots for MA");

	zroots( MA_param, q, MA_rt, true);

	// debug: record the results of the root finder and send back?
	//        (note: this requires extra setup on the python side as well)
	//recording_record(0, AR_rt, sizeinbytes(AR_rt));
	//recording_record(0, MA_rt, sizeinbytes(MR_rt));

	uint8_t returnval = 0;  // 0.0f;

//	print_value_rf(EPSS, buffer);
//	log_info("EPSS = %s", buffer);
//	print_value_rf(EPS, buffer);
//	log_info("EPS = %s", buffer);

	// test for root magnitude <= 1 and if so return a fail result  // n+1 coefficients, n roots
	for(i=1; i <= p; i++) { // 0 or 1 for start point?
		if( cabsf(AR_rt[i]) <= ONE ) returnval = ROOT_FAIL; // return ROOT_FAIL;  // REAL_CONST( ROOT_FAIL );
	}
	//		print_value_rf(cabsf(AR_rt[i]), buffer);
	//		log_info("root magnitude (AR) = %s", buffer);

	for(i=1; i <= q; i++)  // 0 or 1 for start point?
		if( cabsf(MA_rt[i]) <= ONE ) returnval = ROOT_FAIL;  // REAL_CONST( ROOT_FAIL );

	// At this point send returnval to normal vertex
//	log_info("ack_key = 0x%08x", ack_key);
//	print_value_rf(returnval, buffer);
//	log_info("check after magnitude test, returnval = %s", buffer);
//	log_info("check returnval using accum = %k", returnval);

	// wait until packet is acknowledged/sent...
	while (!spin1_send_mc_packet(ack_key, returnval, WITH_PAYLOAD)) {
		spin1_delay_us(1);
	}

// if all conditions have been passed then return a pass result
	// return ZERO;  // REAL_CONST( 0.0 ); here send ZERO to sdram
//	log_info("ROOT FINDER: returnval sent to ARMA");

	//log_info("ROOT FINDER: at end of run(), returnval = 0x%08x", returnval);

	// DEBUG: exit here to prevent too much writing to iobuf - we just
	//        want to test whether the values coming back out are sensible
//    spin1_exit(0);
}

//void multicast_callback(uint key, uint payload) {
//	use(key);
//	parameter_rec_ptr[0] = payload;
////    log_info("ROOT FINDER: multicast_callback");
////    spin1_callback_off(TIMER_TICK);
////    spin1_schedule_callback(run, 0, 0, 2);
//}

//	// this needs to be edited - it's not the data sequence that comes in here
//	// but the params
//
//
////    uint sequence = key & parameters.sequence_mask;
////    if (sequence == next_sequence) {
////        data_receive_ptr[0] = payload;
////        data_receive_ptr++;
////        last_sequence = sequence;
////        next_sequence = (sequence + 1) & parameters.sequence_mask;
////    }
//}
//
//void timer_callback(uint time, uint unused) {
//	// again here, not data sequence, but params
//
////    spin1_delay_us(parameters.timer >> 1);
////    use(time);
////    use(unused);
////    spin1_send_mc_packet(
////        parameters.acknowledge_key, last_sequence, WITH_PAYLOAD);
//}
//
//void empty_multicast_callback(uint key, uint payload) {
//    use(key);
//    use(payload);
//}

//void trigger_run(uint unused0, uint unused1) {
void trigger_run(uint key, uint payload) {
//    use(unused0);
//    use(unused1);
	use(key);
	parameter_rec_ptr[0] = payload;
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
	// Get the ack_key from more_parameters
	address_t data_address = data_specification_get_data_address();
	address_t more_parameters_address = data_specification_get_region(
	        PARAMETERS, data_address);
	struct more_parameters *more_sdram_params =
			(struct more_parameters *) more_parameters_address;
	spin1_memcpy(&more_parameters, more_sdram_params,
			sizeof(struct more_parameters));

	ack_key = more_parameters.acknowledge_key;

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
