
// compile with gcc -std=c99 -O3 -fcx-limited-range root_finder.c -l m -o rf

#include "mcmc_root_finder.h"

#include <spin1_api.h>
#include <stdint.h>
#include <stdbool.h>
#include <debug.h>
#include <data_specification.h>

#define CALC_TYPE float
#define ROOT_FAIL -100.0

#ifndef use
#define use(x) do {} while ((x)!=(x))
#endif

// Timeout between sending and receiving results in number of timer ticks
#define TIMEOUT 3

// The parameters to be read from memory
enum params {
    DATA_SIZE = 0,
    N_CHIPS,
    KEY,
    WINDOW_SIZE,
    SEQUENCE_MASK,
    TIMER,
    DATA
};

// An array of sequence numbers received on each core
// Note that this potentially will be in SDRAM with enough cores
uint *sequence_received;

// The list of keys for each of the cores - ordered for quick searching
uint *chip_keys;

// The size of the remaining data to be sent
uint data_size;

// The number of chips the data is to be sent to
uint n_chips;

// The base key to send the data with
uint key;

// The window size for sending the data
uint window_size;

// The mask which indicates the sequence number
uint sequence_mask;

// Pointer to the start of the data still to be sent and confirmed
uint *data;

// The timer tick at which the send will have timed out
uint send_timeout;

// The next sequence number to be sent
uint next_sequence = 0;

// The end sequence ignoring the wrap around of sequences
uint next_end_sequence_unwrapped = 0xFFFFFFFF;

// The end sequence with wrap around of sequences
uint next_end_sequence = 0;

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
		h = g2 - RCmul( 2.0f, f / b );
		sq = csqrtf( RCmul( (float) (m-1), RCmul((float) m, h ) - g2 ));
		gp = g + sq;
		gm = g - sq;
		abp = cabsf( gp );
		abm = cabsf( gm );

		if ( abp < abm ) gp = gm;

		dx = FMAX( abp, abm ) > 0.0f ? ( ((float) m ) + 0.0f * I ) / gp :
		RCmul( 1.0f + abx, cosf((float)iter) + sinf((float)iter) * I );

		x1 = *x - dx;

		if ( crealf(*x) == crealf(x1) && cimagf(*x) == cimagf(x1) ) return; // converged

		if ( iter % MT )
			*x = x1;
		else
			*x = *x - RCmul( frac[ iter / MT ], dx ); // occasionally take fractional step to break rare limit cycle
		}

	printf("Too many iterations in laguerre_poly_root()"); exit(1); // very unusual and only for complex roots

	return;
}

#undef EPSS
#undef MR
#undef MT
#undef MAXIT

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

		x = 0.0f + 0.0f * I;  		// start at zero to favour convergence to smallest root

		laguerre_poly_root( ad, j, &x, &its );

		if ( fabs(cimagf(x)) <= 2.0f * EPS * fabs(crealf(x)))
			x = crealf(x) + 0.0f * I; // set imaginary part to zero

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

#undef EPS
#undef MAXM

void run(uint unused0, uint unused1) {
    use(unused0);
    use(unused1);

    uint8_t i, p, q;
    CALC_TYPE *parameters;

    // Set up the data somewhere in here


	// This is probably the point where the executables need to be separated;
	// So if sigma > ZERO then we need to write the parameters to SDRAM?
	// Then the other executable reads these parameters back from SDRAM,
	// sets the data structure up and tests for the root magnitudes, sending
	// back ZERO or ROOT_FAIL... ?
	// I'm assuming here that somehow I can write to SDRAM and then the other
	// executable can read it... there needs to be some form of control
	// whereby this is possible...

    // Need to get p and q from somewhere...

	// set up data structures for the coefficients of a polynomial
	// characteristic equation for AR and MA model respectively.
	// C99 compile flag required for complex type and variable length arrays
	CALC_TYPE AR_eq[p+1], MA_eq[q+1];
	complex float AR_param[p+1], MA_param[q+1], AR_rt[p+1], MA_rt[q+1];

	// This command reverses the sequence of the first p parameters,
	// get their negative values and adds an element 1 at the end
	// to form the coeffients of AR characteristic equation.
	// The characteristic equation is -a_p*x^p-a_(p-1)*x^(p-1)-...+1=0;
	AR_eq[p] = REAL_CONST( 1.0 );  // ONE;
	for(i=0; i < p; i++) {
		AR_eq[i] = -parameters[p-i-1];
	}

	// This command reverses the sequence from the p+1 to q elements
	// of the parameters, get their negative values and add an element 1
	// at the end to form the coeffients of MA characteristic equation.
	// The characteristic equation is -b_q*x^q-b_(q-1)*x^(q-1)-...+1=0;
	MA_eq[q] = REAL_CONST( 1.0 );  // ONE:
	for(i=0; i < q; i++) {
		MA_eq[i] = -parameters[p+(q-i-1)];
	}

	// read parameters into complex vectors so that we can calculate roots - from 0?
	for(i=0; i <= p; i++) {
		// only real part relevant - so imaginary part = 0
		AR_param[i] = (float)AR_eq[i] + 0.0f * I;
	}

	for( i = 0; i <= q; i++ ) {
		// only real part relevant - so imaginary part = 0
		MA_param[i] = (float)MA_eq[i] + 0.0f * I;
	}

	// this function finds the complex roots in each case
	zroots( AR_param, p, AR_rt, true);
	zroots( MA_param, q, MA_rt, true);

	CALC_TYPE returnval = 0.0f;  // ZERO;

	// test for root magnitude <= 1 and if so return a fail result  // n+1 coefficients, n roots
	for(i=1; i <= p; i++)  // 0 or 1 for start point?
		if( cabsf(AR_rt[i]) <= 1.0f ) returnval = ROOT_FAIL; // return ROOT_FAIL;  // REAL_CONST( ROOT_FAIL );

	for(i=1; i <= q; i++)  // 0 or 1 for start point?
		if( cabsf(MA_rt[i]) <= 1.0f ) returnval = ROOT_FAIL;  // REAL_CONST( ROOT_FAIL );

	// At this point send returnval to SDRAM

// if all conditions have been passed then return a pass result
	// return ZERO;  // REAL_CONST( 0.0 ); here send ZERO to sdram
}

void send_callback(uint send_time, uint unused) {
    use(unused);

    // If the data has all been sent, send a start message and quit
    if (data_size == 0) {
        log_info("All data has been sent and confirmed");
        while (!spin1_send_mc_packet(key, 0, NO_PAYLOAD)) {
            spin1_delay_us(1);
        }
        spin1_exit(0);
        return;
    }

    // Send the next packets
    uint packets_to_send = window_size;
    if (window_size > data_size) {
        packets_to_send = data_size;
    }
    for (uint i = 0; i < packets_to_send; i++) {
        uint sequence = (next_sequence + i) & sequence_mask;
        while (!spin1_send_mc_packet(
                key + sequence, data[i], WITH_PAYLOAD)) {
            spin1_delay_us(1);
        }
    }
    next_end_sequence_unwrapped = next_sequence + packets_to_send - 1;
    next_end_sequence = next_end_sequence_unwrapped & sequence_mask;
    send_timeout = send_time + TIMEOUT;
}


void timer_callback(uint time, uint unused) {
    use(unused);

    // Only do anything if there is sequence to check
    if (next_end_sequence_unwrapped == 0xFFFFFFFF) {
        return;
    }

    // Determine if the timeout would expire
    uint timed_out = 0;
    if (time >= send_timeout) {
        timed_out = 1;
    }

    // Find the smallest sequence number received
    uint smallest_sequence = 0xFFFFFFFF;
    for (uint i = 0; i < n_chips; i++) {
        if (sequence_received[i] < smallest_sequence) {
            smallest_sequence = sequence_received[i];

            // If one of the cores hasn't yet got the latest sequence
            // and we haven't timed out, we may as well give up
            if ((smallest_sequence < next_end_sequence_unwrapped) &&
                    !timed_out) {
                break;
            }
        }
    }

    // If all chips have the data, or we have timed out, start sending again
    if ((smallest_sequence == next_end_sequence_unwrapped) || timed_out) {

        // Adjust the sequence to the next sequence to send
        if (smallest_sequence >= next_sequence &&
                smallest_sequence <= next_end_sequence_unwrapped) {

            // Work out how many words have been sent
            uint words_sent = (smallest_sequence - next_sequence) + 1;

            // Adjust the pointers to the next bit of data to send
            data = &(data[words_sent]);
            data_size -= words_sent;

            next_sequence = (smallest_sequence + 1) & sequence_mask;
        }

        // Avoid sending more data in case the send takes too long
        next_end_sequence_unwrapped = 0xFFFFFFFF;

        // Send the data
        spin1_schedule_callback(send_callback, time, 0, 1);
    }
}

void multicast_callback(uint key, uint payload) {

    // Find the core that this packet is from using binary search
    uint imin = 0;
    uint imax = n_chips;
    while (imin < imax) {

        uint imid = (imax + imin) >> 1;

        // If the key is found, update the sequence with the payload
        if (chip_keys[imid] == key) {

            uint sequence = payload;
            uint current_sequence = sequence_received[imid];

            // If the range is wrapped, adjust the sequence to ignore the wrap
            if (next_end_sequence < next_sequence) {
                if (sequence < next_sequence) {
                    sequence += sequence_mask + 1;
                }
            }


            // Only update if the sequence is in the expected range
            if (sequence >= next_sequence &&
                    sequence <= next_end_sequence_unwrapped) {

                // If the current sequence is out of range, update
                if (current_sequence < next_sequence ||
                        current_sequence > next_end_sequence_unwrapped ||
                        sequence > current_sequence) {
                    sequence_received[imid] = sequence;
                }
            }
            break;
        }

        if (chip_keys[imid] < key) {
            imin = imid + 1;
        } else {
            imax = imid;
        }
    }
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

void c_main() {

    address_t data_address = data_specification_get_data_address();
    address_t params = data_specification_get_region(0, data_address);

    // Get the size of the data in words
    data_size = params[DATA_SIZE];
    log_info("Data size = %d", data_size);

    // Get a count of the chips to load on to
    n_chips = params[N_CHIPS];
    log_info("N chips = %d", n_chips);

    // Get the key to send the data with
    key = params[KEY];
    log_info("Key = 0x%08x", key);

    // Get the number of packets to be sent at the same time
    window_size = params[WINDOW_SIZE];
    log_info("Window size = %d", window_size);

    // Get the total number of window spaces available
    sequence_mask = params[SEQUENCE_MASK];
    log_info("Sequence mask = 0x%08x", sequence_mask);

    // Get the timer tick
    uint timer = params[TIMER];
    log_info("Timer = %d", timer);

    // Get a pointer to the data - not worth copying at present
    data = (uint *) &(params[DATA]);

    // Get a pointer to the keys for each chip which will verify
    // the reception of data
    uint *chip_data = (uint *) &(data[data_size]);

    // Try to allocate an array of keys in DTCM
    chip_keys = (uint *) spin1_malloc(n_chips * sizeof(uint));
    if (chip_keys == NULL) {
        log_warning("Could not allocate chip keys in DTCM - using SDRAM");
        chip_keys = chip_data;
    } else {
        spin1_memcpy(chip_keys, chip_data, n_chips * sizeof(uint));
    }

    // Try to allocate an array of last sequences in DTCM
    sequence_received = (uint *) spin1_malloc(n_chips * sizeof(uint));
    if (sequence_received == NULL) {
        log_warning("Could not allocate sequences in DTCM - using SDRAM");
        sequence_received = (uint *) sark_xalloc(
                sv->sdram_heap, n_chips * sizeof(uint), 0, ALLOC_LOCK);
        if (sequence_received == NULL) {
            log_error("Could not allocate sequences in SDRAM");
            rt_error(RTE_SWERR);
        }
    }
    for (uint i = 0; i < n_chips; i++) {
        sequence_received[i] = sequence_mask;
    }
    send_timeout = TIMEOUT;

    // Set up the timer
    spin1_set_timer_tick(timer);

    // Register for the start message
    spin1_callback_on(MC_PACKET_RECEIVED, trigger_run, -1);

    // Setup callback on multicast packet with payload, for acknowledge packets
    spin1_callback_on(MCPL_PACKET_RECEIVED, multicast_callback, -1);
    spin1_callback_on(MC_PACKET_RECEIVED, empty_multicast_callback, -1);

    // Setup callback on timer to timeout sent packets and send next packets
    spin1_callback_on(TIMER_TICK, timer_callback, 2);

    // Schedule the first call of the send callback
    spin1_schedule_callback(send_callback, 0, 0, 1);

    // Start in sync with all the cores, to ensure they are ready
    // to receive packets
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
