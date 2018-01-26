

#include "spinn_real_type.h"  // REAL -> float, double, accum or whatever
#include "root_finder.h"

#define ROOT_FAIL -100.0		// how bad is a root failure


/*

	Log prior function for ARIMA wind modelling

	Penalises infeasible ARMA parameter sets

	Brought together from work in C by MH and work in MATLAB by JY

 This function defines the log of the priors belief over the parameters;
 z(t)=a1*z(t-1)+a2*z(t-2)+a3*z(t-3)+...a_p*z(t-p)+e(t)-b1*e(t-1)-b2*e(t-2)-...-b_q*e(t-q)+mu

 order=[p,q] defines the order of ARMA(p,q);
 parameters are the variables involved in the structure.
 The number of parameters equals p+q+2;
 parameters=[a1,...a_p,b1,...b_q,mu,sigma]

*/
REAL logprior( uint8_t order[2], REAL parameters )
{
	void zroots( complex float a[], int m, complex float roots[], bool polish );

	uint8_t p = order[0], q = order[1], i, index;	// read in AR and MA parameter dimensions

	REAL sigma = parameters[ p+q+1 ]; // last entry in vector
	if( sigma <= REAL_CONST( 0.0 ) ) return REAL_CONST( ROOT_FAIL );  // first fail condition can provide early exit


// set up data structures for the coefficients of a polynomial characteristic equation for AR and MA model respectively.
// C99 compile flag required for complex type and variable length arrays
	REAL AR_eq[p+1], MA_eq[q+1];
	complex float AR_param[p+1], MA_param[q+1], AR_rt[p+1], MA_rt[q+1];


// This command reverses the sequence of the first p parameters, get their negative values and add an element 1 in the end
// to form the coeffients of AR characteristic equation.
// The characteristic equation is -a_p*x^p-a_(p-1)*x^(p-1)-...+1=0;
//	index = p;
// check all loop indexing! probably wrong ATM
	AR_eq[p] = REAL_CONST( 1.0 );
	for( i = 0; i < p; i++ )
		AR_eq[i] = -parameters[ p - i - 1 ];

// This command reverses the sequence from the p+1 to q elements of the parameters, get their negative values and add an
// element 1 in the end to form the coeffients of MA characteristic equation.
// The characteristic equation is -b_q*x^q-b_(q-1)*x^(q-1)-...+1=0;
// check all loop indexing! probably wrong ATM
	MA_eq[q] = REAL_CONST( 1.0 );
	for( i = 0; i < q; i++ )
		MA_eq[i] = -parameters[ p+ (q - i - 1) ];


// read parameters into complex vectors so that we can calculate roots - from 0?
	for( i = 0; i <= p; i++ )
		AR_param[i] = (float)AR_eq[i] + 0.0f * I;  // only real part relevant - so imaginary part = 0

	for( i = 0; i <= q; i++ )
		MA_param[i] = (float)MA_eq[i] + 0.0f * I;  // only real part relevant - so imaginary part = 0


// this function finds the complex roots in each case
	zroots( AR_param, p, AR_rt, true ); // check length values (2nd argument)

	zroots( MA_param, q, MA_rt, true ); // " "


// test for root magnitude > 1 and if so return a fail result  // n+1 coefficients, n roots
	for( i = 1; i <= p; i++ )  // 0 or 1 for start point?
		if( cabsf( AR_rt[i] ) <= 1.0f ) return REAL_CONST( ROOT_FAIL ); // less than ? (failure if any root INSIDE unit circle...)

	for( i = 1; i <= q; i++ )  // 0 or 1 for start point?
		if( cabsf( MA_rt[i] ) <= 1.0f ) return REAL_CONST( ROOT_FAIL ); // sim. to above?


// if all conditions have been passed then return a pass result
	return REAL_CONST( 0.0 );
}
