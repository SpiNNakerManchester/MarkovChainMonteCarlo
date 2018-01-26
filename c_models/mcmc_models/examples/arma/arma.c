#include "../../mcmc_model.h"
#include "arma.h"

#define ROOT_FAIL CONCAT(-100.000000, SUFFIX)		// how bad is a root failure
#define REAL float

uint32_t mcmc_model_get_params_n_bytes() {
    return sizeof(struct mcmc_params);
}

uint32_t mcmc_model_get_state_n_bytes() {
    return sizeof(struct mcmc_state);
}

/*
 ARMA wind log likelihood
  - data, size n_pts
  - use state->order_p and state->order_q which are the orders of an ARMA(p,q) model,
    which has the form z(t)=a_1*z(t-1)+...+a_p*z(t-p)-b_1*e(t-1)-...-b_q*e(t-q)+mu
  - params is [a1,...,ap,b1,...,bq,mu,sigma]
 */
CALC_TYPE mcmc_model_likelihood(
        CALC_TYPE *data, uint32_t n_pts, mcmc_params_pointer_t params,
		mcmc_state_pointer_t state) {

	// read in AR and MA parameter dimensions
	uint8_t p = state->order_p;
	uint8_t q = state->order_q;
	uint8_t i, j;
	uint32_t N = n_pts;

	uint32_t error_length = N + q;  // check this and all indexing below!
	CALC_TYPE err[error_length];

	// get parameters array from struct
	CALC_TYPE *parameters = params->parameters;

	// temp storage for dot products in later loop
	CALC_TYPE tempdotp = ZERO; // REAL_CONST( 0.0 );
	CALC_TYPE tempdotq = ZERO; // REAL_CONST( 0.0 );

	// mu is the penultimate element of the parameters
	// sigma is the last element of the parameters
	CALC_TYPE mu = parameters[p+q];
	CALC_TYPE sigma = parameters[p+q+1];

/*
err is all zeros with size of (N+q)*1
err=zeros(N+q,1); % this vector will become non-zero from the p+q+1
*/

	for(i=0; i < error_length; i++)
		err[i] = ZERO;  // REAL_CONST( 0.0 );

/*
	Y=zeros(N,1);% this is the predicted output
*/
	REAL Y[N]; // C99 compile flag required

	// compiler issue, could add mu here rather than in the main loop

/*
MATLAB code that works:

for i=p+1:N
    Y(i)=parameters(1:p)*data(i-1:-1:i-p)-parameters(p+1:p+q)*err(i+q-1:-1:i)+mu;
    err(i+q)=data(i)-Y(i);
end

// parameters(1:p) refers to the elements from the 1st to the pth.
// data(i-1:-1:i-p) refers to the elements from the i-1 backward one by one to i-p
// parameters(1:p) refers to the elements from the (p+1)th to the (p+q)th.
// err(i+q-1:-1:i) refers to the elements from i+q-1 backward one by one to i
*/

	// OK this is the part that is slightly complex to understand - very careful about indexing here
	for(i=p; i < N; i++) {
		// loop over p for parameters * data dot product
		for (j=0; j < p; j++) {
			tempdotp = parameters[j] * data[(i-1)-j]; // check this
		}
		// loop over q for parameters * err dot product
		for (j=0; j < q; j++) {
			tempdotq = parameters[p+j] * err[(i+q-1)-j]; // check this
		}

		// Add two results together plus mean
		Y[i] = tempdotp + tempdotq + mu;
		// Error is data minus predicted data
		err[i+q] = data[i] - Y[i];

		// Reset dot products
		tempdotp = ZERO;
		tempdotq = ZERO;
	}

/*

% The first part of log of likelihood equals the negative sum of square error and then divided
% by double of square sigma. The second part is half the number of (N-p) times the log of square sigma.
% The substraction between the first part and the second part equals the lglikelihood.

lglikelihood=-sum(err(p+q+1:end).^2)/(2*sigma^2)-0.5*(N-p)*log(sigma^2);

*/
	CALC_TYPE sum = ZERO; // REAL_CONST( 0.0 );
	CALC_TYPE temp;
	CALC_TYPE divisor = ONE / (TWO*SQR(sigma) - HALF*(N-p)*LN(SQR(sigma)));
	// Check value of divisor here (similar to in transition_jump??

	//float divisor = 1.0f / (( 2.0f * sigma * sigma ) - 0.5f * (N-p) * logf( sigma * sigma )); // to avoid divide, create one floating point inverse & multiply later
	for(i=p+q; i < N+q; i++) { // check
		temp = err[i];
		sum += temp * temp;
	}

	return -sum * (CALC_TYPE) divisor;  // possibly add some error checking in case divisor is stupid value
}

/*
 prior probability for the parameters, using state->order_p and state->order_q

 prior condition is that each root is outside the unit circle
 */
CALC_TYPE mcmc_model_prior_prob(
        mcmc_params_pointer_t params, mcmc_state_pointer_t state) {
	// read in AR and MA parameter dimensions
	uint8_t p = state->order_p;
	uint8_t q = state->order_q;
//	uint8_t i;

//	REAL sigma = params->sigma; // could plausibly use this rather than a parameters array?
	CALC_TYPE *parameters = params->parameters;
	CALC_TYPE sigma = parameters[p+q+1];  // last entry in vector
	if( sigma <= ZERO ) return ROOT_FAIL;  // first fail condition can provide early exit

//	// This is probably the point where the executables need to be separated;
//	// So if sigma > ZERO then we need to write the parameters to SDRAM?
//	// Then the other executable reads these parameters back from SDRAM,
//	// sets the data structure up and tests for the root magnitudes, sending
//	// back ZERO or ROOT_FAIL... ?
//	// I'm assuming here that somehow I can write to SDRAM and then the other
//	// executable can read it... there needs to be some form of control
//	// whereby this is possible...
//
//	// set up data structures for the coefficients of a polynomial
//	// characteristic equation for AR and MA model respectively.
//	// C99 compile flag required for complex type and variable length arrays
//	CALC_TYPE AR_eq[p+1], MA_eq[q+1];
//	complex float AR_param[p+1], MA_param[q+1], AR_rt[p+1], MA_rt[q+1];
//
//	// This command reverses the sequence of the first p parameters,
//	// get their negative values and adds an element 1 at the end
//	// to form the coeffients of AR characteristic equation.
//	// The characteristic equation is -a_p*x^p-a_(p-1)*x^(p-1)-...+1=0;
//	AR_eq[p] = ONE; // REAL_CONST( 1.0 );
//	for(i=0; i < p; i++) {
//		AR_eq[i] = -parameters[p-i-1];
//	}
//
//	// This command reverses the sequence from the p+1 to q elements
//	// of the parameters, get their negative values and add an element 1
//	// at the end to form the coeffients of MA characteristic equation.
//	// The characteristic equation is -b_q*x^q-b_(q-1)*x^(q-1)-...+1=0;
//	MA_eq[q] = ONE; // REAL_CONST( 1.0 );
//	for(i=0; i < q; i++) {
//		MA_eq[i] = -parameters[p+(q-i-1)];
//	}
//
//	// read parameters into complex vectors so that we can calculate roots - from 0?
//	for(i=0; i <= p; i++) {
//		// only real part relevant - so imaginary part = 0
//		AR_param[i] = (float)AR_eq[i] + 0.0f * I;
//	}
//
//	for( i = 0; i <= q; i++ ) {
//		// only real part relevant - so imaginary part = 0
//		MA_param[i] = (float)MA_eq[i] + 0.0f * I;
//	}
//
//	// this function finds the complex roots in each case
//	zroots( AR_param, p, AR_rt, true);
//	zroots( MA_param, q, MA_rt, true);
//
//	// test for root magnitude <= 1 and if so return a fail result  // n+1 coefficients, n roots
//	for(i=1; i <= p; i++)  // 0 or 1 for start point?
//		if( cabsf(AR_rt[i]) <= 1.0f ) return ROOT_FAIL;  // REAL_CONST( ROOT_FAIL );
//
//	for(i=1; i <= q; i++)  // 0 or 1 for start point?
//		if( cabsf(MA_rt[i]) <= 1.0f ) return ROOT_FAIL;  // REAL_CONST( ROOT_FAIL );
//
//// if all conditions have been passed then return a pass result
//	return ZERO;  // REAL_CONST( 0.0 );
}


/*
 generate jumps from MH transition distribution which is uncorrelated bivariate
 t with degrees_freedom degrees of freedom
 */
void mcmc_model_transition_jump(
        mcmc_params_pointer_t params, mcmc_state_pointer_t state,
        mcmc_state_pointer_t new_state) {
	// loop over parameters and apply relevant jump_scale
	// - it'll look something like this...
	CALC_TYPE *parameters = params->parameters;
	uint32_t p = state->order_p;
	uint32_t q = state->order_q;
	unsigned int i;
	for (i=0; i < p; i++) {
		parameters[i] += t_deviate() * params->p_jump_scale;
	}
	for (i=p; i < p+q; i++) {
		parameters[i] += t_deviate() * params->q_jump_scale;
	}

//    new_state->order_p = state->order_p + (t_deviate() * params->alpha_jump_scale);
//    new_state->order_q = state->order_q + (t_deviate() * params->beta_jump_scale);
}

// Root finder - it may be the case that the ITCM limit forces us to
//               move everything below into a separate executable

// multiply complex by scalar
//complex float RCmul(float x, complex float a)
//{
//	return (x*crealf(a)) + (x*cimagf(a))*I;
//}
//
//
//#define EPSS 1.0e-7
//#define MR 8
//#define MT 10
//#define MAXIT (MT*MR)
///*
//	Given degree m and the m+1 complex coefficients a[0..m] of the polynomial
//	a[i]*x^i and a complex value x, this function improves x by Laguerre's
//	method until it converges - within the achievable roundoff limit -to a root
//	of the given polynomial. Number of iterations taken is returned as its.
//*/
//void laguerre_poly_root(
//		complex float a[], int m, complex float *x, int *its)
//{
//	int iter,j;
//	float abx, abp, abm, err;
//	complex float dx, x1, b, d, f, g, h, sq, gp, gm, g2;
//	static float frac[MR+1] = { 0.0, 0.5, 0.25, 0.75, 0.13, 0.38, 0.62, 0.88, 1.0 };
//
//	for (iter=1; iter <= MAXIT; iter++) {  // loop over iterations up to maximum
//		*its = iter;
//		b = a[m];
//		err = cabsf( b );
//		d = f = ZERO + ZERO * I; // 0.0f + 0.0f * I;
//		abx = cabsf( *x );
//		for (j=m-1; j>=0; j--) {	// efficient computation of polynomial and first two derivatives
//			f = ( *x * f ) + d;
//			d = ( *x * d ) + b;
//			b = ( *x * b ) + a[j];
//			err = cabsf( b ) + abx * err;
//		}
//
//		err *= EPSS;	// estimate of roundoff error in evaluating polynomial
//
//		if (cabsf(b) <= err) return;	// we are on the root
//
//		g = d / b;								// the generic case so use Laguerre's formula
//		g2 = SQR(g);  // g * g;
//		h = g2 - RCmul( TWO, f / b );
//		sq = csqrtf(RCmul((float) (m-1), RCmul((float) m, h) - g2));
//		gp = g + sq;
//		gm = g - sq;
//		abp = cabsf( gp );
//		abm = cabsf( gm );
//
//		if (abp < abm) gp = gm;
//
////		dx = FMAX( abp, abm ) > 0.0f ? (((float) m ) + 0.0f * I) / gp :
//		dx = FMAX( abp, abm ) > ZERO ? (((float) m ) + ZERO * I) / gp :
////		RCmul(1.0f + abx, cosf((float)iter) + sinf((float)iter) * I);
//		RCmul(ONE + abx, cosf((float)iter) + sinf((float)iter) * I);
//
//		x1 = *x - dx;
//
//		if (crealf(*x) == crealf(x1) && cimagf(*x) == cimagf(x1)) return; // converged
//
//		if (iter % MT)
//			*x = x1;
//		else
//			*x = *x - RCmul(frac[iter/MT], dx); // occasionally take fractional step to break rare limit cycle
//	}
//
//	// printf("Too many iterations in laguerre_poly_root()"); exit(1); // very unusual and only for complex roots
//
//	return;
//}
//
//#undef EPSS
//#undef MR
//#undef MT
//#undef MAXIT
//
///* from NR in C errata
//*** 40,42 ****
//  		dx=((FMAX(abp,abm) > 0.0 ? Cdiv(Complex((float) m,0.0),gp)
// 			: RCmul(exp(log(1+abx)),Complex(cos((float)iter),sin((float)iter)))));
//  		x1=Csub(*x,dx);
//--- 40,42 ----
//  		dx=((FMAX(abp,abm) > 0.0 ? Cdiv(Complex((float) m,0.0),gp)
// 			: RCmul(1+abx,Complex(cos((float)iter),sin((float)iter)))));
//  		x1=Csub(*x,dx);
//
//*/
//
//#define EPS 2.0e-6
//#define MAXM 100
///*
//	Given the degree m and m+1 complex coefficients a[0..m] of the polynomial a[i]*x^i this function
//	successively calls laguerre_poly_root() and finds all m complex roots in roots[1..m]. The bool
//	polish should be input as true if polishing is desired (i.e. almost always)
//*/
//void zroots(
//		complex float a[], int m, complex float roots[], bool polish)
//{
//	int i, its, j, jj;
//	complex float x, b, c, ad[MAXM];
//
//	for (j=0; j<=m; j++) ad[j] = a[j]; // copy coeffs for successive deflation
//
//	for (j=m; j>=1; j--) {	// loop over each root to be found
//
//		x = ZERO + ZERO * I;  		// start at zero to favour convergence to smallest root
//
//		laguerre_poly_root(ad, j, &x, &its);
//
//		if (fabs(cimagf(x)) <= TWO * EPS * fabs(crealf(x)))
//			x = crealf(x) + ZERO * I; // set imaginary part to zero
//
//		roots[j] = x;
//
//		b = ad[j]; // forward deflation
//
//		for (jj=j-1; jj>=0; jj--) {
//			c = ad[jj];
//			ad[jj] = b;
//			b = (x * b) + c;
//		}
//	}
//
//	if (polish)
//		for (j=1; j<=m; j++)  	// polish roots using undeflated coeffs
//			laguerre_poly_root(a, m, &roots[j], &its);
//
//	for (j=2; j <= m; j++) {		// sort roots by their real parts
//
//		x = roots[j];
//
//		for (i=j-1; i>=1; i--) {
//			if (crealf(roots[i]) <= crealf(x)) break;
//			roots[i+1] = roots[i];
//		}
//		roots[i+1] = x;
//	}
//}
//
//#undef EPS
//#undef MAXM
//

