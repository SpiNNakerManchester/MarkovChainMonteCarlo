//#include <complex.h>
//#include <stdbool.h>
//#include <stdtype.h> // not sure what this is or is supposed to be
#include "../../mcmc_model.h"

#define PPOLYORDER 9
#define QPOLYORDER 9

struct mcmc_params {

    // scaling of transition distribution for MH jumps
    CALC_TYPE p_jump_scale[PPOLYORDER];
    CALC_TYPE q_jump_scale[QPOLYORDER];
    CALC_TYPE mu_jump_scale;
    CALC_TYPE sigma_jump_scale;

    // polynomial orders for each polynomial
//    uint8_t order_p;
//    uint8_t order_q;

};

struct mcmc_state {
    // parameters has size order_p (p poly) + order_q (q poly) + 2 (mu and sigma)
    CALC_TYPE parameters[PPOLYORDER+QPOLYORDER+2];
};

//enum regions {
//    RECORDING,
//    PARAMETERS,
//    MODEL_PARAMETERS,
//    MODEL_STATE
//};

// Put the root finder stuff here

// some utility functions
// float MAX()
//static float maxarg1,maxarg2;
//#define FMAX( a, b ) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ? (maxarg1) : (maxarg2))
//
//// multiply complex by scalar
//complex float RCmul(float x, complex float a);
//
///*
//	Given degree m and the m+1 complex coefficients a[0..m] of the polynomial a[i]*x^i and a complex value x,
//	this function improves x by Laguerre's method until it converges - within the achievable roundoff limit -
//	to a root of the given polynomial. Number of iterations taken is returned as its
//*/
//void laguerre_poly_root(
//		complex float a[], int m, complex float *x, int *its);
//
///*
//	Given the degree m and m+1 complex coefficients a[0..m] of the polynomial a[i]*x^i this function
//	successively calls laguerre_poly_root() and finds all m complex roots in roots[1..m]. The bool
//	polish should be input as true if polishing is desired (i.e. almost always)
//*/
//void zroots(
//		complex float a[], int m, complex float roots[], bool polish);
