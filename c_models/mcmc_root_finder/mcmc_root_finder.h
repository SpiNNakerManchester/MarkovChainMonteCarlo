// compile with gcc -std=c99 -O3 -fcx-limited-range root_finder.c -l m -o rf

#include "mcmc_spinn_real_type.h"

// some utility functions
static float maxarg1,maxarg2;
#define FMAX( a, b ) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ? (maxarg1) : (maxarg2))

// multiply complex by scalar
complex float RCmul( float x, complex float a );

/*
	Given degree m and the m+1 complex coefficients a[0..m] of the polynomial
	a[i]*x^i and a complex value x, this function improves x by Laguerre's
	method until it converges - within the achievable roundoff limit - to a
	root of the given polynomial. Number of iterations taken is returned as its
*/
void laguerre_poly_root( complex float a[], int m, complex float *x, int *its );

/*
	Given the degree m and m+1 complex coefficients a[0..m] of the polynomial
	a[i]*x^i this function successively calls laguerre_poly_root() and finds
	all m complex roots in roots[1..m]. The bool polish should be input as
	true if polishing is desired (i.e. almost always)
*/
void zroots( complex float a[], int m, complex float roots[], bool polish );
