
// compile with gcc -std=c99 -O3 -fcx-limited-range root_finder.c -l m -o rf

#include "root_finder.h"


// multiply complex by scalar
complex float RCmul( float x, complex float a )
{
	return ( x * crealf( a )) + ( x * cimagf( a )) * I;
}


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


#define ORDER 18

/*
	test program - set up some polynomial coefficients and find their complex roots
	seems to works up to high orders but needs testing on some real ARMA coefficients
*/
//int test_main( void )
int main( void )
{
	int  i, m = ORDER;
	bool pos = true;  // alternator for coeff signs
	complex float a[ORDER+1], roots[ORDER+1];

	for( i = 0; i <= ORDER; i++ ) {  // set up some polynomial coefficients with alternating signs
		if( pos ) {
			a[i] = ( (float)(ORDER - i) + 0.5f ) + 0.0f * I;
			pos = false;
			}
		else {
			a[i] = -( (float)(ORDER - i) + 0.3f ) + 0.0f * I;
			pos = true;
			}
		printf("\n real coeff of order:%3d  = %9.3f", i,  crealf(a[i]) );
		}

	printf("\n\n");

	zroots( a, m, roots, true );  // this function does the actual work

	for( i = 1; i <= ORDER; i++ )	// output complex roots
		printf("\n root:%3d  real = %12.8f   imag = %12.8f", i,  crealf(roots[i]),  cimagf(roots[i]) );

	printf("\n\n");

	return 0;
}
