// compile with gcc -std=c99 -O3 -fcx-limited-range root_finder.c -l m -o rf
#include "../mcmc_models/examples/arma/arma.h"
#include <stdbool.h>

#define FOR( i, n ) for( i = 0; i < (n); i++ )
#define LA_TYPE double  // float
#define LA_ZERO 0.0  // 0.0f
#define LA_ONE 1.0  // 1.0f
#define LA_QUARTER 0.25  // 0.25f
#define LA_SMALL 1.0e-300  // 1.0e-30f
#define LA_SQRT sqrt  // sqrtf
#define NCOVSAMPLES 5000

// Matrices and vectors are mainly defined by the number of state parameters
typedef LA_TYPE Mat[PPOLYORDER+QPOLYORDER+2][PPOLYORDER+QPOLYORDER+2];
typedef LA_TYPE Vec[PPOLYORDER+QPOLYORDER+2];
// Except for this one which is defined by the number of samples
// This needs a rethink if (approx) NCOVSAMPLES > 800 as it won't fit in DTCM
// MH suggests (and JH's Matlab code uses) a value of NCOVSAMPLES = 5000
//typedef float DataMat[NCOVSAMPLES][PPOLYORDER+QPOLYORDER+2];

// functions required to do cholesky decomposition
void zero_upper_triang(Mat A, uint32_t size);

void cholesky(Mat A, const uint32_t size, bool zero_upper);

//void vec_times_mat_scaled(const Vec restrict vec, const Mat restrict mat,
//		Vec res, const float scale, const uint32_t size);
void vec_times_mat_scaled(const Vec vec, const Mat mat, Vec res,
		const LA_TYPE scale, const uint32_t size);

void mean_covar_of_mat_n(float **data, Vec mean, Mat cov,
		const uint32_t n, const uint32_t d);
