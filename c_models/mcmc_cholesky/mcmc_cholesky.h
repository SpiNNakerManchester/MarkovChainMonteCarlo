// compile with gcc -std=c99 -O3 -fcx-limited-range root_finder.c -l m -o rf
#include "../mcmc_models/examples/arma/arma.h"
#include <stdbool.h>

#define FOR( i, n ) for( i = 0; i < (n); i++ )

// Guessing that I need to define Mat and Vec somewhere here
#define NCOVSAMPLES 700

// Matrices and vectors are mainly defined by the number of state parameters
typedef float Mat[PPOLYORDER+QPOLYORDER+2][PPOLYORDER+QPOLYORDER+2];
typedef float Vec[PPOLYORDER+QPOLYORDER+2];
// Except for this one which is defined by the number of samples
typedef float DataMat[NCOVSAMPLES][PPOLYORDER+QPOLYORDER+2];

// functions required to do cholesky decomposition
void zero_upper_triang(Mat A, uint32_t size);

void cholesky(Mat A, const uint32_t size, bool zero_upper);

//void vec_times_mat_scaled(const Vec restrict vec, const Mat restrict mat,
//		Vec res, const float scale, const uint32_t size);
void vec_times_mat_scaled(const Vec vec, const Mat mat, Vec res,
		const float scale, const uint32_t size);

void mean_covar_of_mat_n(const DataMat data, Vec mean, Mat cov,
		const uint32_t n, const uint32_t d);
