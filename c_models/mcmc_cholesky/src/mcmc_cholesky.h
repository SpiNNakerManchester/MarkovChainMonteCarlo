#include "../../mcmc_models/examples/arma/src/arma.h"
#include <stdbool.h>

//#define FOR( i, n ) for ( i = 0; i < (n); i++ )
#define LA_TYPE double  // float
#define LA_ZERO 0.0  // 0.0f
#define LA_ONE 1.0  // 1.0f
#define LA_QUARTER 0.25  // 0.25f
#define LA_SMALL 1.0e-300  // 1.0e-30f
#define LA_SQRT sqrt  // sqrtf
#define LA_SCALE 0.25  // LA_QUARTER
#define NCOVSAMPLES 5000  // 30000
#define MATDIM PPOLYORDER+QPOLYORDER+2

// Matrices and vectors are defined by the number of state parameters
typedef LA_TYPE Vec[MATDIM];
typedef LA_TYPE Mat[MATDIM][MATDIM];

// functions required to do cholesky decomposition
void zero_upper_triang(Mat A, uint32_t size);

void cholesky(Mat A, const uint32_t size, bool zero_upper);

//void vec_times_mat_scaled(const Vec restrict vec, const Mat restrict mat,
//		Vec res, const float scale, const uint32_t size);
void vec_times_mat_scaled(const Vec vec, const Mat mat, Vec res,
		const LA_TYPE scale, const uint32_t size);

void mean_covar_of_mat_n(const uint32_t n, const uint32_t d);
