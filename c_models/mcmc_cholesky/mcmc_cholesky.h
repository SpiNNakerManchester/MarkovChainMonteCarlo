// compile with gcc -std=c99 -O3 -fcx-limited-range root_finder.c -l m -o rf

#define FOR( i, n ) for( i = 0; i < (n); i++ )

// Guessing that I need to define Mat and Vec somewhere too

// functions required to do cholesky decomposition
void zero_upper_triang(Mat A, uint32_t size);

void cholesky(Mat A, const uint32_t size, bool zero_upper);

void vec_times_mat_scaled(const Vec restrict vec, const Mat restrict mat,
		Vec res, const float scale, const uint32_t size);

void mean_covar_of_mat_n(const Mat data, Vec mean, Mat cov,
		const uint32_t n, const uint32_t d);

