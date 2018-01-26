
#include "spinn_real_type.h"  // REAL -> float, double, accum or whatever
#include "root_finder.h"

/*

This function defines the log of likelihood giventhe data, the model
structure and the parameters.

data is the processed time series data

 order=[p, q] are the orders of ARMA(p,q) model as follows:
 z(t)=a1*z(t-1)+a2*z(t-2)+a3*z(t-3)+...a_p*z(t-p)+e(t)-b1*e(t-1)-b2*e(t-2)-...-b_q*e(t-q)+mu

 parameters are the variables involved in the structure.
 The number of parameters equals p+q+2;
 parameters=[a1,...a_p,b1,...b_q,mu,sigma]

 N is the length of the dataset. Error vector is of length N+q?  Check this

*/
REAL loglikelihood( REAL data[], uint32_t N, REAL err[], uint8_t order[2], REAL parameters )
{
	uint8_t p = order[0], q = order[1], i, index;	// read in AR and MA parameter dimensions

	REAL tempdotp = REAL_CONST( 0.0 ), tempdotq = REAL_CONST( 0.0 );
// mu=parameters(end-1); % mu is the second last element of the parameter.
// sigma=parameters(end); % sigma is the last element of the parameter
	REAL sigma = parameters[ p+q+1 ], mu = parameters[ p+q ];

/*

// err is all zeros with size of (N+q)*1
err=zeros(N+q,1); % this vector will become non-zero from the p+q+1
*/
	uint32_t error_length = N + q;  // check this and all indexing below!

	for( i = 0; i < error_length; i++ )
		err[i] = REAL_CONST( 0.0 );

/*
	Y=zeros(N,1);% this is the predicted output
*/
	REAL Y[N]; // C99 compile flag required

	// compiler issue, could add mu here rather than in the main loop

/*
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
	for( i = p; i < N; i++ ) {
		// loop over p for parameters * data dot product
		for ( j = 0; j < p; j++ ) {
			tempdotp = parameters[j] * data[(i-1)-j]; // check this
		}
		// loop over q for parameters * err dot product
		for ( j = 0; j < q; j++ ) {
			tempdotq = parameters[p+j] * err[(i+q-1)-j]; // check this
		}

		Y[i] = tempdotp + tempdotq + mu;
		err[i+q] = data[i] - Y[i];
		}

/*

% The first part of log of likelihood equals the negative sum of square error and then divided
% by double of square sigma. The second part is half the number of (N-p) times the log of square sigma.
% The substraction between the first part and the second part equals the lglikelihood.

lglikelihood=-sum(err(p+q+1:end).^2)/(2*sigma^2)-0.5*(N-p)*log(sigma^2);

*/
	REAL sum = REAL_CONST( 0.0 ), temp;
	float divisor = 1.0f / (( 2.0f * sigma * sigma ) - 0.5f * (N-p) * logf( sigma * sigma )); // to avoid divide, create one floating point inverse & multiply later

	for( i = p+q; i < N+q; i++ ) { // check
		temp = err[i];
		sum += temp * temp;
		}

	return -sum * (REAL) divisor;  // possibly add some error checking in case divisor is stupid value
}
