#include "../../mcmc_model.h"
#include "arma.h"
#include <spin1_api.h>
#include <debug.h>
#include <data_specification.h>

#define ROOT_FAIL CONCAT(-1000.000000, SUFFIX)		// how bad is a root failure
#define REAL float

enum regions {
    RECORDING,
    PARAMETERS,
    MODEL_PARAMETERS,
    MODEL_STATE
};

// The type of the seed
typedef uint32_t uniform_seed[5];

struct parameters {

    // no of MCMC transitions to reach apparent equilibrium before generating
    // inference samples
    uint32_t burn_in;

    // subsequent MCMC samples are correlated, so thin the chain to avoid this
    uint32_t thinning;

    // number of posterior samples required
    uint32_t n_samples;

    // number of data points
    uint32_t n_data_points;

    // data receive window size - 0 if not a receiver
    uint32_t data_window_size;

    // Sequence mask for data reception - 0 if not a receiver
    uint32_t sequence_mask;

    // Acknowledge key for data reception - 0 if not a receiver
    uint32_t acknowledge_key;

    // Tag which is allocated for the data
    uint32_t data_tag;

    // Timer for data acknowledgement - 0 if not receiver
    uint32_t timer;

    // Key for sending parameters
    uint32_t key;

    // The random seed
    uniform_seed seed;

    // The number of degrees of freedom to jump around with
    CALC_TYPE degrees_of_freedom;
};

// The general parameters
struct parameters parameters;

// Define spin1_wfi
extern void spin1_wfi();

// Key
uint32_t key;

// Model state address
address_t model_state_address;

// size of model state
uint32_t state_n_bytes;

uint32_t mcmc_model_get_params_n_bytes() {
    return sizeof(struct mcmc_params);
}

uint32_t mcmc_model_get_state_n_bytes() {
    return sizeof(struct mcmc_state);
}

//#if TYPE_SELECT != 2
//struct double_uint {
//    uint first_word;
//    uint second_word;
//};
//
//union double_to_ints {
//    CALC_TYPE double_value;
//    struct double_uint int_values;
//};
//
//void print_value_arma(CALC_TYPE d_value, char *buffer) {
//    union double_to_ints converter;
//    converter.double_value = d_value;
//    io_printf(
//        buffer, "0x%08x%08x",
//        converter.int_values.second_word, converter.int_values.first_word);
//}
//#endif

uint8_t result_value;

void result_callback(uint key, uint payload) {
//	log_info("ARMA: result_callback");
	use(key);
	result_value = payload;
}

/*
 ARMA wind log likelihood
  - data, size n_pts
  - p and q are the orders of an ARMA(p,q) model, which has the form
    z(t)=a_1*z(t-1)+...+a_p*z(t-p)-b_1*e(t-1)-...-b_q*e(t-q)+mu
  - params is [a1,...,ap,b1,...,bq,mu,sigma]
 */
CALC_TYPE mcmc_model_likelihood(
        CALC_TYPE *data, uint32_t n_pts, mcmc_params_pointer_t params,
		mcmc_state_pointer_t state) {
	use(params);

//	char buffer[1024];

	// read in AR and MA parameter dimensions
	uint8_t p = PPOLYORDER;  // state->order_p;
	uint8_t q = QPOLYORDER;  // state->order_q;
	uint8_t j;
	uint32_t i;
	uint32_t N = n_pts;

	uint32_t error_length = N + q;  // check this and all indexing below!
	CALC_TYPE err[error_length];

	// get parameters array from state struct
	CALC_TYPE *state_parameters = state->parameters;

	// temp storage for dot products in later loop
	CALC_TYPE tempdotp = ZERO; // REAL_CONST( 0.0 );
	CALC_TYPE tempdotq = ZERO; // REAL_CONST( 0.0 );

	// mu is the penultimate element of the parameters
	// sigma is the last element of the parameters
	CALC_TYPE mu = state_parameters[p+q];
	CALC_TYPE sigma = state_parameters[p+q+1];

/*
err is all zeros with size of (N+q)*1
err=zeros(N+q,1); % this vector will become non-zero from the p+q+1
*/
	for(i=0; i < error_length; i++) {
		err[i] = ZERO;  // REAL_CONST( 0.0 );
	}

/*
	Y=zeros(N,1);% this is the predicted output
*/
	//REAL Y[N]; // C99 compile flag required
	CALC_TYPE Y;  // this is only used inside the upcoming loop

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
			tempdotp += state_parameters[j] * data[(i-1)-j]; // check this
		}
		// loop over q for parameters * err dot product
		for (j=0; j < q; j++) {
			tempdotq += state_parameters[p+j] * err[(i+q-1)-j]; // check this
		}

		// Combine two results together plus mean
		Y = tempdotp - tempdotq + mu;
		// Error is data minus predicted data
		err[i+q] = data[i] - Y;

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
	CALC_TYPE denominator = TWO*SQR(sigma);
	CALC_TYPE Nminusp = (CALC_TYPE) (N-p);
	CALC_TYPE temp2 = HALF*Nminusp*LN(SQR(sigma));
	CALC_TYPE divisor = ONE / denominator;
	// Check value of divisor here (similar to in transition_jump??

	// Loop to do sum
	for(i=p+q; i < N+q; i++) { // check
		temp = err[i];
		sum += temp * temp;  // might be better using SQR here?
	}

	// possibly add some error checking in case divisor is stupid value
	return (-sum * divisor) - temp2;
}

/*
 prior probability for the parameters, using state->order_p and state->order_q

 prior condition is that each root is outside the unit circle
 */
CALC_TYPE mcmc_model_prior_prob(
        mcmc_params_pointer_t params, mcmc_state_pointer_t state) {
	// debug for writing values
//	char buffer[1024];
	// set result_value in order to wait for result from root_finder
	result_value = 2;  // 1.0f;

	// read in AR and MA parameter dimensions
	uint8_t p = PPOLYORDER;  // state->order_p;
	uint8_t q = QPOLYORDER;  // state->order_q;

	// Get the sigma value
	CALC_TYPE *state_parameters;
	state_parameters = state->parameters;
	CALC_TYPE sigma = state_parameters[p+q+1];  // last entry in vector

	// If sigma is less than zero then we can exit without using root_finder
	if( sigma <= ZERO ) return ROOT_FAIL;

    // If we're still going, we need to copy state parameters to sdram
	spin1_memcpy(model_state_address, state_parameters, state_n_bytes);

	// Send mc packet with payload of address value to wake up root_finder
	spin1_send_mc_packet(key, model_state_address, WITH_PAYLOAD);

	// Wait here for the result to come back
	while (result_value==2) {
		spin1_wfi();
	}

	// Read the result and return it
	CALC_TYPE returnval = ZERO;
	if (result_value==1) returnval = ROOT_FAIL;
	return returnval;

}


/*
 generate jumps from MH transition distribution which is uncorrelated bivariate
 t with degrees_freedom degrees of freedom
 */
void mcmc_model_transition_jump(
        mcmc_params_pointer_t params, mcmc_state_pointer_t state,
        mcmc_state_pointer_t new_state) {
	// loop over parameters and apply relevant jump_scale
//	char buffer[1024];

	// NOTE (16/2/18): this needs to use Cholesky decomposition to
	// update the jump scale parameters

	CALC_TYPE *parameters = state->parameters;
	uint32_t p = PPOLYORDER;  // state->order_p;
	uint32_t q = QPOLYORDER;  // state->order_q;
	unsigned int i;

	// Update polynomial coefficients
	for (i=0; i < p; i++) {
		new_state->parameters[i] = parameters[i] +
				(t_deviate() * params->p_jump_scale[i]);
	}
	for (i=p; i < p+q; i++) {
		new_state->parameters[i] = parameters[i] +
				(t_deviate() * params->q_jump_scale[i-p]);
	}

	// Update mu and sigma
	new_state->parameters[p+q] = parameters[p+q] +
			(t_deviate() * params->mu_jump_scale);
	new_state->parameters[p+q+1] = parameters[p+q+1] +
			(t_deviate() * params->sigma_jump_scale);
}

/*
 exit function: send message to root finder to exit, then exit
 */
void mcmc_exit_function() {
	// send message to root finder
	spin1_send_mc_packet(key, 0, 0);

	// exit for this core
	spin1_exit(0);
}

/*
 set up addresses for sending information to root finder
 */
void mcmc_get_address_and_key() {
	// get address from data spec
	address_t data_address = data_specification_get_data_address();

	// get key
	address_t parameters_address = data_specification_get_region(
	        PARAMETERS, data_address);
	struct parameters *sdram_params = (struct parameters *) parameters_address;
	spin1_memcpy(&parameters, sdram_params, sizeof(struct parameters));
	key = parameters.key;

	// store size of state parameters and model state address
	state_n_bytes = mcmc_model_get_state_n_bytes();
	model_state_address = data_specification_get_region(
	        MODEL_STATE, data_address);

	// set up callback for results from root finder
	spin1_callback_on(MCPL_PACKET_RECEIVED, result_callback, -1);

}
