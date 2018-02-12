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

uint32_t mcmc_model_get_params_n_bytes() {
    return sizeof(struct mcmc_params);
}

uint32_t mcmc_model_get_state_n_bytes() {
    return sizeof(struct mcmc_state);
}

#if TYPE_SELECT != 2
struct double_uint {
    uint first_word;
    uint second_word;
};

union double_to_ints {
    CALC_TYPE double_value;
    struct double_uint int_values;
};

void print_value_arma(CALC_TYPE d_value, char *buffer) {
    union double_to_ints converter;
    converter.double_value = d_value;
    io_printf(
        buffer, "0x%08x%08x",
        converter.int_values.second_word, converter.int_values.first_word);
}
#endif

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
//void print_value(CALC_TYPE d_value, char *buffer) {
//    union double_to_ints converter;
//    converter.double_value = d_value;
//    io_printf(
//        buffer, "0x%08x%08x",
//        converter.int_values.second_word, converter.int_values.first_word);
//}

volatile CALC_TYPE result_value;

void result_callback(uint key, uint payload) {
//	log_info("ARMA: result_callback");
	result_value = payload;
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
	use(params);

	char buffer[1024];

	// read in AR and MA parameter dimensions
	uint8_t p = PPOLYORDER;  // state->order_p;
	uint8_t q = QPOLYORDER;  // state->order_q;
	uint8_t i, j;
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
			tempdotp = state_parameters[j] * data[(i-1)-j]; // check this
		}
		// loop over q for parameters * err dot product
		for (j=0; j < q; j++) {
			tempdotq = state_parameters[p+j] * err[(i+q-1)-j]; // check this
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

//	print_value_arma(divisor, buffer);
//	log_info("ARMA likelihood, divisor = %s", divisor);

	//float divisor = 1.0f / (( 2.0f * sigma * sigma ) - 0.5f * (N-p) * logf( sigma * sigma )); // to avoid divide, create one floating point inverse & multiply later
	for(i=p+q; i < N+q; i++) { // check
		temp = err[i];
		sum += temp * temp;
	}

//	print_value_arma(sum, buffer);
//	log_info("ARMA likelihood, sum = %s", sum);

	return -sum * (CALC_TYPE) divisor;  // possibly add some error checking in case divisor is stupid value
}

/*
 prior probability for the parameters, using state->order_p and state->order_q

 prior condition is that each root is outside the unit circle
 */
CALC_TYPE mcmc_model_prior_prob(
        mcmc_params_pointer_t params, mcmc_state_pointer_t state) {
	// debug for writing values
	char buffer[1024];
	result_value = 2.0f;  // 1.0f;

	// read in AR and MA parameter dimensions
	uint8_t p = PPOLYORDER;  // state->order_p;
	uint8_t q = QPOLYORDER;  // state->order_q;
//	uint8_t i;

//	REAL sigma = params->sigma; // could plausibly use this rather than a parameters array?
	CALC_TYPE *state_parameters;
	state_parameters = state->parameters;
	CALC_TYPE sigma = state_parameters[p+q+1];  // last entry in vector

//	log_info("ARMA: state_parameters[0] = 0x%08x", state_parameters[0]);

	if( sigma <= ZERO ) return ROOT_FAIL;  // first fail condition can provide early exit

	address_t data_address = data_specification_get_data_address();
	address_t parameters_address = data_specification_get_region(
	        PARAMETERS, data_address);
	struct parameters *sdram_params = (struct parameters *) parameters_address;
	spin1_memcpy(&parameters, sdram_params, sizeof(struct parameters));
	uint32_t key = parameters.key;
//	log_info("Key = 0x%08x", key);

	// send key to wake up root finder
//	spin1_send_mc_packet(key, 0, 0);

	// Pointer to receive the data with
//	uint32_t *param_receive_ptr;

	// write parameters to sdram
	uint32_t state_n_bytes = mcmc_model_get_state_n_bytes();
//	log_info("ARMA: state_n_bytes = %d", state_n_bytes);

//	log_info("ARMA: state_n_bytes: %d", state_n_bytes);

	address_t model_state_address = data_specification_get_region(
	        MODEL_STATE, data_address);

    // copy state parameters to sdram
	spin1_memcpy(model_state_address, state_parameters, state_n_bytes);

	// send mc packet to wake up root finder
	// Send mc packet with payload of address value
	spin1_send_mc_packet(key, model_state_address, WITH_PAYLOAD);

//	log_info("ARMA: model_state_address = %d", model_state_address);

//	print_value(sigma, buffer);
//	log_info("ARMA: sigma = %s", buffer);

//	log_info("ARMA: parameters_address = %d", parameters_address);


//	param_receive_ptr = (uint32_t *) sark_xalloc(
//			sv->sdram_heap, (PPOLYORDER + QPOLYORDER + 2) * sizeof(CALC_TYPE),
//	        0, ALLOC_LOCK);

//	log_info("ARMA: model_state_address = %d", model_state_address);

	uint mode = spin1_int_enable();

	// callback for result here
	spin1_callback_off(MCPL_PACKET_RECEIVED);
	//spin1_callback_off(MC_PACKET_RECEIVED);
	spin1_callback_on(MCPL_PACKET_RECEIVED, result_callback, -1);

	// wait for result to come back
	//spin1_callback_on(MCPL_PACKET_RECEIVED, result_callback, -1);

//	print_value_arma(result_value, buffer);
//	log_info("callback turned on, now going in to wait... %s", buffer);

	// do we need to sit here and wait and do nothing until result is here?
	while (result_value==2.0f) {
//		print_value(result_value, buffer);
//		log_info("while loop, result_value"); //  = %s", buffer);
//		log_info("what is in result_value... %f", result_value);
		spin1_wfi();
	}

//	print_value_arma(result_value, buffer);
//	log_info("ARMA: returning result from root finder: %s", buffer);

	spin1_mode_restore(mode);

	// read result and return it
	CALC_TYPE returnval = result_value;
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
	// - it'll look something like this...
	CALC_TYPE *parameters = state->parameters;
//	CALC_TYPE *new_parameters = new_state->parameters;
	uint32_t p = PPOLYORDER;  // state->order_p;
	uint32_t q = QPOLYORDER;  // state->order_q;
	unsigned int i;
	for (i=0; i < p; i++) {
		new_state->parameters[i] = parameters[i] +
				(t_deviate() * params->p_jump_scale[i]);
	}
	for (i=p; i < p+q; i++) {
		new_state->parameters[i] = parameters[i] +
				(t_deviate() * params->q_jump_scale[i]);
	}
	new_state->parameters[p+q] = parameters[p+q] +
			(t_deviate() * params->mu_jump_scale);
	new_state->parameters[p+q+1] = parameters[p+q+1] +
			(t_deviate() * params->sigma_jump_scale);
}