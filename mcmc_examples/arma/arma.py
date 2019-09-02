import sys
import numpy
from mcmc import mcmc_framework
# from mcmc_examples.arma.arma_model import ARMAModel
from mcmc_examples.arma.arma_float_model import ARMAFloatModel
# from mcmc_examples.lighthouse.lighthouse_fixed_point_model \
#     import ARMAFixedPointModel
from six import iteritems

# Data to use for 1000 data points (read from file)
data_10000 = numpy.loadtxt("data_10000.csv", delimiter=",")

# Data to use for 1000 data points (read from file)
# data_1000 = numpy.loadtxt("data_1000.csv", delimiter=",")

# Edit this number if you want to use less of the data that you've loaded
n_test_points = 10000  # 1000  # 5000
data_points = data_10000[0:n_test_points]
# data_points = data_1000[0:n_test_points]

seed = None  # set this if you want to use a different seed on each core
# seed = [  # use this for the same seed on each core
#    123456789, 234567891, 345678912, 456789123, 0
# ]

# set number of posterior samples to get and number of boards to use
n_samples = 100  # 20000
n_boards = 3

# get n_samples and n_boards from command line arguments if specified
if (len(sys.argv) == 2):
    n_samples = int(sys.argv[1])
elif (len(sys.argv) == 3):
    n_samples = int(sys.argv[1])
    n_boards = int(sys.argv[2])

print("Running ARMA MCMC on ", n_boards, " boards, and collecting ",
      n_samples, " samples")

# mu and sigma values
mu = 0.1
sigma = 0.08

# jump scale values
mu_jump_scale = 0.001  # 0.01
sigma_jump_scale = 0.0001  # 0.001

# size of p and q arrays
np = 9
nq = 9

# set up parameters array with p polynomial
# and scaling of t transition distribution for MH jumps in p direction
parameters = []
jump_scale = []
for i in range(0, np):
    parameters.append(0.01)
    jump_scale.append(0.0001)

# add q polynomial to parameters array
# scaling of t transition distribution for MH jumps in q direction
for i in range(0, nq):
    parameters.append(0.01)
    jump_scale.append(0.0001)

parameters.append(mu)
jump_scale.append(mu_jump_scale)
parameters.append(sigma)
jump_scale.append(sigma_jump_scale)

print('data: ', data_points)
print('initial state parameter set: ', parameters)
print('initial jump scale set: ', jump_scale)

# Run and get the samples
# model = ARMAModel(parameters, jump_scale)
# model = ARMAFloatModel(parameters, jump_scale)
model = ARMAFloatModel(parameters, jump_scale)  # note: this sets both True
# p_jump_scale, q_jump_scale, mu_jump_scale, sigma_jump_scale)
# model = ARMAFixedPointModel(
#    alpha_jump_scale, alpha_min, alpha_max, beta_jump_scale, beta_min,
#    beta_max)

# Uncomment below for a spinn-5 run
# samples = mcmc_framework.run_mcmc(
#    model, data_points, n_samples, burn_in=5000, thinning=50,
#    degrees_of_freedom=6.0, seed=seed, n_boards=n_boards)

# Uncomment below for a spinn-3 run
samples = mcmc_framework.run_mcmc(
    model, data_points, n_samples, burn_in=5000, thinning=50,
    degrees_of_freedom=6.0, seed=seed, n_boards=1)

# print('samples: ', samples)

# Save the results
for coord, sample in iteritems(samples):
    fname = "results_board_x"+str(coord[0])+"_y"+str(
        coord[1])+"_n_boards"+str(n_boards)+"_n_samples"+str(n_samples)
    numpy.save(fname+".npy", sample)
    numpy.savetxt(fname+".csv", sample, fmt="%f", delimiter=",")
