import numpy
from mcmc import mcmc_framework
#from mcmc_examples.arma.arma_model import ARMAModel
from mcmc_examples.arma.arma_float_model import ARMAFloatModel
#from mcmc_examples.lighthouse.lighthouse_fixed_point_model \
#     import ARMAFixedPointModel

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

# number of posterior samples required per core
n_samples = 100

# mu and sigma values
mu = 0.1
sigma = 0.08

mu_jump_scale = 0.001  # 0.01
sigma_jump_scale = 0.0001 # 0.001

# size of p and q arrays
np = 9
nq = 9

# set up parameters array with p polynomial
# and scaling of t transition distribution for MH jumps in p direction
parameters = []
p_jump_scale = []
for i in range(0,np):
    parameters.append(0.01)
    p_jump_scale.append(0.0001)

# add q polynomial to parameters array
# scaling of t transition distribution for MH jumps in q direction
q_jump_scale = []
for i in range(0,nq):
    parameters.append(0.01)
    q_jump_scale.append(0.0001)

parameters.append(mu)
parameters.append(sigma)

print 'data: ', data_points

print 'state parameter set: ', parameters

# Run and get the samples
#model = ARMAModel(parameters, p_jump_scale, q_jump_scale)
model = ARMAFloatModel(parameters, p_jump_scale, q_jump_scale, mu_jump_scale,
                       sigma_jump_scale)
#model = ARMAFixedPointModel(
#    alpha_jump_scale, alpha_min, alpha_max, beta_jump_scale, beta_min,
#    beta_max)
samples = mcmc_framework.run_mcmc(
    model, data_points, n_samples,
    degrees_of_freedom=3.0, seed=seed, n_chips=3*43,
    root_finder=True)  # n_chips=23*48)

print 'samples: ', samples

# Save the results
numpy.save("results.npy", samples)
# numpy.savetxt("results.csv", samples, fmt="%f", delimiter=",")
numpy.savetxt("results.csv", samples, fmt="%f", delimiter=",")
