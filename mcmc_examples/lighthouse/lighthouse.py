import numpy
from mcmc import mcmc_framework
from mcmc_examples.lighthouse.lighthouse_model import LightHouseModel
#from mcmc_examples.lighthouse.lighthouse_float_model \
#     import LightHouseFloatModel
#from mcmc_examples.lighthouse.lighthouse_fixed_point_model \
#     import LightHouseFixedPointModel

# Data to use for 50 data points
data_points = [
    2.37706, 4.51142, 0.406605, 0.909418, 0.642899,
    1.21925, 1.47647, -2.95771, -0.801802, -1.86529,
    -75.6865, -1.94398, 1.57055, 2.53382, 0.783884,
    1.16725, 1.16995, 0.367477, -1.24639, -2.29897,
    0.461939, -0.126669, 0.0965992, 1.56107, 0.747027,
    66.7057, 1.07821, -0.125864, -0.693059, 7.48744,
    1.94184, -0.439164, 3.64695, -63.2296, 0.783037,
    2.26351, 1.30222, 0.542981, 3.78199, -37.1692,
    1.54959, 0.485336, 1.02509, -0.204211, 0.164426,
    -13.1977, 0.650243, 0.671339, 2.93511, 0.788114
]

# Data to use for 10 data points
#data_points = [
#    2.37706, 4.51142, 0.406605, 0.909418, 0.642899,
#    1.21925, 1.47647, -2.95771, -0.801802, -1.86529
#]

seed = None  # set this if you want to use a different seed on each core
#seed = [  # use this for the same seed on each core
#    123456789, 234567891, 345678912, 456789123, 0
#]

# number of posterior samples required per core
n_samples = 100

# scaling of t transition distribution for MH jumps in alpha direction
alpha_jump_scale = 0.8

# scaling of t transition distribution for MH jumps in beta direction
beta_jump_scale = 0.25

# specification of prior knowledge about lighthouse position
#
# you know that the lighthouse is no further than 3.0 units along shore
# from the reference zero position
#
# lighthouse cannot be closer than 0.2 from the shore because of rocks, or
# further than 2.0 because of shipping lanes
alpha_min = -3.0
alpha_max = 3.0
beta_min = 0.2
beta_max = 2.0

# Run and get the samples
model = LightHouseModel(
    alpha_jump_scale, alpha_min, alpha_max, beta_jump_scale, beta_min,
    beta_max)
#model = LightHouseFloatModel(
#    alpha_jump_scale, alpha_min, alpha_max, beta_jump_scale, beta_min,
#    beta_max)
#model = LightHouseFixedPointModel(
#    alpha_jump_scale, alpha_min, alpha_max, beta_jump_scale, beta_min,
#    beta_max)
samples = mcmc_framework.run_mcmc(
    model, data_points, n_samples,
    degrees_of_freedom=3.0, seed=seed, n_chips=3*44)  # n_chips=23*48)

# Save the results
numpy.save("results.npy", samples)
numpy.savetxt("results.csv", samples, fmt="%f", delimiter=",")
