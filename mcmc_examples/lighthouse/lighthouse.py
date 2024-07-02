# Copyright (c) 2016 The University of Manchester
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     https://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import sys
import os
from time import gmtime, strftime
import pathos.multiprocessing  # type: ignore
import numpy
from mcmc import mcmc_framework
# from mcmc_examples.lighthouse.lighthouse_model import LightHouseModel
from mcmc_examples.lighthouse.lighthouse_float_model \
     import LightHouseFloatModel
# from mcmc_examples.lighthouse.lighthouse_fixed_point_model \
#      import LightHouseFixedPointModel

# Data to use for 50 data points
data_50 = [
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
data_10 = [
    2.37706, 4.51142, 0.406605, 0.909418, 0.642899,
    1.21925, 1.47647, -2.95771, -0.801802, -1.86529
]

# Data to use for 1000 data points (for testing purposes to check DMA method)
# (This is simply the 50 data points repeated 20 times).
data_1000 = [
    2.37706, 4.51142, 0.406605, 0.909418, 0.642899,
    1.21925, 1.47647, -2.95771, -0.801802, -1.86529,
    -75.6865, -1.94398, 1.57055, 2.53382, 0.783884,
    1.16725, 1.16995, 0.367477, -1.24639, -2.29897,
    0.461939, -0.126669, 0.0965992, 1.56107, 0.747027,
    66.7057, 1.07821, -0.125864, -0.693059, 7.48744,
    1.94184, -0.439164, 3.64695, -63.2296, 0.783037,
    2.26351, 1.30222, 0.542981, 3.78199, -37.1692,
    1.54959, 0.485336, 1.02509, -0.204211, 0.164426,
    -13.1977, 0.650243, 0.671339, 2.93511, 0.788114,
    2.37706, 4.51142, 0.406605, 0.909418, 0.642899,
    1.21925, 1.47647, -2.95771, -0.801802, -1.86529,
    -75.6865, -1.94398, 1.57055, 2.53382, 0.783884,
    1.16725, 1.16995, 0.367477, -1.24639, -2.29897,
    0.461939, -0.126669, 0.0965992, 1.56107, 0.747027,
    66.7057, 1.07821, -0.125864, -0.693059, 7.48744,
    1.94184, -0.439164, 3.64695, -63.2296, 0.783037,
    2.26351, 1.30222, 0.542981, 3.78199, -37.1692,
    1.54959, 0.485336, 1.02509, -0.204211, 0.164426,
    -13.1977, 0.650243, 0.671339, 2.93511, 0.788114,
    2.37706, 4.51142, 0.406605, 0.909418, 0.642899,
    1.21925, 1.47647, -2.95771, -0.801802, -1.86529,
    -75.6865, -1.94398, 1.57055, 2.53382, 0.783884,
    1.16725, 1.16995, 0.367477, -1.24639, -2.29897,
    0.461939, -0.126669, 0.0965992, 1.56107, 0.747027,
    66.7057, 1.07821, -0.125864, -0.693059, 7.48744,
    1.94184, -0.439164, 3.64695, -63.2296, 0.783037,
    2.26351, 1.30222, 0.542981, 3.78199, -37.1692,
    1.54959, 0.485336, 1.02509, -0.204211, 0.164426,
    -13.1977, 0.650243, 0.671339, 2.93511, 0.788114,
    2.37706, 4.51142, 0.406605, 0.909418, 0.642899,
    1.21925, 1.47647, -2.95771, -0.801802, -1.86529,
    -75.6865, -1.94398, 1.57055, 2.53382, 0.783884,
    1.16725, 1.16995, 0.367477, -1.24639, -2.29897,
    0.461939, -0.126669, 0.0965992, 1.56107, 0.747027,
    66.7057, 1.07821, -0.125864, -0.693059, 7.48744,
    1.94184, -0.439164, 3.64695, -63.2296, 0.783037,
    2.26351, 1.30222, 0.542981, 3.78199, -37.1692,
    1.54959, 0.485336, 1.02509, -0.204211, 0.164426,
    -13.1977, 0.650243, 0.671339, 2.93511, 0.788114,
    2.37706, 4.51142, 0.406605, 0.909418, 0.642899,
    1.21925, 1.47647, -2.95771, -0.801802, -1.86529,
    -75.6865, -1.94398, 1.57055, 2.53382, 0.783884,
    1.16725, 1.16995, 0.367477, -1.24639, -2.29897,
    0.461939, -0.126669, 0.0965992, 1.56107, 0.747027,
    66.7057, 1.07821, -0.125864, -0.693059, 7.48744,
    1.94184, -0.439164, 3.64695, -63.2296, 0.783037,
    2.26351, 1.30222, 0.542981, 3.78199, -37.1692,
    1.54959, 0.485336, 1.02509, -0.204211, 0.164426,
    -13.1977, 0.650243, 0.671339, 2.93511, 0.788114,
    2.37706, 4.51142, 0.406605, 0.909418, 0.642899,
    1.21925, 1.47647, -2.95771, -0.801802, -1.86529,
    -75.6865, -1.94398, 1.57055, 2.53382, 0.783884,
    1.16725, 1.16995, 0.367477, -1.24639, -2.29897,
    0.461939, -0.126669, 0.0965992, 1.56107, 0.747027,
    66.7057, 1.07821, -0.125864, -0.693059, 7.48744,
    1.94184, -0.439164, 3.64695, -63.2296, 0.783037,
    2.26351, 1.30222, 0.542981, 3.78199, -37.1692,
    1.54959, 0.485336, 1.02509, -0.204211, 0.164426,
    -13.1977, 0.650243, 0.671339, 2.93511, 0.788114,
    2.37706, 4.51142, 0.406605, 0.909418, 0.642899,
    1.21925, 1.47647, -2.95771, -0.801802, -1.86529,
    -75.6865, -1.94398, 1.57055, 2.53382, 0.783884,
    1.16725, 1.16995, 0.367477, -1.24639, -2.29897,
    0.461939, -0.126669, 0.0965992, 1.56107, 0.747027,
    66.7057, 1.07821, -0.125864, -0.693059, 7.48744,
    1.94184, -0.439164, 3.64695, -63.2296, 0.783037,
    2.26351, 1.30222, 0.542981, 3.78199, -37.1692,
    1.54959, 0.485336, 1.02509, -0.204211, 0.164426,
    -13.1977, 0.650243, 0.671339, 2.93511, 0.788114,
    2.37706, 4.51142, 0.406605, 0.909418, 0.642899,
    1.21925, 1.47647, -2.95771, -0.801802, -1.86529,
    -75.6865, -1.94398, 1.57055, 2.53382, 0.783884,
    1.16725, 1.16995, 0.367477, -1.24639, -2.29897,
    0.461939, -0.126669, 0.0965992, 1.56107, 0.747027,
    66.7057, 1.07821, -0.125864, -0.693059, 7.48744,
    1.94184, -0.439164, 3.64695, -63.2296, 0.783037,
    2.26351, 1.30222, 0.542981, 3.78199, -37.1692,
    1.54959, 0.485336, 1.02509, -0.204211, 0.164426,
    -13.1977, 0.650243, 0.671339, 2.93511, 0.788114,
    2.37706, 4.51142, 0.406605, 0.909418, 0.642899,
    1.21925, 1.47647, -2.95771, -0.801802, -1.86529,
    -75.6865, -1.94398, 1.57055, 2.53382, 0.783884,
    1.16725, 1.16995, 0.367477, -1.24639, -2.29897,
    0.461939, -0.126669, 0.0965992, 1.56107, 0.747027,
    66.7057, 1.07821, -0.125864, -0.693059, 7.48744,
    1.94184, -0.439164, 3.64695, -63.2296, 0.783037,
    2.26351, 1.30222, 0.542981, 3.78199, -37.1692,
    1.54959, 0.485336, 1.02509, -0.204211, 0.164426,
    -13.1977, 0.650243, 0.671339, 2.93511, 0.788114,
    2.37706, 4.51142, 0.406605, 0.909418, 0.642899,
    1.21925, 1.47647, -2.95771, -0.801802, -1.86529,
    -75.6865, -1.94398, 1.57055, 2.53382, 0.783884,
    1.16725, 1.16995, 0.367477, -1.24639, -2.29897,
    0.461939, -0.126669, 0.0965992, 1.56107, 0.747027,
    66.7057, 1.07821, -0.125864, -0.693059, 7.48744,
    1.94184, -0.439164, 3.64695, -63.2296, 0.783037,
    2.26351, 1.30222, 0.542981, 3.78199, -37.1692,
    1.54959, 0.485336, 1.02509, -0.204211, 0.164426,
    -13.1977, 0.650243, 0.671339, 2.93511, 0.788114,
    2.37706, 4.51142, 0.406605, 0.909418, 0.642899,
    1.21925, 1.47647, -2.95771, -0.801802, -1.86529,
    -75.6865, -1.94398, 1.57055, 2.53382, 0.783884,
    1.16725, 1.16995, 0.367477, -1.24639, -2.29897,
    0.461939, -0.126669, 0.0965992, 1.56107, 0.747027,
    66.7057, 1.07821, -0.125864, -0.693059, 7.48744,
    1.94184, -0.439164, 3.64695, -63.2296, 0.783037,
    2.26351, 1.30222, 0.542981, 3.78199, -37.1692,
    1.54959, 0.485336, 1.02509, -0.204211, 0.164426,
    -13.1977, 0.650243, 0.671339, 2.93511, 0.788114,
    2.37706, 4.51142, 0.406605, 0.909418, 0.642899,
    1.21925, 1.47647, -2.95771, -0.801802, -1.86529,
    -75.6865, -1.94398, 1.57055, 2.53382, 0.783884,
    1.16725, 1.16995, 0.367477, -1.24639, -2.29897,
    0.461939, -0.126669, 0.0965992, 1.56107, 0.747027,
    66.7057, 1.07821, -0.125864, -0.693059, 7.48744,
    1.94184, -0.439164, 3.64695, -63.2296, 0.783037,
    2.26351, 1.30222, 0.542981, 3.78199, -37.1692,
    1.54959, 0.485336, 1.02509, -0.204211, 0.164426,
    -13.1977, 0.650243, 0.671339, 2.93511, 0.788114,
    2.37706, 4.51142, 0.406605, 0.909418, 0.642899,
    1.21925, 1.47647, -2.95771, -0.801802, -1.86529,
    -75.6865, -1.94398, 1.57055, 2.53382, 0.783884,
    1.16725, 1.16995, 0.367477, -1.24639, -2.29897,
    0.461939, -0.126669, 0.0965992, 1.56107, 0.747027,
    66.7057, 1.07821, -0.125864, -0.693059, 7.48744,
    1.94184, -0.439164, 3.64695, -63.2296, 0.783037,
    2.26351, 1.30222, 0.542981, 3.78199, -37.1692,
    1.54959, 0.485336, 1.02509, -0.204211, 0.164426,
    -13.1977, 0.650243, 0.671339, 2.93511, 0.788114,
    2.37706, 4.51142, 0.406605, 0.909418, 0.642899,
    1.21925, 1.47647, -2.95771, -0.801802, -1.86529,
    -75.6865, -1.94398, 1.57055, 2.53382, 0.783884,
    1.16725, 1.16995, 0.367477, -1.24639, -2.29897,
    0.461939, -0.126669, 0.0965992, 1.56107, 0.747027,
    66.7057, 1.07821, -0.125864, -0.693059, 7.48744,
    1.94184, -0.439164, 3.64695, -63.2296, 0.783037,
    2.26351, 1.30222, 0.542981, 3.78199, -37.1692,
    1.54959, 0.485336, 1.02509, -0.204211, 0.164426,
    -13.1977, 0.650243, 0.671339, 2.93511, 0.788114,
    2.37706, 4.51142, 0.406605, 0.909418, 0.642899,
    1.21925, 1.47647, -2.95771, -0.801802, -1.86529,
    -75.6865, -1.94398, 1.57055, 2.53382, 0.783884,
    1.16725, 1.16995, 0.367477, -1.24639, -2.29897,
    0.461939, -0.126669, 0.0965992, 1.56107, 0.747027,
    66.7057, 1.07821, -0.125864, -0.693059, 7.48744,
    1.94184, -0.439164, 3.64695, -63.2296, 0.783037,
    2.26351, 1.30222, 0.542981, 3.78199, -37.1692,
    1.54959, 0.485336, 1.02509, -0.204211, 0.164426,
    -13.1977, 0.650243, 0.671339, 2.93511, 0.788114,
    2.37706, 4.51142, 0.406605, 0.909418, 0.642899,
    1.21925, 1.47647, -2.95771, -0.801802, -1.86529,
    -75.6865, -1.94398, 1.57055, 2.53382, 0.783884,
    1.16725, 1.16995, 0.367477, -1.24639, -2.29897,
    0.461939, -0.126669, 0.0965992, 1.56107, 0.747027,
    66.7057, 1.07821, -0.125864, -0.693059, 7.48744,
    1.94184, -0.439164, 3.64695, -63.2296, 0.783037,
    2.26351, 1.30222, 0.542981, 3.78199, -37.1692,
    1.54959, 0.485336, 1.02509, -0.204211, 0.164426,
    -13.1977, 0.650243, 0.671339, 2.93511, 0.788114,
    2.37706, 4.51142, 0.406605, 0.909418, 0.642899,
    1.21925, 1.47647, -2.95771, -0.801802, -1.86529,
    -75.6865, -1.94398, 1.57055, 2.53382, 0.783884,
    1.16725, 1.16995, 0.367477, -1.24639, -2.29897,
    0.461939, -0.126669, 0.0965992, 1.56107, 0.747027,
    66.7057, 1.07821, -0.125864, -0.693059, 7.48744,
    1.94184, -0.439164, 3.64695, -63.2296, 0.783037,
    2.26351, 1.30222, 0.542981, 3.78199, -37.1692,
    1.54959, 0.485336, 1.02509, -0.204211, 0.164426,
    -13.1977, 0.650243, 0.671339, 2.93511, 0.788114,
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

data_points = data_50

print('data_points: ', data_points)

seed = None  # set this if you want to use a different seed on each core
# seed = [  # use this for the same seed on each core
#    123456789, 234567891, 345678912, 456789123, 0
# ]

# set number of posterior samples to get and number of boards to use
# and number of threads to run (so this will run n_threads jobs each using
# n_boards boards, and collect n_samples samples)
n_samples = 100  # 100 is the "default"
n_boards = 1
n_threads = 1

# get n_samples and n_boards from command line arguments if specified
if (len(sys.argv) == 2):
    if sys.argv[1] != 'test_scripts.py':
        n_samples = int(sys.argv[1])
elif (len(sys.argv) == 3):
    n_samples = int(sys.argv[1])
    n_boards = int(sys.argv[2])
elif (len(sys.argv) == 4):
    n_samples = int(sys.argv[1])
    n_boards = int(sys.argv[2])
    n_threads = int(sys.argv[3])

print("Running MCMC lighthouse on ", n_boards, " boards, and collecting ",
      n_samples, " samples")

# scaling of t transition distribution for jumps in alpha direction
alpha_jump_scale = 0.8

# scaling of t transition distribution for jumps in beta direction
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
# model = LightHouseModel(
#    alpha_jump_scale, alpha_min, alpha_max, beta_jump_scale, beta_min,
#    beta_max)
model = LightHouseFloatModel(
   alpha_jump_scale, alpha_min, alpha_max, beta_jump_scale, beta_min,
   beta_max)
# model = LightHouseFixedPointModel(
#     alpha_jump_scale, alpha_min, alpha_max, beta_jump_scale, beta_min,
#     beta_max)


def run_job(_thread_id, _model=model, _data_points=None,
            _n_samples=n_samples, _seed=seed):
    """
    Main method to run once or in multiple threads
    """
    if _data_points is None:
        _data_points = list(data_points)
    samples = mcmc_framework.run_mcmc(
        _model, _data_points, _n_samples,
        degrees_of_freedom=3.0, seed=_seed, n_boards=n_boards)

    print('samples: ', samples)

    dirpath = (f'results_{strftime("%Y-%m-%d_%H:%M:%S", gmtime())}'
               f'_nboards{n_boards}_nsamples{_n_samples}')
    os.mkdir(dirpath)
    for coord, sample in samples.items():
        fname = (f"{dirpath}/results_th{_thread_id[0]}"
                 f"_board_x{coord[0]}_y{coord[1]}_nboards{n_boards}"
                 f"_nsamples{_n_samples}")
        numpy.save(fname+".npy", sample)
        numpy.savetxt(fname+".csv", sample, fmt="%f", delimiter=",")


# run threaded if requested
if (n_threads == 1):
    # simply call the function run_job, don't run with threads
    run_job([0], model, data_points, n_samples, seed)
else:
    connection_threads = [[n, model, data_points, n_samples, seed]
                          for n in range(n_threads)]

    pool = pathos.multiprocessing.Pool(processes=n_threads)
    pool.map(func=run_job, iterable=connection_threads)

    print("exit main thread")
