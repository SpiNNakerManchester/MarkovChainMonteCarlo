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

from typing import List
from spinn_utilities.overrides import overrides

from mcmc.mcmc_model import MCMCModel
from mcmc.mcmc_parameter import MCMCParameter
from mcmc.mcmc_state_variable import MCMCStateVariable

import numpy


class LightHouseFloatModel(MCMCModel):
    """ MCMC Model for the lighthouse problem, using float(32)
    """

    def __init__(
            self, alpha_jump_scale, alpha_min, alpha_max, beta_jump_scale,
            beta_min, beta_max, root_finder=False, cholesky=False):
        """
        :param alpha_jump_scale:\
            scaling of t transition distribution for MH jumps in alpha\
            direction
        :param alpha_min: The minimum value of alpha
        :param alpha_max: The maximum value of alpha
        :param beta_jump_scale:\
            scaling of t transition distribution for MH jumps in beta direction
        :param beta_min: The minimum value of beta
        :param beta_max: The maximum value of beta
        """

        self._alpha_jump_scale = alpha_jump_scale
        self._alpha_min = alpha_min
        self._alpha_max = alpha_max
        self._beta_jump_scale = beta_jump_scale
        self._beta_min = beta_min
        self._beta_max = beta_max
        self._root_finder = root_finder
        self._cholesky = cholesky

    @overrides(MCMCModel.get_binary_name)
    def get_binary_name(self) -> str:
        return "lighthouse.aplx"

    @overrides(MCMCModel.get_parameters)
    def get_parameters(self) -> List[MCMCParameter]:
        return [
            MCMCParameter(self._alpha_jump_scale, numpy.float32),
            MCMCParameter(self._beta_jump_scale, numpy.float32),
            MCMCParameter(self._alpha_min, numpy.float32),
            MCMCParameter(self._alpha_max, numpy.float32),
            MCMCParameter(self._beta_min, numpy.float32),
            MCMCParameter(self._beta_max, numpy.float32)
        ]

    @overrides(MCMCModel.get_state_variables)
    def get_state_variables(self) -> List[MCMCStateVariable]:
        return [
            MCMCStateVariable("alpha", 0.0, numpy.float32),
            MCMCStateVariable("beta", 1.0, numpy.float32)
        ]

    @property
    def root_finder(self):
        return self._root_finder

    @property
    def cholesky(self):
        return self._cholesky
