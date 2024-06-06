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


class ARMAFloatModel(MCMCModel):
    """ MCMC Model for the ARMA problem, using single-point (float)
    """

    def __init__(
            self, parameters, jump_scale, root_finder=True, cholesky=True):
        """
        :param parameters:\
            array of coefficients of polynomials, plus mu and sigma
        :param jump_scale:\
            array of jump scale values for parameters
        """
        self._parameters = parameters
        self._jump_scale = jump_scale
        self._root_finder = root_finder
        self._cholesky = cholesky

    @overrides(MCMCModel.get_binary_name)
    def get_binary_name(self) -> str:
        return "arma.aplx"

    @overrides(MCMCModel.get_parameters)
    def get_parameters(self) -> List[MCMCParameter]:
        # Best here to convert the arrays into individual values
        return_params = []
        for i in range(len(self._jump_scale)):
            return_params.append(
                MCMCParameter(self._jump_scale[i], numpy.float32))

        return return_params

    @overrides(MCMCModel.get_state_variables)
    def get_state_variables(self) -> List[MCMCStateVariable]:
        # Best here to convert the arrays into individual values
        return_state_vars = []
        for i in range(len(self._parameters)):
            return_state_vars.append(
                MCMCStateVariable("param_"+str(i),
                                  self._parameters[i], numpy.float32))

        return return_state_vars

    @property
    def root_finder(self):
        return self._root_finder

    @property
    def cholesky(self):
        return self._cholesky
