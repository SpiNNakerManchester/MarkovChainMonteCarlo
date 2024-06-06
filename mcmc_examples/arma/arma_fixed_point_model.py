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

from spinn_front_end_common.interface.ds import DataType

from mcmc.mcmc_model import MCMCModel
from mcmc.mcmc_parameter import MCMCParameter
from mcmc.mcmc_state_variable import MCMCStateVariable
# pylint: disable=wrong-spelling-in-comment


class ARMAFixedPointModel(MCMCModel):
    """ MCMC Model for the ARMA problem, using fixed point
    """

    def __init__(
            self, parameters, p_jump_scale, q_jump_scale):
        """
        :param parameters: polynomial coefficients, mu, sigma
        :param p_jump_scale:\
            scaling of t transition distribution for MH jumps in p direction
        :param q_jump_scale:\
            scaling of t transition distribution for MH jumps in q direction
        """

        self._parameters = parameters
        self._p_jump_scale = p_jump_scale
        self._q_jump_scale = q_jump_scale

    @overrides(MCMCModel.get_binary_name)
    def get_binary_name(self) -> str:
        return "arma.aplx"

    @overrides(MCMCModel.get_parameters)
    def get_parameters(self) -> List[MCMCParameter]:
        return [
            MCMCParameter(self._parameters, DataType.S1615),  # array ?
            MCMCParameter(self._p_jump_scale, DataType.S1615),
            MCMCParameter(self._q_jump_scale, DataType.S1615)
        ]

    @overrides(MCMCModel.get_state_variables)
    def get_state_variables(self) -> List[MCMCStateVariable]:
        return [
            MCMCStateVariable("order_p", 10, DataType.S1615),  # check type
            MCMCStateVariable("order_q", 10, DataType.S1615)
        ]
