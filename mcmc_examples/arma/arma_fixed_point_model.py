# Copyright (c) 2016-2021 The University of Manchester
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

from spinn_utilities.overrides import overrides

from data_specification.enums.data_type import DataType

from mcmc.mcmc_model import MCMCModel
from mcmc.mcmc_parameter import MCMCParameter
from mcmc.mcmc_state_variable import MCMCStateVariable


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
    def get_binary_name(self):
        return "arma.aplx"

    @overrides(MCMCModel.get_parameters)
    def get_parameters(self):
        return [
            MCMCParameter(self._parameters, DataType.S1615),  # array ?
            MCMCParameter(self._p_jump_scale, DataType.S1615),
            MCMCParameter(self._q_jump_scale, DataType.S1615)
        ]

    @overrides(MCMCModel.get_state_variables)
    def get_state_variables(self):
        return [
            MCMCStateVariable("order_p", 10, DataType.S1615),  # check type
            MCMCStateVariable("order_q", 10, DataType.S1615)
        ]
