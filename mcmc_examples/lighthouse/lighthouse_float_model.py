# Copyright (c) 2016 The University of Manchester
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
    def get_binary_name(self):
        return "lighthouse.aplx"

    @overrides(MCMCModel.get_parameters)
    def get_parameters(self):
        return [
            MCMCParameter(self._alpha_jump_scale, numpy.float32),
            MCMCParameter(self._beta_jump_scale, numpy.float32),
            MCMCParameter(self._alpha_min, numpy.float32),
            MCMCParameter(self._alpha_max, numpy.float32),
            MCMCParameter(self._beta_min, numpy.float32),
            MCMCParameter(self._beta_max, numpy.float32)
        ]

    @overrides(MCMCModel.get_state_variables)
    def get_state_variables(self):
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
