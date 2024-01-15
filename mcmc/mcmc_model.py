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
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

from typing import List
from spinn_utilities.abstract_base import AbstractBase
from spinn_utilities.abstract_base import abstractmethod
from mcmc.mcmc_parameter import MCMCParameter
from mcmc.mcmc_state_variable import MCMCStateVariable


class MCMCModel(object, metaclass=AbstractBase):

    @abstractmethod
    def get_parameters(self) -> List[MCMCParameter]:
        """ Get the parameters of the model

        :rtype: list of :py:class:`mcmc.mcmc_parameter.MCMCParameter`
        """
        raise NotImplementedError

    @abstractmethod
    def get_state_variables(self) -> List[MCMCStateVariable]:
        """ Get the state variables of the model

        :rtype: list of :py:class:`mcmc.mcmc_state_variable.MCMCStateVariable`
        """
        raise NotImplementedError

    @abstractmethod
    def get_binary_name(self) -> str:
        """ Get the name of the binary compiled with this model

        :rtype: str
        """
        raise NotImplementedError
