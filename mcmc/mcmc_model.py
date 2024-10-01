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
from spinn_utilities.abstract_base import AbstractBase
from spinn_utilities.abstract_base import abstractmethod
from mcmc.mcmc_parameter import MCMCParameter
from mcmc.mcmc_state_variable import MCMCStateVariable


class MCMCModel(object, metaclass=AbstractBase):
    """
    Base class frr Markov chain Monte Carlo models
    """

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
