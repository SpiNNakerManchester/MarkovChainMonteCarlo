from spinn_utilities.abstract_base import AbstractBase
from spinn_utilities.abstract_base import abstractmethod


class MCMCModel(object, metaclass=AbstractBase):

    @abstractmethod
    def get_parameters(self):
        """ Get the parameters of the model

        :rtype: list of :py:class:`mcmc.mcmc_parameter.MCMCParameter`
        """

    @abstractmethod
    def get_state_variables(self):
        """ Get the state variables of the model

        :rtype: list of :py:class:`mcmc.mcmc_state_variable.MCMCStateVariable`
        """

    @abstractmethod
    def get_binary_name(self):
        """ Get the name of the binary compiled with this model

        :rtype: str
        """
