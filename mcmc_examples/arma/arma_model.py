from spinn_utilities.overrides import overrides

from mcmc.mcmc_model import MCMCModel
from mcmc.mcmc_parameter import MCMCParameter
from mcmc.mcmc_state_variable import MCMCStateVariable

import numpy


class ARMAModel(MCMCModel):
    """ MCMC Model for the ARMA problem, using double (float64)
    """

    def __init__(
            self, parameters, p_jump_scale, q_jump_scale):
        """
        :param parameters:\
            array of coefficients of polynomials, plus mu and sigma
        :param p_jump_scale:\
            array of scale values for polynomial coefficients for p
        :param q_jump_scale:\
            array of scale values for polynomial coefficients for q
        """

        self._parameters = parameters
        self._p_jump_scale = p_jump_scale
        self._q_jump_scale = q_jump_scale

    @overrides(MCMCModel.get_binary_name)
    def get_binary_name(self):
        return "arma.aplx"

    @overrides(MCMCModel.get_parameters)
    def get_parameters(self):
        # It's probably best here to convert the arrays into individual values?
        return_params = []
        for i in range(len(self._p_jump_scale)):
            return_params.append(
                MCMCParameter(self._p_jump_scale[i], numpy.float64))
        for i in range(len(self._q_jump_scale)):
            return_params.append(
                MCMCParameter(self._q_jump_scale[i], numpy.float64))

        return return_params

    @overrides(MCMCModel.get_state_variables)
    def get_state_variables(self):
        # It's probably best here to convert the arrays into individual values?
        return_state_vars = []
        for i in range(len(self._parameters)):
            return_state_vars.append(
                MCMCStateVariable("param_"+str(i),
                                  self._parameters[i], numpy.float64))

        return return_state_vars
