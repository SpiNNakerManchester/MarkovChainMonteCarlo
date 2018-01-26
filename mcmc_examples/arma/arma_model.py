from spinn_utilities.overrides import overrides

from mcmc.mcmc_model import MCMCModel
from mcmc.mcmc_parameter import MCMCParameter
from mcmc.mcmc_state_variable import MCMCStateVariable

import numpy


class ARMAModel(MCMCModel):
    """ MCMC Model for the ARMA problem, using single-point (float)
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

        return [
            MCMCParameter(self._parameters, numpy.float64), # it's an array?
            MCMCParameter(self._p_jump_scale, numpy.float64),
            MCMCParameter(self._q_jump_scale, numpy.float64)
        ]

    @overrides(MCMCModel.get_state_variables)
    def get_state_variables(self):
        return [
            MCMCStateVariable("order_p", 9, numpy.uint32), # edit as appropriate
            MCMCStateVariable("order_q", 9, numpy.uint32)
        ]
