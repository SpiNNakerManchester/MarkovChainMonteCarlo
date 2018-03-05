class MCMCParameter(object):
    """ A Parameter of the MCMC simulation
    """

    def __init__(self, value, data_type):
        """

        :param value: The value of the parameter
        :param data_type: The numpy data type of the parameter
        """
        self._value = value
        self._data_type = data_type

    @property
    def value(self):
        return self._value

    @property
    def data_type(self):
        return self._data_type
