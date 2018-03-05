class MCMCStateVariable(object):
    """ A State variable of the MCMC simulation
    """

    def __init__(self, name, initial_value, data_type):
        """

        :param name: The name of the variable
        :param initial_value: The initial value of the variable
        :param data_type: The numpy data type of the variable
        """
        self._name = name
        self._initial_value = initial_value
        self._data_type = data_type

    @property
    def name(self):
        return self._name

    @property
    def initial_value(self):
        return self._initial_value

    @property
    def data_type(self):
        return self._data_type
