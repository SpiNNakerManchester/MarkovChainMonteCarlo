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
