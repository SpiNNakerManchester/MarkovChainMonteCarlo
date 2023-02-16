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
