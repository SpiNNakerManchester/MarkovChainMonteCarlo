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
