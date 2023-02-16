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

from pacman.model.graphs.machine import MachineVertex
from pacman.model.resources import ConstantSDRAM
from spinn_utilities.overrides import overrides

from spinn_front_end_common.abstract_models.abstract_has_associated_binary \
    import AbstractHasAssociatedBinary
from spinn_front_end_common.abstract_models\
    .abstract_generates_data_specification \
    import AbstractGeneratesDataSpecification
from spinn_front_end_common.utilities.utility_objs.executable_type \
    import ExecutableType

from enum import Enum


class MCMCRootFinderRegions(Enum):
    """ Regions in the MCMCRootFinder Data
    """
    PARAMETERS = 0


class MCMCRootFinderVertex(
        MachineVertex, AbstractHasAssociatedBinary,
        AbstractGeneratesDataSpecification):
    """ A vertex that runs the (MCMC) root finder algorithm
    """

    def __init__(
            self, vertex, model):
        """
        :param vertex: The MCMC vertex associated with this vertex
        :param model: The model being simulated
        """

        MachineVertex.__init__(self, label="MCMC RF Node")
        self._vertex = vertex
        self._model = model

        vertex.coordinator.register_processor(self)

        # Other parameters need to be added here
        self._n_other_params = 1

        self._n_parameter_bytes = self._n_other_params * 4

    @property
    @overrides(MachineVertex.sdram_required)
    def sdram_required(self):
        return ConstantSDRAM(self._n_parameter_bytes)

    @overrides(AbstractHasAssociatedBinary.get_binary_file_name)
    def get_binary_file_name(self):
        return "mcmc_root_finder.aplx"

    @overrides(AbstractHasAssociatedBinary.get_binary_start_type)
    def get_binary_start_type(self):
        return ExecutableType.SYNC

    @overrides(
        AbstractGeneratesDataSpecification.generate_data_specification)
    def generate_data_specification(self, spec, placement):

        # Reserve and write the parameters region
        spec.reserve_memory_region(
            MCMCRootFinderRegions.PARAMETERS.value, self._n_parameter_bytes)
        spec.switch_write_focus(MCMCRootFinderRegions.PARAMETERS.value)

        # Write the acknowledge key
        spec.write_value(self._vertex.get_result_key(placement))

        # End the specification
        spec.end_specification()

    def read_samples(self, buffer_manager, placement):
        # pylint: disable=unused-argument
        """ Read back the samples (dummy call)
        """
#        print 'There are no samples to read back on a root finder vertex'
        return None
