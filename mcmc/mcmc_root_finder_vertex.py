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

from pacman.model.graphs.machine import MachineVertex
from pacman.model.resources import ConstantSDRAM
from spinn_utilities.overrides import overrides

from pacman.model.placements import Placement

from spinnman.model.enums import ExecutableType

from spinn_front_end_common.abstract_models.abstract_has_associated_binary \
    import AbstractHasAssociatedBinary
from spinn_front_end_common.abstract_models\
    .abstract_generates_data_specification \
    import AbstractGeneratesDataSpecification
from spinn_front_end_common.interface.ds import DataSpecificationGenerator

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
    def sdram_required(self) -> ConstantSDRAM:
        return ConstantSDRAM(self._n_parameter_bytes)

    @overrides(AbstractHasAssociatedBinary.get_binary_file_name)
    def get_binary_file_name(self) -> str:
        return "mcmc_root_finder.aplx"

    @overrides(AbstractHasAssociatedBinary.get_binary_start_type)
    def get_binary_start_type(self) -> ExecutableType:
        return ExecutableType.SYNC

    @overrides(
        AbstractGeneratesDataSpecification.generate_data_specification)
    def generate_data_specification(
            self, spec: DataSpecificationGenerator, placement: Placement):

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
