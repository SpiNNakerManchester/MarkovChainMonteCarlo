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

from typing import List
import random
import numpy

from spinn_utilities.overrides import overrides
from spinn_utilities.progress_bar import ProgressBar

from spinnman.model.enums import ExecutableType

from pacman.model.graphs.machine import MachineVertex
from pacman.model.placements import Placement
from pacman.model.resources import ConstantSDRAM

from spinn_front_end_common.abstract_models.abstract_has_associated_binary \
    import AbstractHasAssociatedBinary
from spinn_front_end_common.abstract_models\
    .abstract_generates_data_specification \
    import AbstractGeneratesDataSpecification
from spinn_front_end_common.data import FecDataView
from spinn_front_end_common.interface.ds import (
    DataSpecificationGenerator, DataType)


class MCMCCoordinatorVertex(
        MachineVertex, AbstractHasAssociatedBinary,
        AbstractGeneratesDataSpecification):
    """ A vertex that runs the MCMC algorithm
    """

    # The number of bytes for the parameters
    # (in mcmc_coordinator.c, 6*4 not including data)
    _N_PARAMETER_BYTES = 24

    # The data type of the data count
    _DATA_COUNT_TYPE = DataType.UINT32

    # The data type of the keys
    _KEY_ELEMENT_TYPE = DataType.UINT32

    def __init__(
            self, model, data, n_samples, burn_in, thinning,
            degrees_of_freedom, seed=None,
            send_timer=1000, receive_timer=1000, window_size=1024,
            n_sequences=2048, data_partition_name="MCMCData",
            acknowledge_partition_name="MCMCDataAck", data_tag=1):
        """

        :param model: The model being simulated
        :param data: The data to sample
        :param n_samples: The number of samples to generate
        :param burn_in:\
            no of MCMC transitions to reach apparent equilibrium before\
            generating inference samples
        :param thinning:\
            sampling rate i.e. 5 = 1 sample for 5 generated steps
        :param degrees_of_freedom:\
            The number of degrees of freedom to jump around with
        :param seed: The random seed to use
        """

        MachineVertex.__init__(self, label="MCMC Coordinator Node")
        self._model = model
        self._data = data
        self._n_samples = n_samples
        self._burn_in = burn_in
        self._thinning = thinning
        self._degrees_of_freedom = degrees_of_freedom
        self._seed = seed
        self._send_timer = send_timer
        self._receive_timer = receive_timer
        self._window_size = window_size
        self._n_sequences = n_sequences
        self._data_partition_name = data_partition_name
        self._acknowledge_partition_name = acknowledge_partition_name
        self._data_tag = data_tag

        # The data type of each data element
        if (self._model.get_parameters()[0].data_type is numpy.float64):
            self._data_element_type = DataType.FLOAT_64
        elif (self._model.get_parameters()[0].data_type is numpy.float32):
            self._data_element_type = DataType.FLOAT_32
        elif (self._model.get_parameters()[0].data_type is DataType.S1615):
            self._data_element_type = DataType.S1615

        # The numpy data type of each data element
        if (self._model.get_parameters()[0].data_type is numpy.float64):
            self._numpy_data_element_type = numpy.float64
        elif (self._model.get_parameters()[0].data_type is numpy.float32):
            self._numpy_data_element_type = numpy.float32
        elif (self._model.get_parameters()[0].data_type is DataType.S1615):
            self._numpy_data_element_type = numpy.uint32

        self._data_size = (
            (len(self._data) * self._data_element_type.size) +
            self._DATA_COUNT_TYPE.size
        )
        self._sdram_usage = (
            self._N_PARAMETER_BYTES + self._data_size
        )

        self._mcmc_vertices = list()
        self._mcmc_placements = list()
        self._data_receiver = dict()

    def register_processor(self, mcmc_vertex):
        """
        Records the vertex to its list of vertices
        """
        self._mcmc_vertices.append(mcmc_vertex)

    @property
    def n_samples(self):
        """
        n_samples as passed into the init.
        """
        return self._n_samples

    @property
    def burn_in(self):
        """
        burn_in as passed inot the init
        """
        return self._burn_in

    @property
    def thinning(self):
        """
        thinning as passed into the init
        """
        return self._thinning

    @property
    def degrees_of_freedom(self):
        """
        degrees_of_freedom as passed into the init
        """
        return self._degrees_of_freedom

    @property
    def n_data_points(self):
        """
        The length of the data passed inot the init
        """
        return len(self._data)

    def _is_receiver_placement(self, placement):
        """
        Checks if this is the first placement seen with this X, Y
        """
        x = placement.x
        y = placement.y
        if (x, y) not in self._data_receiver:
            self._data_receiver[x, y] = placement.p
            return True
        return self._data_receiver[(x, y)] == placement.p

    def get_data_window_size(self, placement):
        """
        Gets the window size passed into the init if this is the first \
        placement seen with this X, Y otherwise 0
        """
        if self._is_receiver_placement(placement):
            return self._window_size
        return 0

    def get_sequence_mask(self, placement, routing_info):
        """
        Gets the mask for the data_partition_name vertex

        Only if this is the first placement seen with this X, Y otherwise 0
        """
        if self._is_receiver_placement(placement):
            mask = routing_info.get_routing_info_from_pre_vertex(
                self, self._data_partition_name).mask
            return ~mask & 0xFFFFFFFF
        return 0

    def get_acknowledge_key(self, placement, routing_info):
        """
        Gets the key for the acknowledge_partition_name vertex

        Only if this is the first placement seen with this X, Y otherwise 0
        """
        if self._is_receiver_placement(placement):
            key = routing_info.get_first_key_from_pre_vertex(
                placement.vertex, self._acknowledge_partition_name)
            return key
        return 0

    @property
    def data_tag(self):
        """
        Gets the data_tag value passed into the init
        """
        return self._data_tag

    @property
    def acknowledge_timer(self):
        """
        Gets the receive_timer passed inot the init
        """
        return self._receive_timer

    @property
    def seed(self):
        """
        Get a consistent random seed

        Uses the one passed into the init or creates one on the first call
        """
        if self._seed is None:
            return [random.randint(0, 0xFFFFFFFF) for _ in range(5)]
        return self._seed

    @property
    def data_partition_name(self):
        """
        The data_partition_name passed into the init.
        """
        return self._data_partition_name

    @property
    def acknowledge_partition_name(self):
        """
        acknowledge_partition_name passed inot the init.
        """
        return self._acknowledge_partition_name

    @property
    @overrides(MachineVertex.sdram_required)
    def sdram_required(self) -> ConstantSDRAM:
        sdram = self._N_PARAMETER_BYTES + self._data_size
        sdram += len(self._mcmc_vertices) * self._KEY_ELEMENT_TYPE.size
        return ConstantSDRAM(sdram)

    @overrides(AbstractHasAssociatedBinary.get_binary_file_name)
    def get_binary_file_name(self) -> str:
        return "mcmc_coordinator.aplx"

    @overrides(AbstractHasAssociatedBinary.get_binary_start_type)
    def get_binary_start_type(self) -> ExecutableType:
        return ExecutableType.SYNC

    @overrides(
        AbstractGeneratesDataSpecification.generate_data_specification)
    def generate_data_specification(
            self, spec: DataSpecificationGenerator, placement: Placement):
        routing_info = FecDataView.get_routing_infos()

        # Reserve and write the parameters region
        region_size = self._N_PARAMETER_BYTES + self._data_size
        region_size += len(self._mcmc_vertices) * self._KEY_ELEMENT_TYPE.size
        spec.reserve_memory_region(0, region_size)
        spec.switch_write_focus(0)

        # Get the placement of the vertices and find out how many chips
        # are needed
        keys: List[int] = list()
        for vertex in self._mcmc_vertices:
            mcmc_placement = FecDataView.get_placement_of_vertex(vertex)
            self._mcmc_placements.append(mcmc_placement)
            if self._is_receiver_placement(mcmc_placement):
                key = routing_info.get_first_key_from_pre_vertex(
                    vertex, self._acknowledge_partition_name)
                if key is not None:
                    keys.append(key)
        keys.sort()

        # Write the data size in words
        spec.write_value(
            int(len(self._data) * (float(self._data_element_type.size) / 4.0)),
            data_type=self._DATA_COUNT_TYPE)

        # Write the number of chips
        spec.write_value(len(keys), data_type=DataType.UINT32)

        # Write the key
        vertex_routing_info = routing_info.get_routing_info_from_pre_vertex(
            self, self._data_partition_name)
        assert vertex_routing_info is not None
        spec.write_value(vertex_routing_info.key, data_type=DataType.UINT32)

        # Write the window size
        spec.write_value(self._window_size, data_type=DataType.UINT32)

        # Write the sequence mask
        spec.write_value(
            ~vertex_routing_info.mask & 0xFFFFFFFF, data_type=DataType.UINT32)

        # Write the timer
        spec.write_value(self._send_timer, data_type=DataType.UINT32)

        # Write the data - Arrays must be 32-bit values, so convert
        if (self._model.get_parameters()[0].data_type is DataType.S1615):
            data_convert = [int(x * float(DataType.S1615.scale))
                            for x in self._data]
            data = numpy.array(data_convert,
                               dtype=self._numpy_data_element_type)
        else:
            data = numpy.array(self._data,
                               dtype=self._numpy_data_element_type)

        spec.write_array(data.view(numpy.uint32))

        # Write the keys
        spec.write_array(keys)

        # End the specification
        spec.end_specification()

    @overrides(MachineVertex.get_n_keys_for_partition)
    def get_n_keys_for_partition(self, partition_id: str) -> int:
        return self._n_sequences

    def read_samples(self, buffer_manager):
        """ Read back the samples
        """
        progress = ProgressBar(len(self._mcmc_placements), "Reading results")
        samples = list()
        for placement in self._mcmc_placements:
            # Read the data recorded
            sample = placement.vertex.read_samples(buffer_manager, placement)
            if sample is not None:
                samples.append(sample)
            progress.update()
        progress.end()

        # Merge all the arrays
        return numpy.hstack(samples)
