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

from enum import Enum
from typing import Sequence
import numpy

from spinn_utilities.overrides import overrides

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
from spinn_front_end_common.interface.buffer_management.buffer_models\
    .abstract_receive_buffers_to_host import AbstractReceiveBuffersToHost
from spinn_front_end_common.interface.ds import (
    DataSpecificationGenerator, DataType)
from spinn_front_end_common.utilities import helpful_functions
from spinn_front_end_common.interface.buffer_management \
    import recording_utilities


class MCMCRegions(Enum):
    """ Regions in the MCMC Data
    """
    RECORDING = 0
    PARAMETERS = 1
    MODEL_PARAMS = 2
    MODEL_STATE = 3


class MCMCVertex(
        MachineVertex, AbstractHasAssociatedBinary,
        AbstractGeneratesDataSpecification, AbstractReceiveBuffersToHost):
    """ A vertex that runs the MCMC algorithm
    """

    def __init__(self, coordinator, model,
                 parameter_partition_name="MCMCParameter",
                 result_partition_name="MCMCResultAck",
                 cholesky_partition_name="MCMCCholeskyParameter",
                 cholesky_result_partition_name="MCMCCholeskyResultAck"):
        """

        :param coordinator: The coordinator vertex
        :param model: The model being simulated
        """

        MachineVertex.__init__(self, label="MCMC Node")
        self._coordinator = coordinator
        self._model = model
        self._parameter_partition_name = parameter_partition_name
        self._result_partition_name = result_partition_name
        self._cholesky_partition_name = cholesky_partition_name
        self._cholesky_result_partition_name = cholesky_result_partition_name

        self._coordinator.register_processor(self)

        state = self._get_model_state_array()
        self._recording_size = self._coordinator.n_samples * len(state) * 4

        params = self._get_model_parameters_array()

        # The number of bytes for the parameters
        # (11 * uint32) + (5 * seed array) + (1 * d.o.f.)
        self._n_parameter_bytes = 0
        if (model.get_parameters()[0].data_type is numpy.float64):
            self._n_parameter_bytes = (11 * 4) + (5 * 4) + (1 * 8)
        else:
            self._n_parameter_bytes = (11 * 4) + (5 * 4) + (1 * 4)

        self._sdram_usage = (
            self._n_parameter_bytes + self._recording_size +
            recording_utilities.get_recording_header_size(1) +
            (len(params) * 4)
            )

        self._data_receiver = dict()
        self._cholesky_data_receiver = dict()

    def _get_model_parameters_array(self):
        """
        Convert the models parameters into a numpy array
        """
        parameters = self._model.get_parameters()
        numpy_format = list()
        numpy_values = list()
        for i, param in enumerate(parameters):
            if (param.data_type is numpy.float64):
                numpy_format.append((f"f{{{i}}}", param.data_type))
                numpy_values.append(param.value)
            elif (param.data_type is numpy.float32):
                numpy_format.append((f"f{{{i}}}", param.data_type))
                numpy_values.append(param.value)
            elif (param.data_type is DataType.S1615):
                numpy_format.append((f"f{{{i}}}", numpy.uint32))
                numpy_values.append(
                    int(param.value * float(DataType.S1615.scale)))
            else:
                # throw exception for unknown data type
                raise TypeError(
                    "Error: unsupported data type used for model parameters")
        return numpy.array(
            [tuple(numpy_values)], dtype=numpy_format).view("uint32")

    def _get_model_state_array(self):
        """
        Convert the models state variables into a numpy array
        """
        state = self._model.get_state_variables()
        numpy_format = list()
        numpy_values = list()
        for i, param in enumerate(state):
            if (param.data_type is numpy.float64):
                numpy_format.append((f"f{{{i}}}", param.data_type))
                numpy_values.append(param.initial_value)
            elif (param.data_type is numpy.float32):
                numpy_format.append((f"f{{{i}}}", param.data_type))
                numpy_values.append(param.initial_value)
            elif (param.data_type is numpy.uint32):
                numpy_format.append((f"f{{{i}}}", param.data_type))
                numpy_values.append(param.initial_value)
            elif (param.data_type is DataType.S1615):
                numpy_format.append((f"f{{{i}}}", numpy.uint32))
                numpy_values.append(
                    int(param.initial_value * float(DataType.S1615.scale)))
            else:
                # throw exception for unknown data type
                raise TypeError(
                    "Error: unsupported data type for model state params")
        return numpy.array(
            [tuple(numpy_values)], dtype=numpy_format).view("uint32")

    def get_result_key(self, placement):
        """
        Get the key for the given placement if it is the first on this X Y

        param Placement placement: the placement to get the result key for
        """
        if self._is_receiver_placement(placement):
            routing_info = FecDataView.get_routing_infos()
            key = routing_info.get_first_key_from_pre_vertex(
                placement.vertex, self._result_partition_name)
            return key
        return 0

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

    def get_cholesky_result_key(self, placement):
        """
        Get the key for the given placement if it is the first on this X Y

        param Placement placement: the placement to get the result key for
        """
        if self._is_cholesky_receiver_placement(placement):
            routing_info = FecDataView.get_routing_infos()
            key = routing_info.get_first_key_from_pre_vertex(
                placement.vertex, self._cholesky_result_partition_name)
            return key
        return 0

    def _is_cholesky_receiver_placement(self, placement):
        """
        Checks if this is the first placement seen with this X, Y
        """
        x = placement.x
        y = placement.y
        if (x, y) not in self._cholesky_data_receiver:
            self._cholesky_data_receiver[x, y] = placement.p
            return True
        return self._cholesky_data_receiver[(x, y)] == placement.p

    @property
    def coordinator(self):
        """
        Returns the coordinator vertex passed intoi the init
        """
        return self._coordinator

    @property
    def parameter_partition_name(self):
        """
        The parameter_partition_name passed into the init
        """
        return self._parameter_partition_name

    @property
    def result_partition_name(self):
        """
        The result_partition_name passed into the init
        """
        return self._result_partition_name

    @property
    def cholesky_partition_name(self):
        """
        The cholesky_partition_name passed into the init
        """
        return self._cholesky_partition_name

    @property
    def cholesky_result_partition_name(self):
        """
        The cholesky_result_partition_name passed into the init
        """
        return self._cholesky_result_partition_name

    @property
    @overrides(MachineVertex.sdram_required)
    def sdram_required(self) -> ConstantSDRAM:
        return ConstantSDRAM(self._sdram_usage)

    @overrides(AbstractHasAssociatedBinary.get_binary_file_name)
    def get_binary_file_name(self) -> str:
        return self._model.get_binary_name()

    @overrides(AbstractHasAssociatedBinary.get_binary_start_type)
    def get_binary_start_type(self) -> ExecutableType:
        return ExecutableType.SYNC

    @overrides(
        AbstractGeneratesDataSpecification.generate_data_specification)
    def generate_data_specification(
            self, spec: DataSpecificationGenerator, placement: Placement):

        routing_info = FecDataView.get_routing_infos()
        # Reserve and write the recording regions
        spec.reserve_memory_region(
            MCMCRegions.RECORDING.value,
            recording_utilities.get_recording_header_size(1))
        spec.switch_write_focus(MCMCRegions.RECORDING.value)
        spec.write_array(recording_utilities.get_recording_header_array(
            [self._recording_size]))

        # Reserve and write the parameters region
        spec.reserve_memory_region(
            MCMCRegions.PARAMETERS.value, self._n_parameter_bytes)
        spec.switch_write_focus(MCMCRegions.PARAMETERS.value)

        # Write the burn-in
        spec.write_value(self._coordinator.burn_in, data_type=DataType.UINT32)

        # Write the thinning value
        spec.write_value(self._coordinator.thinning, data_type=DataType.UINT32)

        # Write the number of samples
        spec.write_value(
            self._coordinator.n_samples, data_type=DataType.UINT32)

        # Write the number of data points
        spec.write_value(self._coordinator.n_data_points)

        # Write the data window size
        spec.write_value(self._coordinator.get_data_window_size(placement))

        # Write the sequence mask
        spec.write_value(self._coordinator.get_sequence_mask(
            placement, routing_info))

        # Write the acknowledge key
        spec.write_value(self._coordinator.get_acknowledge_key(
            placement, routing_info))

        # Write the data tag
        spec.write_value(self._coordinator.data_tag)

        # Write the timer value
        spec.write_value(self._coordinator.acknowledge_timer)

        # Write the (first) key for sending parameter data, if needed
        if (self._model.root_finder):
            routing_info_rf = routing_info.get_routing_info_from_pre_vertex(
                self, self._parameter_partition_name)
            assert routing_info_rf is not None
            spec.write_value(routing_info_rf.key, data_type=DataType.UINT32)
        else:
            spec.write_value(0, data_type=DataType.UINT32)

        if (self._model.cholesky):
            routing_info_ch = routing_info.get_routing_info_from_pre_vertex(
                self, self._cholesky_partition_name)
            assert routing_info_ch is not None
            spec.write_value(routing_info_ch.key, data_type=DataType.UINT32)
        else:
            spec.write_value(0, data_type=DataType.UINT32)

        # Write the seed = 5 32-bit random numbers
        if (self._model.get_parameters()[0].data_type is DataType.S1615):
            seed = [int(x * DataType.S1615.scale)
                    for x in self._coordinator.seed]
            spec.write_array(seed)
        else:
            spec.write_array(self._coordinator.seed)

        # Write the degrees of freedom
        if (self._model.get_parameters()[0].data_type is numpy.float64):
            spec.write_value(self._coordinator.degrees_of_freedom,
                             data_type=DataType.FLOAT_64)
        elif (self._model.get_parameters()[0].data_type is numpy.float32):
            spec.write_value(self._coordinator.degrees_of_freedom,
                             data_type=DataType.FLOAT_32)
        elif (self._model.get_parameters()[0].data_type is DataType.S1615):
            degrees_of_freedom = int(
                self._coordinator.degrees_of_freedom * float(
                    DataType.S1615.scale))
            spec.write_value(degrees_of_freedom, data_type=DataType.UINT32)

        # Reserve and write the model parameters
        params = self._get_model_parameters_array()
        spec.reserve_memory_region(
            MCMCRegions.MODEL_PARAMS.value, len(params) * 4)
        spec.switch_write_focus(MCMCRegions.MODEL_PARAMS.value)
        spec.write_array(params)

        # Reserve and write the model state
        state = self._get_model_state_array()
        spec.reserve_memory_region(
            MCMCRegions.MODEL_STATE.value, len(state) * 4)
        spec.switch_write_focus(MCMCRegions.MODEL_STATE.value)
        spec.write_array(state)

        # End the specification
        spec.end_specification()

    def read_samples(self, buffer_manager, placement):
        """ Read back the samples
        """

        # Read the data recorded
        data, _ = buffer_manager.get_data_by_placement(placement, 0)

        numpy_format = list()
        output_format = list()
        for var in self._model.get_state_variables():
            if (var.data_type is DataType.S1615):
                numpy_format.append((var.name, numpy.int32))
                output_format.append((var.name, numpy.float32))
            else:
                numpy_format.append((var.name, var.data_type))

        # Convert the data into an array of state variables
        if (self._model.get_parameters()[0].data_type is DataType.S1615):
            data_view = numpy.array(data, dtype=numpy.uint8).view(numpy_format)
            convert = numpy.zeros_like(
                data_view, dtype=numpy.float64).view(output_format)

            for i in range(data_view.size):
                for j in range(len(numpy_format)):
                    convert[i][j] = float(
                        data_view[i][j]) / float(DataType.S1615.scale)

            return convert
        else:
            return numpy.array(data, dtype="uint8").view(numpy_format)

    @overrides(AbstractReceiveBuffersToHost.get_recorded_region_ids)
    def get_recorded_region_ids(self) -> Sequence[int]:
        return [0]

    @overrides(AbstractReceiveBuffersToHost.get_recording_region_base_address)
    def get_recording_region_base_address(self, placement):
        return helpful_functions.locate_memory_region_for_placement(
            placement, MCMCRegions.RECORDING.value)
