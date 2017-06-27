from pacman.model.graphs.machine import MachineVertex
from pacman.model.resources.resource_container import ResourceContainer
from pacman.model.resources.dtcm_resource import DTCMResource
from pacman.model.resources.sdram_resource import SDRAMResource
from pacman.model.resources.cpu_cycles_per_tick_resource \
    import CPUCyclesPerTickResource
from pacman.model.decorators.overrides import overrides
from pacman.executor.injection_decorator import inject_items

from data_specification.enums.data_type import DataType

from spinn_front_end_common.abstract_models.abstract_has_associated_binary \
    import AbstractHasAssociatedBinary
from spinn_front_end_common.abstract_models\
    .abstract_generates_data_specification \
    import AbstractGeneratesDataSpecification
from spinn_front_end_common.interface.buffer_management.buffer_models\
    .abstract_receive_buffers_to_host import AbstractReceiveBuffersToHost
from spinn_front_end_common.utilities import constants, helpful_functions
from spinn_front_end_common.interface.buffer_management \
    import recording_utilities
from spinn_front_end_common.interface.simulation import simulation_utilities
from spinn_front_end_common.utilities.utility_objs.executable_start_type \
    import ExecutableStartType

from enum import Enum
import numpy
import random


class MCMCRegions(Enum):
    """ Regions in the MCMC Data
    """
    RECORDING = 0
    PARAMETERS = 1


class MCMCVertex(
        MachineVertex, AbstractHasAssociatedBinary,
        AbstractGeneratesDataSpecification, AbstractReceiveBuffersToHost):
    """ A vertex that runs the MCMC algorithm
    """

    # The number of bytes for the parameters
    _N_PARAMETER_BYTES = 116

    # The size of each sample in bytes (2 doubles)
    _SAMPLE_SIZE = (8 * 2)

    def __init__(self, coordinator):
        """

        :param coordinator: The coordinator vertex
        """

        MachineVertex.__init__(self, label="MCMC Node", constraints=None)
        self._coordinator = coordinator
        self._coordinator.register_processor(self)
        self._recording_size = self._coordinator.n_samples * self._SAMPLE_SIZE
        self._sdram_usage = (
            self._N_PARAMETER_BYTES + self._recording_size +
            recording_utilities.get_recording_header_size(1)
            )


    @property
    @overrides(MachineVertex.resources_required)
    def resources_required(self):

        resources = ResourceContainer(
            dtcm=DTCMResource(0),
            sdram=SDRAMResource(self._sdram_usage),
            cpu_cycles=CPUCyclesPerTickResource(0),
            iptags=[], reverse_iptags=[])
        return resources

    @overrides(AbstractHasAssociatedBinary.get_binary_file_name)
    def get_binary_file_name(self):
        return "mcmc.aplx"

    @overrides(AbstractHasAssociatedBinary.get_binary_start_type)
    def get_binary_start_type(self):
        return ExecutableStartType.SYNC

    @inject_items({
        "routing_info": "MemoryRoutingInfos",
        "tags": "MemoryTags",
    })
    @overrides(
        AbstractGeneratesDataSpecification.generate_data_specification,
        additional_arguments=["routing_info", "tags"])
    def generate_data_specification(
            self, spec, placement, routing_info, tags):

        # Reserve and write the recording regions
        spec.reserve_memory_region(
            MCMCRegions.RECORDING.value, recording_utilities.get_recording_header_size(1))
        spec.switch_write_focus(MCMCRegions.RECORDING.value)
        ip_tags = tags.get_ip_tags_for_vertex(self) or []
        spec.write_array(recording_utilities.get_recording_header_array(
            [self._recording_size], ip_tags=ip_tags))

        # Reserve and write the parameters region
        spec.reserve_memory_region(
            MCMCRegions.PARAMETERS.value, self._N_PARAMETER_BYTES)
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

        # Write a padding word
        spec.write_value(0, data_type=DataType.UINT32)

        # Write the alpha and beta jump scales, ranges and degrees of freedom
        values = numpy.array(
            [self._coordinator.alpha_jump_scale,
             self._coordinator.beta_jump_scale,
             self._coordinator.alpha_min, self._coordinator.alpha_max,
             self._coordinator.beta_min, self._coordinator.beta_max,
             self._coordinator.degrees_of_freedom],
            dtype=numpy.float64).view(numpy.uint32)
        spec.write_array(values)

        # Write the seed = 5 32-bit random numbers
        spec.write_array(self._coordinator.seed)

        # End the specification
        spec.end_specification()

    def read_samples(self, buffer_manager, placement):
        """ Read back the samples
        """

        # Read the data recorded
        data_values, missing = buffer_manager.get_data_for_vertex(placement, 0)
        data = data_values.read_all()

        # Convert the data into an array of 2-doubles
        return numpy.array(
            data, dtype=numpy.uint8).view([
                ("alpha", numpy.float64),
                ("beta", numpy.float64)])

    def get_minimum_buffer_sdram_usage(self):
        return 1024

    def get_n_timesteps_in_buffer_space(self, buffer_space, machine_time_step):
        return recording_utilities.get_n_timesteps_in_buffer_space(
            buffer_space, 4)

    def get_recorded_region_ids(self):
        return [0]

    def get_recording_region_base_address(self, txrx, placement):
        return helpful_functions.locate_memory_region_for_placement(
            placement, MCMCRegions.RECORDING.value, txrx)
