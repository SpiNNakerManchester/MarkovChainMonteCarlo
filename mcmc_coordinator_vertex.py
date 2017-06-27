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
from spinn_front_end_common.abstract_models\
    .abstract_provides_n_keys_for_partition \
    import AbstractProvidesNKeysForPartition
from spinn_front_end_common.utilities.utility_objs.executable_start_type \
    import ExecutableStartType

from spinn_machine.utilities.progress_bar import ProgressBar

import numpy
import random


class MCMCCoordinatorVertex(
        MachineVertex, AbstractHasAssociatedBinary,
        AbstractGeneratesDataSpecification,
        AbstractProvidesNKeysForPartition):
    """ A vertex that runs the MCMC algorithm
    """

    # The number of bytes for the parameters
    _N_PARAMETER_BYTES = 24

    # The data type of the data count
    _DATA_COUNT_TYPE = DataType.UINT32

    # The data type of each data element
    _DATA_ELEMENT_TYPE = DataType.FLOAT_64

    # The numpy data type of each data element
    _NUMPY_DATA_ELEMENT_TYPE = numpy.double

    # The data type of the keys
    _KEY_ELEMENT_TYPE = DataType.UINT32

    def __init__(
            self, data, n_samples, burn_in, thinning,
            alpha_jump_scale, alpha_min, alpha_max,
            beta_jump_scale, beta_min, beta_max,
            degrees_of_freedom, seed=None,
            send_timer=1000, receive_timer=1000, window_size=1024,
            n_sequences=2048, data_partition_name="MCMCData",
            acknowledge_partition_name="MCMCDataAck", data_tag=1):
        """

        :param data: The data to sample
        :param n_samples: The number of samples to generate
        :param burn_in:\
            no of MCMC transitions to reach apparent equilibrium before\
            generating inference samples
        :param thinning:\
            sampling rate i.e. 5 = 1 sample for 5 generated steps
        :param alpha_jump_scale:\
            scaling of t transition distribution for MH jumps in alpha\
            direction
        :param alpha_min: The minimum value of alpha
        :param alpha_max: The maximum value of alpha
        :param beta_jump_scale:\
            scaling of t transition distribution for MH jumps in beta direction
        :param beta_min: The minimum value of beta
        :param beta_max: The maximum value of beta
        :param degrees_of_freedom:\
            The number of degrees of freedom to jump around with
        """

        MachineVertex.__init__(self, label="MCMC Node", constraints=None)
        self._data = data
        self._n_samples = n_samples
        self._burn_in = burn_in
        self._thinning = thinning
        self._alpha_jump_scale = alpha_jump_scale
        self._alpha_min = alpha_min
        self._alpha_max = alpha_max
        self._beta_jump_scale = beta_jump_scale
        self._beta_min = beta_min
        self._beta_max = beta_max
        self._degrees_of_freedom = degrees_of_freedom
        self._seed = seed
        self._send_timer = send_timer
        self._receive_timer = receive_timer
        self._window_size = window_size
        self._n_sequences = n_sequences
        self._data_partition_name = data_partition_name
        self._acknowledge_partition_name = acknowledge_partition_name
        self._data_tag = data_tag

        self._data_size = (
            (len(self._data) * self._DATA_ELEMENT_TYPE.size) +
            self._DATA_COUNT_TYPE.size
        )
        self._sdram_usage = (
            self._N_PARAMETER_BYTES + self._data_size
        )

        self._mcmc_vertices = list()
        self._mcmc_placements = list()
        self._data_receiver = dict()

    def register_processor(self, mcmc_vertex):
        self._mcmc_vertices.append(mcmc_vertex)

    @property
    def n_samples(self):
        return self._n_samples

    @property
    def burn_in(self):
        return self._burn_in

    @property
    def thinning(self):
        return self._thinning

    @property
    def alpha_jump_scale(self):
        return self._alpha_jump_scale

    @property
    def beta_jump_scale(self):
        return self._beta_jump_scale

    @property
    def alpha_min(self):
        return self._alpha_min

    @property
    def alpha_max(self):
        return self._alpha_max

    @property
    def beta_min(self):
        return self._beta_min

    @property
    def beta_max(self):
        return self._beta_max

    @property
    def degrees_of_freedom(self):
        return self._degrees_of_freedom

    @property
    def n_data_points(self):
        return len(self._data)

    def _is_receiver_placement(self, placement):
        x = placement.x
        y = placement.y
        if (x, y) not in self._data_receiver:
            self._data_receiver[x, y] = placement.p
            return True
        return self._data_receiver[(x, y)] == placement.p

    def get_data_window_size(self, placement):
        if self._is_receiver_placement(placement):
            return self._window_size
        return 0

    def get_sequence_mask(self, placement, routing_info):
        if self._is_receiver_placement(placement):
            mask = routing_info.get_routing_info_from_pre_vertex(
                self, self._data_partition_name).first_mask
            return ~mask & 0xFFFFFFFF
        return 0

    def get_acknowledge_key(self, placement, routing_info):
        if self._is_receiver_placement(placement):
            key = routing_info.get_first_key_from_pre_vertex(
                placement.vertex, self._acknowledge_partition_name)
            return key
        return 0

    @property
    def data_tag(self):
        return self._data_tag

    @property
    def acknowledge_timer(self):
        return self._receive_timer

    @property
    def seed(self):
        if self._seed is None:
            return [random.randint(0, 0xFFFFFFFF) for _ in range(5)]
        return self._seed

    @property
    def data_partition_name(self):
        return self._data_partition_name

    @property
    def acknowledge_partition_name(self):
        return self._acknowledge_partition_name

    @property
    @overrides(MachineVertex.resources_required)
    def resources_required(self):
        sdram = self._N_PARAMETER_BYTES + self._data_size
        sdram += len(self._mcmc_vertices) * self._KEY_ELEMENT_TYPE.size

        resources = ResourceContainer(
            dtcm=DTCMResource(0),
            sdram=SDRAMResource(sdram),
            cpu_cycles=CPUCyclesPerTickResource(0),
            iptags=[], reverse_iptags=[])
        return resources

    @overrides(AbstractHasAssociatedBinary.get_binary_file_name)
    def get_binary_file_name(self):
        return "mcmc_coordinator.aplx"

    @overrides(AbstractHasAssociatedBinary.get_binary_start_type)
    def get_binary_start_type(self):
        return ExecutableStartType.SYNC

    @inject_items({
        "placements": "MemoryPlacements",
        "routing_info": "MemoryRoutingInfos"})
    @overrides(
        AbstractGeneratesDataSpecification.generate_data_specification,
        additional_arguments=["placements", "routing_info"])
    def generate_data_specification(
            self, spec, placement, placements, routing_info):

        # Reserve and write the parameters region
        region_size = self._N_PARAMETER_BYTES + self._data_size
        region_size += len(self._mcmc_vertices) * self._KEY_ELEMENT_TYPE.size
        spec.reserve_memory_region(0, region_size)
        spec.switch_write_focus(0)

        # Get the placement of the vertices and find out how many chips
        # are needed
        keys = list()
        for vertex in self._mcmc_vertices:
            mcmc_placement = placements.get_placement_of_vertex(vertex)
            self._mcmc_placements.append(mcmc_placement)
            if self._is_receiver_placement(mcmc_placement):
                key = routing_info.get_first_key_from_pre_vertex(
                    vertex, self._acknowledge_partition_name)
                keys.append(key)
        keys.sort()

        # Write the data size in words
        spec.write_value(
            len(self._data) * (float(self._DATA_ELEMENT_TYPE.size) / 4.0),
            data_type=self._DATA_COUNT_TYPE)

        # Write the number of chips
        spec.write_value(len(keys), data_type=DataType.UINT32)

        # Write the key
        routing_info = routing_info.get_routing_info_from_pre_vertex(
            self, self._data_partition_name)
        spec.write_value(routing_info.first_key, data_type=DataType.UINT32)

        # Write the window size
        spec.write_value(self._window_size, data_type=DataType.UINT32)

        # Write the sequence mask
        spec.write_value(
            ~routing_info.first_mask & 0xFFFFFFFF, data_type=DataType.UINT32)

        # Write the timer
        spec.write_value(self._send_timer, data_type=DataType.UINT32)

        # Write the data - Arrays must be 32-bit values, so convert
        data = numpy.array(self._data, dtype=self._NUMPY_DATA_ELEMENT_TYPE)
        spec.write_array(data.view(numpy.uint32))

        # Write the keys
        spec.write_array(keys)

        # End the specification
        spec.end_specification()

    @overrides(AbstractProvidesNKeysForPartition.get_n_keys_for_partition)
    def get_n_keys_for_partition(self, partition, graph_mapper):
        return self._n_sequences

    def read_samples(self, buffer_manager):
        """ Read back the samples
        """
        progress = ProgressBar(len(self._mcmc_placements), "Reading results")
        samples = list()
        for placement in self._mcmc_placements:

            # Read the data recorded
            samples.append(
                placement.vertex.read_samples(buffer_manager, placement))
            progress.update()
        progress.end()

        # Merge all the arrays
        return numpy.hstack(samples)
