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
from spinn_front_end_common.utilities import helpful_functions
from spinn_front_end_common.interface.buffer_management \
    import recording_utilities
from spinn_front_end_common.utilities.utility_objs.executable_type \
    import ExecutableType

from spinn_machine.utilities.progress_bar import ProgressBar

from enum import Enum

class MCMCCholeskyRegions(Enum):
    """ Regions in the MCMCCholesky Data
    """
    PARAMETERS = 0
    # Note: recording region may also be necessary in order to read history


class MCMCCholeskyVertex(
        MachineVertex, AbstractHasAssociatedBinary,
        AbstractGeneratesDataSpecification):  #, AbstractReceiveBuffersToHost):
#        AbstractProvidesNKeysForPartition):
    """ A vertex that runs the (MCMC) Cholesky algorithm
    """

    def __init__(
            self, vertex, model):
        """
        :param vertex: The MCMC vertex associated with this vertex
        :param model: The model being simulated
        """

        MachineVertex.__init__(self, label="MCMC Cholesky Node",
                               constraints=None)
        self._vertex = vertex
        self._model = model

        vertex.coordinator.register_processor(self)

        # Other parameters may need to be added here
        self._n_other_params = 1

        self._n_parameter_bytes = self._n_other_params * 4

        self._sdram_usage = (self._n_parameter_bytes)

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
        return "mcmc_root_finder.aplx"

    @overrides(AbstractHasAssociatedBinary.get_binary_start_type)
    def get_binary_start_type(self):
        return ExecutableType.SYNC

    @inject_items({
        "routing_info": "MemoryRoutingInfos"})
    @overrides(
        AbstractGeneratesDataSpecification.generate_data_specification,
        additional_arguments=["routing_info"])
    def generate_data_specification(
            self, spec, placement, routing_info):

        # Reserve and write the parameters region
        spec.reserve_memory_region(
            MCMCCholeskyRegions.PARAMETERS.value, self._n_parameter_bytes)
        spec.switch_write_focus(MCMCCholeskyRegions.PARAMETERS.value)

        # We may need the recording region in this instance as well

        # Write the acknowledge key
        spec.write_value(self._vertex.get_result_key(
            placement, routing_info))

        # End the specification
        spec.end_specification()

    def read_samples(self, buffer_manager, placement):
        """ Read back the samples (dummy call)
        """
#        print 'There are no samples to read back on a Cholesky vertex'
        return None
