from pacman.model.graphs.machine import MachineVertex
from pacman.model.resources import ResourceContainer, ConstantSDRAM
from pacman.model.decorators.overrides import overrides
from pacman.executor.injection_decorator import inject_items

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

        MachineVertex.__init__(self, label="MCMC RF Node", constraints=None)
        self._vertex = vertex
        self._model = model

        vertex.coordinator.register_processor(self)

        # Other parameters need to be added here
        self._n_other_params = 1

        self._n_parameter_bytes = self._n_other_params * 4

        self._sdram_usage = (self._n_parameter_bytes)

    @property
    @overrides(MachineVertex.resources_required)
    def resources_required(self):
        resources = ResourceContainer(
            sdram=ConstantSDRAM(self._sdram_usage))
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
            MCMCRootFinderRegions.PARAMETERS.value, self._n_parameter_bytes)
        spec.switch_write_focus(MCMCRootFinderRegions.PARAMETERS.value)

        # Write the acknowledge key
        spec.write_value(self._vertex.get_result_key(
            placement, routing_info))

        # End the specification
        spec.end_specification()

    def read_samples(self, buffer_manager, placement):
        """ Read back the samples (dummy call)
        """
#        print 'There are no samples to read back on a root finder vertex'
        return None
