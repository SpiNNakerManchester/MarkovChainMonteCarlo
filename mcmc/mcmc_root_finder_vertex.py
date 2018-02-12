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
# from spinn_front_end_common.abstract_models\
#     .abstract_provides_n_keys_for_partition \
#     import AbstractProvidesNKeysForPartition
from spinn_front_end_common.utilities import helpful_functions
from spinn_front_end_common.interface.buffer_management \
    import recording_utilities
from spinn_front_end_common.utilities.utility_objs.executable_type \
    import ExecutableType

from spinn_machine.utilities.progress_bar import ProgressBar

from enum import Enum
import numpy
# import random

class MCMCRootFinderRegions(Enum):
    """ Regions in the MCMCRootFinder Data
    """
    PARAMETERS = 0
#     RECORDING = 0  # not sure this is needed?
#     PARAMETERS = 1
#     MODEL_PARAMS = 2
#     MODEL_STATE = 3


class MCMCRootFinderVertex(
        MachineVertex, AbstractHasAssociatedBinary,
        AbstractGeneratesDataSpecification, AbstractReceiveBuffersToHost):
#        AbstractProvidesNKeysForPartition):
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

        # Do we need these functions in order to work out "recording sizes"?
#         state = self._get_model_state_array()
#         self._recording_size = len(state) * 4  # not really recording size?
        # params = self._get_model_parameters_array()

        # (We definitely need to use the "state" array, that's for sure...)
        # Other parameters need to be added here
        self._n_other_params = 1

        self._n_parameter_bytes = self._n_other_params * 4

        # here we use n_other_params, plus 1 d.o.f which is type dependent
#         if (model.get_parameters()[0].data_type is numpy.float64):
#             self._n_parameter_bytes = (self._n_other_params * 4) + (1 * 8)
#         else:
#             self._n_parameter_bytes = (self._n_other_params + 1) * 4

        self._sdram_usage = (
            self._n_parameter_bytes )  #  + self._recording_size)


#     def _get_model_parameters_array(self):
#         parameters = self._model.get_parameters()
#         numpy_format = list()
#         numpy_values = list()
#         for i, param in enumerate(parameters):
#             if (param.data_type is numpy.float64):
#                 numpy_format.append(('f{}'.format(i), param.data_type))
#                 numpy_values.append(param.value)
#             elif (param.data_type is numpy.float32):
#                 numpy_format.append(('f{}'.format(i), param.data_type))
#                 numpy_values.append(param.value)
#             elif (param.data_type is DataType.S1615):
#                 numpy_format.append(('f{}'.format(i), numpy.uint32))
#                 numpy_values.append(
#                     int(param.value * float(DataType.S1615.scale)))
#             else:
#                 # throw exception for unknown data type
#                 raise Exception(
#                     "Error: unsupported data type used for model parameters")
#         return numpy.array(
#             [tuple(numpy_values)], dtype=numpy_format).view("uint32")


#     def _get_model_state_array(self):
#         state = self._model.get_state_variables()
#         numpy_format = list()
#         numpy_values = list()
#         for i, param in enumerate(state):
#             if (param.data_type is numpy.float64):
#                 numpy_format.append(('f{}'.format(i), param.data_type))
#                 numpy_values.append(param.initial_value)
#             elif (param.data_type is numpy.float32):
#                 numpy_format.append(('f{}'.format(i), param.data_type))
#                 numpy_values.append(param.initial_value)
#             elif (param.data_type is numpy.uint32):
#                 numpy_format.append(('f{}'.format(i), param.data_type))
#                 numpy_values.append(param.initial_value)
#             elif (param.data_type is DataType.S1615):
#                 numpy_format.append(('f{}'.format(i), numpy.uint32))
#                 numpy_values.append(
#                     int(param.initial_value * float(DataType.S1615.scale)))
#             else:
#                 # throw exception for unknown data type
#                 raise Exception(
#                     "Error: unsupported data type for model state params")
#         return numpy.array(
#             [tuple(numpy_values)], dtype=numpy_format).view("uint32")

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
            self, spec, placement, routing_info):  # , tags):

        # Reserve and write the parameters region
        spec.reserve_memory_region(
            MCMCRootFinderRegions.PARAMETERS.value, self._n_parameter_bytes)
        spec.switch_write_focus(MCMCRootFinderRegions.PARAMETERS.value)

        # Write the acknowledge key
        spec.write_value(self._vertex.get_result_key(
            placement, routing_info))

#         # Reserve and write the model state
#         state = self._get_model_state_array()
#         spec.reserve_memory_region(
#             MCMCRootFinderRegions.MODEL_STATE.value, len(state) * 4)
#         spec.switch_write_focus(MCMCRootFinderRegions.MODEL_STATE.value)
#         spec.write_array(state)

        # End the specification
        spec.end_specification()

    def read_samples(self, buffer_manager, placement):
        """ Read back the samples (dummy call)
        """
        print('There are no samples to read back on a root finder vertex')

    def get_minimum_buffer_sdram_usage(self):
        return 1024

    def get_n_timesteps_in_buffer_space(self, buffer_space, machine_time_step):
        return recording_utilities.get_n_timesteps_in_buffer_space(
            buffer_space, 4)

    def get_recorded_region_ids(self):
        return [0]

    def get_recording_region_base_address(self, txrx, placement):
        return helpful_functions.locate_memory_region_for_placement(
            placement, MCMCRootFinderRegions.PARAMETERS.value, txrx)
