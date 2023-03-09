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
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

import spinnaker_graph_front_end as g

from spinnman.exceptions import SpinnmanException
from .mcmc_vertex import MCMCVertex
from .mcmc_coordinator_vertex import MCMCCoordinatorVertex
from .mcmc_root_finder_vertex import MCMCRootFinderVertex
from .mcmc_cholesky_vertex import MCMCCholeskyVertex
from . import model_binaries

from pacman.model.graphs.machine import MachineEdge

from spinnman.model.enums.cpu_state import CPUState

from spinn_front_end_common.data import FecDataView

import logging
import time

# timing
start_time = time.time()
logger = logging.getLogger(__name__)


def run_mcmc(
        model, data, n_samples, burn_in=2000, thinning=5,
        degrees_of_freedom=3.0, seed=None, n_chips=None, n_boards=None):
    """ Executes an MCMC model, returning the received samples

    :param model: The MCMCModel to be used
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
    :param n_chips: The number of chips to run the model on
    :param root_finder: Use the root finder by adding root finder vertices
    :param cholesky: Use the Cholesky algorithm by adding Cholesky vertices

    :return: The samples read as
        a numpy array with fields for each model state variable
    :rtype: dict(tuple(int,int),~numpy.ndarray)
    """

    # Set up the simulation
    g.setup(n_boards_required=n_boards, n_chips_required=n_chips,
            model_binary_module=model_binaries)

    # Get the number of cores available for use
    n_cores = 0
    machine = g.machine()

    # Create a coordinator for each board
    coordinators = dict()
    boards = dict()
    for chip in machine.ethernet_connected_chips:

        # Create a coordinator
        coordinator = MCMCCoordinatorVertex(
            model, data, n_samples, burn_in, thinning,
            degrees_of_freedom, seed)
        g.add_machine_vertex_instance(coordinator)

        # Put the coordinator on the Ethernet chip
        coordinator.set_fixed_location(chip.x, chip.y)
        coordinators[chip.x, chip.y] = coordinator
        boards[chip.x, chip.y] = chip.ip_address

    # Go through all the chips and add the workhorses
    n_chips_on_machine = machine.n_chips
    n_workers = 0
    if (model.root_finder):
        n_root_finders = 0
    if (model.cholesky):
        n_cholesky = 0
    for chip in machine.chips:

        # Count the cores in the processor
        # (-1 if this chip also has a coordinator)
        n_cores = len([p for p in chip.processors if not p.is_monitor])
        if (chip.x, chip.y) in coordinators:
            n_cores -= 3  # coordinator and extra_monitor_support (2)
            if (model.root_finder):
                if (model.cholesky):
                    n_cores = n_cores // 3
                else:
                    n_cores = n_cores // 2
        else:
            n_cores -= 1  # just extra_monitor_support
            if (model.root_finder):
                if (model.cholesky):
                    n_cores = n_cores // 3
                else:
                    n_cores = n_cores // 2

        # Find the coordinator for the board (or 0, 0 if it is missing)
        eth_x = chip.nearest_ethernet_x
        eth_y = chip.nearest_ethernet_y
        coordinator = coordinators.get((eth_x, eth_y))
        if coordinator is None:
            print("Warning - couldn't find {}, {} for chip {}, {}".format(
                eth_x, eth_y, chip.x, chip.y))
            coordinator = coordinators[0, 0]
            print("Using coordinator ", coordinator)

        # hard-code remove some cores (chip power monitor etc.) just
        # to see what happens
#        n_cores -= non_worker_cores_per_chip
#        print 'n_cores: ', n_cores

        # Add a vertex for each core
        for _ in range(n_cores):

            # Create the vertex and add it to the graph
            vertex = MCMCVertex(coordinator, model)
            n_workers += 1
            g.add_machine_vertex_instance(vertex)

            # Put the vertex on the same board as the coordinator
            vertex.set_fixed_location(chip.x, chip.y)

            # Add an edge from the coordinator to the vertex, to send the data
            g.add_machine_edge_instance(
                MachineEdge(coordinator, vertex),
                coordinator.data_partition_name)

            # Add an edge from the vertex to the coordinator,
            # to send acknowledgement
            g.add_machine_edge_instance(
                MachineEdge(vertex, coordinator),
                coordinator.acknowledge_partition_name)

            if (model.root_finder):
                # Create a root finder vertex
                rf_vertex = MCMCRootFinderVertex(vertex, model)
                n_root_finders += 1
                g.add_machine_vertex_instance(rf_vertex)

                # put it on the same chip as the standard mcmc vertex?
                # no - put it on a "nearby" chip, however that works
                rf_vertex.set_fixed_location(chip.x, chip.y)

                # Add an edge from mcmc vertex to root finder vertex,
                # to "send" the data - need to work this out
                g.add_machine_edge_instance(
                    MachineEdge(vertex, rf_vertex),
                    vertex.parameter_partition_name)

                # Add edge from root finder vertex back to mcmc vertex
                # to send acknowledgement / result - need to work this out
                g.add_machine_edge_instance(
                    MachineEdge(rf_vertex, vertex),
                    vertex.result_partition_name)

            if (model.cholesky):
                # Create a Cholesky vertex
                cholesky_vertex = MCMCCholeskyVertex(vertex, model)
                n_cholesky += 1
                g.add_machine_vertex_instance(cholesky_vertex)

                # put it on the same chip as the standard mcmc vertex?
                # no - put it on a "nearby" chip, however that works
                cholesky_vertex.set_fixed_location(chip.x, chip.y)

                # Add an edge from mcmc vertex to Cholesky vertex,
                # to "send" the data - need to work this out
                g.add_machine_edge_instance(
                    MachineEdge(vertex, cholesky_vertex),
                    vertex.cholesky_partition_name)

                # Add edge from Cholesky vertex back to mcmc vertex
                # to send acknowledgement / result - need to work this out
                g.add_machine_edge_instance(
                    MachineEdge(cholesky_vertex, vertex),
                    vertex.cholesky_result_partition_name)

    start_computing_time = time.time()

    logger.info("n_chips_on_machine %s", n_chips_on_machine)
    logger.info("Running %s worker cores", n_workers)
    if (model.root_finder):
        logger.info("Running %s root finder cores", n_root_finders)
    if (model.cholesky):
        logger.info("Running %s Cholesky cores", n_cholesky)

    # Run the simulation
    g.run_until_complete()

    mid_computing_time = time.time()

    # Wait for the application to finish
    txrx = FecDataView.get_transceiver()
    app_id = FecDataView.get_app_id()
    logger.info("Running %s worker cores", n_workers)
    if (model.root_finder):
        logger.info("Running %s root finder cores", n_root_finders)
    if (model.cholesky):
        logger.info("Running %s Cholesky cores", n_cholesky)
    logger.info("Waiting for application to finish...")
    running = txrx.get_core_state_count(app_id, CPUState.RUNNING)
    # there are now cores doing extra_monitor etc.
    non_worker_cores = n_chips_on_machine + (2 * len(boards))
    while running > non_worker_cores:
        time.sleep(0.5)
        error = txrx.get_core_state_count(app_id, CPUState.RUN_TIME_EXCEPTION)
        watchdog = txrx.get_core_state_count(app_id, CPUState.WATCHDOG)
        if error > 0 or watchdog > 0:
            error_msg = "Some cores have failed ({} RTE, {} WDOG)".format(
                error, watchdog)
            raise SpinnmanException(error_msg)
        running = txrx.get_core_state_count(app_id, CPUState.RUNNING)
        print('running: ', running)

    finish_computing_time = time.time()

    # Get the data back
    samples = dict()
    for coord, coordinator in coordinators.items():
        samples[coord[0], coord[1]] = coordinator.read_samples(
            g.buffer_manager())

    # Close the machine
    g.stop()

    finish_time = time.time()

    # Note: this timing appears to be incorrect now; needs looking at
    print("Overhead time is %s seconds" % (start_computing_time - start_time))
    print("Computing time is %s seconds"
          % (finish_computing_time - start_computing_time))
    print("run_until_complete takes %s seconds"
          % (mid_computing_time - start_computing_time))
    print("Data collecting time is %s seconds"
          % (finish_time - finish_computing_time))
    print("Overall running time is %s seconds" % (finish_time - start_time))

    return samples
