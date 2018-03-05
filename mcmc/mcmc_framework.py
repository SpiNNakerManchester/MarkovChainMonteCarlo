import spinnaker_graph_front_end as g

from .mcmc_vertex import MCMCVertex
from .mcmc_coordinator_vertex import MCMCCoordinatorVertex
from .mcmc_root_finder_vertex import MCMCRootFinderVertex
from . import model_binaries

from pacman.model.constraints.placer_constraints\
    .chip_and_core_constraint import ChipAndCoreConstraint
from pacman.model.constraints.placer_constraints\
    .same_chip_as_constraint import SameChipAsConstraint
from pacman.model.graphs.machine import MachineEdge

from spinnman.model.enums.cpu_state import CPUState

from spinn_front_end_common.utilities import globals_variables

import numpy
import logging
import time

# timing
start_time = time.time()
logger = logging.getLogger(__name__)


def run_mcmc(
        model, data, n_samples, burn_in=2000, thinning=5,
        degrees_of_freedom=3.0, seed=None, n_chips=None, root_finder=False):
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

    :return: The samples read
    :rtype: A numpy array with fields for each model state variable
    """

    # Set up the simulation
    g.setup(n_chips_required=n_chips, model_binary_module=model_binaries)

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
        coordinator.add_constraint(ChipAndCoreConstraint(chip.x, chip.y))
        coordinators[chip.x, chip.y] = coordinator
        boards[chip.x, chip.y] = chip.ip_address

    # Go through all the chips and add the workhorses
    n_chips_on_machine = machine.n_chips
    n_workers = 0
    if (root_finder):
        n_root_finders = 0
    for chip in machine.chips:

        # Count the cores in the processor
        # (-1 if this chip also has a coordinator)
        n_cores = len([p for p in chip.processors if not p.is_monitor])
        if (chip.x, chip.y) in coordinators:
            n_cores -= 3  # coordinator and extra_monitor_support (2)
            if (root_finder):
                n_cores = n_cores / 2
        else:
            n_cores -= 1  # just extra_monitor_support
            if (root_finder):
                n_cores = n_cores / 2

        # Find the coordinator for the board (or 0, 0 if it is missing)
        eth_x = chip.nearest_ethernet_x
        eth_y = chip.nearest_ethernet_y
        coordinator = coordinators.get((eth_x, eth_y))
        if coordinator is None:
            print "Warning - couldn't find {}, {}".format(eth_x, eth_y)
            coordinator = coordinators[0, 0]

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
            vertex.add_constraint(ChipAndCoreConstraint(chip.x, chip.y))

            # Add an edge from the coordinator to the vertex, to send the data
            g.add_machine_edge_instance(
                MachineEdge(coordinator, vertex),
                coordinator.data_partition_name)

            # Add an edge from the vertex to the coordinator,
            # to send acknowledgement
            g.add_machine_edge_instance(
                MachineEdge(vertex, coordinator),
                coordinator.acknowledge_partition_name)

            if (root_finder):
                # Create a root finder vertex
                rf_vertex = MCMCRootFinderVertex(vertex, model)
                n_root_finders += 1
                g.add_machine_vertex_instance(rf_vertex)

                # put it on the same chip as the standard mcmc vertex?
                # no - put it on a "nearby" chip, however that works
                rf_vertex.add_constraint(ChipAndCoreConstraint(chip.x, chip.y))

                # Add an edge from mcmc vertex to root finder vertex,
                # to "send" the data - need to work this out
                g.add_machine_edge_instance(
                    MachineEdge(vertex, rf_vertex),
                    vertex.parameter_partition_name)  # key set here!

                # Add edge from root finder vertex back to mcmc vertex
                # to send acknowledgement / result - need to work this out
                g.add_machine_edge_instance(
                    MachineEdge(rf_vertex, vertex),
                    vertex.result_partition_name)  # "acknowledge key" set here

    start_computing_time = time.time()

    # Run the simulation
    #g.run(None)
    g.run_until_complete()

    mid_computing_time = time.time()

    # Wait for the application to finish
    txrx = g.transceiver()
    app_id = globals_variables.get_simulator()._app_id
    logger.info("Running {} worker cores".format(n_workers))
    if (root_finder):
        logger.info("Running {} root finder cores".format(n_root_finders))
    logger.info("Waiting for application to finish...")
    running = txrx.get_core_state_count(app_id, CPUState.RUNNING)
    # there are now cores doing extra_monitor etc.
    non_worker_cores = n_chips_on_machine + (2 * len(boards))  # extra_monitor + coordinator
    while running > non_worker_cores:
        time.sleep(0.5)
        error = txrx.get_core_state_count(app_id, CPUState.RUN_TIME_EXCEPTION)
        watchdog = txrx.get_core_state_count(app_id, CPUState.WATCHDOG)
        if error > 0 or watchdog > 0:
            error_msg = "Some cores have failed ({} RTE, {} WDOG)".format(
                error, watchdog)
            raise Exception(error_msg)
        running = txrx.get_core_state_count(app_id, CPUState.RUNNING)
        print 'running: ', running

    finish_computing_time = time.time()

    # Get the data back
    samples = list()
    for coordinator in coordinators.itervalues():
        samples.append(coordinator.read_samples(g.buffer_manager()))

    samples = numpy.hstack(samples)

    # Close the machine
    g.stop()

    finish_time = time.time()

    print("Overhead time is %s seconds" % (start_computing_time - start_time))
    print("Computing time is %s seconds"
          % (finish_computing_time - start_computing_time))
    print("run_until_complete takes %s seconds"
          % (mid_computing_time - start_computing_time))
    print("Data collecting time is %s seconds"
          % (finish_time - finish_computing_time))
    print("Overall running time is %s seconds" % (finish_time - start_time))

    return samples
