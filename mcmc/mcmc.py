import spinnaker_graph_front_end as g

from mcmc.mcmc_vertex import MCMCVertex
from mcmc.mcmc_coordinator_vertex import MCMCCoordinatorVertex
from mcmc import model_binaries

from pacman.model.constraints.placer_constraints\
    .placer_chip_and_core_constraint import PlacerChipAndCoreConstraint
from pacman.model.graphs.machine import MachineEdge

from spinnman.model.enums.cpu_state import CPUState

import numpy
import logging
import time

logger = logging.getLogger(__name__)


def run_mcmc(
        self, model, data, n_samples, burn_in=2000, thinning=5,
        degrees_of_freedom=3.0, seed=None, n_chips=None):
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

    :return: The samples read
    :rtype: A numpy array with fields for each model state variable
    """

    # If the machine is a generated one, set the number of chips to a single
    # 48-node board
    n_chips_required = None
    if g.is_allocated_machine:
        if n_chips is None:
            raise Exception(
                "n_chips must be specified to use an allocated machine")
        n_chips_required = n_chips

    # Set up the simulation
    g.setup(
        n_chips_required=n_chips_required, model_binary_module=model_binaries)

    # Get the number of cores available for use
    n_cores = 0
    machine = g.machine()

    # Create a coordinator for each board
    coordinators = dict()
    boards = dict()
    for chip in machine.ethernet_connected_chips:

        # Create a coordinator
        coordinator = MCMCCoordinatorVertex(
            data, n_samples, burn_in, thinning,
            degrees_of_freedom, seed)
        g.add_machine_vertex_instance(coordinator)

        # Put the coordinator on the Ethernet chip
        coordinator.add_constraint(PlacerChipAndCoreConstraint(chip.x, chip.y))
        coordinators[chip.x, chip.y] = coordinator
        boards[chip.x, chip.y] = chip.ip_address

    # Go through all the chips and add the workhorses
    n_workers = 0
    for chip in machine.chips:

        # Count the cores in the processor
        # (-1 if this chip also has a coordinator)
        n_cores = len([p for p in chip.processors if not p.is_monitor])
        if (chip.x, chip.y) in coordinators:
            n_cores -= 1

        # Find the coordinator for the board (or 0, 0 if it is missing)
        eth_x = chip.nearest_ethernet_x
        eth_y = chip.nearest_ethernet_y
        coordinator = coordinators.get((eth_x, eth_y))
        if coordinator is None:
            print "Warning - couldn't find {}, {}".format(eth_x, eth_y)
            coordinator = coordinators[0, 0]

        # Add a vertex for each core
        for _ in range(n_cores):

            # Create the vertex and add it to the graph
            vertex = MCMCVertex(coordinator, model)
            n_workers += 1
            g.add_machine_vertex_instance(vertex)

            # Put the vertex on the same board as the coordinator
            vertex.add_constraint(PlacerChipAndCoreConstraint(chip.x, chip.y))

            # Add an edge from the coordinator to the vertex, to send the data
            g.add_machine_edge_instance(
                MachineEdge(coordinator, vertex),
                coordinator.data_partition_name)

            # Add an edge from the vertex to the coordinator,
            # to send acknowledgement
            g.add_machine_edge_instance(
                MachineEdge(vertex, coordinator),
                coordinator.acknowledge_partition_name)

    # Run the simulation
    g.run(None)

    # Wait for the application to finish
    txrx = g.transceiver()
    logger.info("Running {} worker cores".format(n_workers))
    logger.info("Waiting for application to finish...")
    running = txrx.get_core_state_count(30, CPUState.RUNNING)
    while running > 0:
        time.sleep(0.5)
        error = txrx.get_core_state_count(30, CPUState.RUN_TIME_EXCEPTION)
        watchdog = txrx.get_core_state_count(30, CPUState.WATCHDOG)
        if error > 0 or watchdog > 0:
            error_msg = "Some cores have failed ({} RTE, {} WDOG)".format(
                error, watchdog)
            raise Exception(error_msg)
        running = txrx.get_core_state_count(30, CPUState.RUNNING)

    # Get the data back
    samples = list()
    for coordinator in coordinators.itervalues():
        samples.append(coordinator.read_samples(g.buffer_manager()))
    samples = numpy.hstack(samples)

    # Close the machine
    g.stop()

    return samples
