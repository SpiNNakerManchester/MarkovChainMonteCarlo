import spinnaker_graph_front_end as g

from mcmc_vertex import MCMCVertex
from mcmc_coordinator_vertex import MCMCCoordinatorVertex

from pacman.model.constraints.placer_constraints\
    .placer_chip_and_core_constraint import PlacerChipAndCoreConstraint
from pacman.model.graphs.machine.impl.machine_edge import MachineEdge

from spinnman.model.cpu_state import CPUState

import numpy
import logging
import time
from spinn_front_end_common.utilities import helpful_functions

logger = logging.getLogger(__name__)

data_points = [
    2.37706, 4.51142, 0.406605, 0.909418, 0.642899,
    1.21925, 1.47647, -2.95771, -0.801802, -1.86529,
    -75.6865, -1.94398, 1.57055, 2.53382, 0.783884,
    1.16725, 1.16995, 0.367477, -1.24639, -2.29897,
    0.461939, -0.126669, 0.0965992, 1.56107, 0.747027,
    66.7057, 1.07821, -0.125864, -0.693059, 7.48744,
    1.94184, -0.439164, 3.64695, -63.2296, 0.783037,
    2.26351, 1.30222, 0.542981, 3.78199, -37.1692,
    1.54959, 0.485336, 1.02509, -0.204211, 0.164426,
    -13.1977, 0.650243, 0.671339, 2.93511, 0.788114
]

# The data to use
# data_points = [
#     2.37706, 4.51142, 0.406605, 0.909418, 0.642899,
#     1.21925, 1.47647, -2.95771, -0.801802, -1.86529
# ]

seed = [
    123456789, 234567891, 345678912, 456789123, 0
]

# no of MCMC transitions to reach apparent equilibrium before generating
# inference samples
burn_in = 2000

# subsequent MCMC samples are correlated, so thin the chain to avoid this
thinning = 5

# number of posterior samples required
n_samples = 100

# scaling of t transition distribution for MH jumps in alpha direction
alpha_jump_scale = 0.8

# scaling of t transition distribution for MH jumps in beta direction
beta_jump_scale = 0.25

# The number of degrees of freedom to jump around with
degrees_of_freedom = 3.0

# specification of prior knowledge about lighthouse position
#
# you know that the lighthouse is no further than 3.0 units along shore
# from the reference zero position
#
# lighthouse cannot be closer than 0.2 from the shore because of rocks, or
# further than 2.0 because of shipping lanes
alpha_min = -3.0
alpha_max = 3.0
beta_min = 0.2
beta_max = 2.0

# If the machine is a generated one, set the number of chips to a single
# 48-node board
n_chips_required = None
if g.is_allocated_machine:
    n_chips_required = 1104

# Set up the simulation
g.setup(n_chips_required=n_chips_required)

# Get the number of cores available for use
n_cores = 0
machine = g.machine()

# Create a coordinator for each board
coordinators = dict()
boards = dict()
for chip in machine.ethernet_connected_chips:

    # Create a coordinator
    coordinator = MCMCCoordinatorVertex(
        data_points, n_samples, burn_in, thinning,
        alpha_jump_scale, alpha_min, alpha_max,
        beta_jump_scale, beta_min, beta_max,
        degrees_of_freedom, seed)
    g.add_machine_vertex_instance(coordinator)

    # Put the coordinator on the Ethernet chip
    coordinator.add_constraint(PlacerChipAndCoreConstraint(chip.x, chip.y))
    coordinators[chip.x, chip.y] = coordinator
    boards[chip.x, chip.y] = chip.ip_address

# Go through all the chips and add the workhorses
n_workers = 0
for chip in machine.chips:

    # Count the cores in the processor (-1 if this chip also has a coordinator)
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
    for core in range(n_cores):

        # Create the vertex and add it to the graph
        vertex = MCMCVertex(coordinator)
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
        error_msg = "Some cores have failed ({} rte, {} wdog)".format(
            error, watchdog)
        print error_msg
        input("PRESS ENTER TO CONTINUE")
        raise Exception(error_msg)
    running = txrx.get_core_state_count(30, CPUState.RUNNING)

# Get the data back
samples = list()
for coordinator in coordinators.itervalues():
    samples.append(coordinator.read_samples(g.buffer_manager()))
samples = numpy.hstack(samples)

# Close the machine
g.stop()

# Save the results
numpy.save("results.npy", samples)
numpy.savetxt("results.csv", samples, fmt="%f", delimiter=",")
