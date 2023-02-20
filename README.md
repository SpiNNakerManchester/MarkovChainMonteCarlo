[![Python Build Status](https://github.com/SpiNNakerManchester/MarkovChainMonteCarlo/workflows/Python%20Actions/badge.svg?branch=master)](https://github.com/SpiNNakerManchester/MarkovChainMonteCarlo/actions?query=workflow%3A%22Python+Actions%22+branch%3Amaster)
[![C Build Status](https://github.com/SpiNNakerManchester/MarkovChainMonteCarlo/workflows/C%20Actions/badge.svg?branch=master)](https://github.com/SpiNNakerManchester/MarkovChainMonteCarlo/actionpiNNFrontEndCommon?branch=master)

# Markov Chain Monte Carlo Simulations On SpiNNaker

This project contains code for running MCMC Simulations on SpiNNaker.  This has been made general, to allow users to add their own MCMC models by providing a few simple components.

## Adding your own model

### C code
1. Create a C header file containing:
    - A structure called ```struct mcmc_params``` which contains every fixed parameter that you require.
    - A structure called ```struct mcmc_state``` which contains each state variable that changes during the execution that you require.
    
    An example can be seen in ```c_models/mcmc_models/examples/lighthouse/lighthouse.h```

1. Create a C source file containing and implementation of each of the functions in ```c_models/mcmc_models/mcmc_model.h```

    An example can be seen in ```c_models/mcmc_models/examples/lighthouse/lighthouse.c```

1. Create a Makefile for your C code, which includes your source file.  Call ```make``` to build your code.

    An example can be seen in ```c_models/mcmc_models/examples/lighthouse/Makefile```

1. Create a Python file containing a class which extends ```mcmc.mcmc_model.MCMCModel``` and implement each of the abstract methods therein.

    An example can be seen in ```mcmc_examples/lighthouse/lighthouse_model.py```

1. Create a Python script for your example, which calls ```mcmc.mcmc_framework.run_mcmc```.  This returns the results of the simulation.

    An example can be seen in ```mcmc_examples/lighthouse/lighthouse.py```
