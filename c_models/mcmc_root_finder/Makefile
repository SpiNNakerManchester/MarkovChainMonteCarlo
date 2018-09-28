# An example Makefile suitable for use in SpiNNaker applications using the
# spinnaker_tools libraries and makefiles.

# The name of the application to be built (binary will be this with a `.aplx`
# extension)
APP = mcmc_root_finder

# Directory to create APLX files in (must include trailing slash)
APP_OUTPUT_DIR = ../../mcmc/model_binaries/

# Directory to place compilation artefacts (must include trailing slash)
BUILD_DIR = ./build/

CFLAGS += ${OTIME} -fcx-limited-range  # check list to add here
LFLAGS += ${OTIME} -fcx-limited-range

SOURCES = mcmc_root_finder.c

# SOURCE_DIRS is more than just src/ in this instance
SOURCE_DIRS := src/ ../mcmc_models/src/ ../mcmc_models/examples/arma/src/

# The spinnaker_tools standard makefile
include $(SPINN_DIRS)/make/local.mk

all: $(APP_OUTPUT_DIR)$(APP).aplx
