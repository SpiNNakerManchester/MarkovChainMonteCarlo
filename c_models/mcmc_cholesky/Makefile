# An example Makefile suitable for use in SpiNNaker applications using the
# spinnaker_tools libraries and makefiles.

# The name of the application to be built (binary will be this with a `.aplx`
# extension)
APP = mcmc_cholesky

# Directory to create APLX files in (must include trailing slash)
APP_OUTPUT_DIR = ../../mcmc/model_binaries/

# Directory to place compilation artefacts (must include trailing slash)
BUILD_DIR = ./build/

CFLAGS += ${OTIME}  # check list to add here
LFLAGS += ${OTIME}

# SOURCE_DIRS = .
SOURCES = mcmc_cholesky.c

# The spinnaker_tools standard makefile
include $(SPINN_DIRS)/make/local.mk

all: $(APP_OUTPUT_DIR)$(APP).aplx
