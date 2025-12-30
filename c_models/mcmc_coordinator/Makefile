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

CUR_DIR := $(abspath $(dir $(lastword $(MAKEFILE_LIST))))
FEC_INSTALL_DIR := $(strip $(if $(FEC_INSTALL_DIR), $(FEC_INSTALL_DIR), $(abspath $(CUR_DIR)/../../../SpiNNFrontEndCommon/c_common/front_end_common_lib)))

# The name of the application to be built (binary will be this with a `.aplx`
# extension)
APP = mcmc_coordinator

# Directory to create APLX files in (must include trailing slash)
APP_OUTPUT_DIR = ../../mcmc/model_binaries/
# key for the database in this APP_OUTPUT_DIR
DATABASE_KEY = M

CFLAGS += ${OTIME}

# SOURCE_DIRS = .
SOURCES = mcmc_coordinator.c

# The spinnaker_tools standard makefile
include $(FEC_INSTALL_DIR)/make/fec.mk

clean:
	$(RM) -r $(CURRENT_DIR)/build
