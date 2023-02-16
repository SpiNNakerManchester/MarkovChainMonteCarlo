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
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

# The name of the application to be built (binary will be this with a `.aplx`
# extension)
APP = arma

# Directory to place compilation artefacts (must include trailing slash)
BUILD_DIR = ./build/

# The sources to be compiled
MAKEFILE_PATH := $(abspath $(lastword $(MAKEFILE_LIST)))
CURRENT_DIR := $(dir $(MAKEFILE_PATH))
SOURCE_DIR := $(abspath $(CURRENT_DIR))/src
SOURCE_DIRS := $(SOURCE_DIR)

SOURCES = arma.c

# will editing this here work to add -std=c99? yes but it doesn't help... 
CFLAGS += -fcx-limited-range $(OSPACE) -l m
# OFLAGS += -fcx-limited-range 

include ../../mcmc_common.mk
