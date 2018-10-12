# The name of the application to be built (binary will be this with a `.aplx`
# extension)
APP = lighthouse

# Directory to place compilation artefacts (must include trailing slash)
BUILD_DIR = ./build/

# The sources to be compiled
MAKEFILE_PATH := $(abspath $(lastword $(MAKEFILE_LIST)))
CURRENT_DIR := $(dir $(MAKEFILE_PATH))
SOURCE_DIR := $(abspath $(CURRENT_DIR))/src
SOURCE_DIRS := $(SOURCE_DIR)

SOURCES = lighthouse.c

CFLAGS += $(OTIME) 

include ../../Makefile.common
