# Copyright (c) 2017 The University of Manchester
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

import os
import unittest
from spinn_utilities.config_holder import run_config_checks
from spinnaker_graph_front_end.config_setup import unittest_setup


class TestCfgChecker(unittest.TestCase):

    def setUp(self):
        unittest_setup()

    def test_cfg_check(self):
        unittests = os.path.dirname(__file__)
        parent = os.path.dirname(unittests)
        mcmc = os.path.join(parent, "mcmc")
        integration_tests = os.path.join(parent, "mcmc_integration_tests")
        mcmc_examples = os.path.join(parent, "mcmc_examples")
        repeaters = ["placer", "router", "info_allocator", "compressor"]
        run_config_checks(
            directories=[mcmc, integration_tests, unittests, mcmc_examples],
            repeaters=repeaters)
