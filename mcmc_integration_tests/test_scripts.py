# Copyright (c) 2019-2021 The University of Manchester
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

from spinnaker_testbase import ScriptChecker
from unittest import SkipTest  # pylint: disable=unused-import


class TestScripts(ScriptChecker):
    """
    This file tests the scripts as configured in script_builder.py

    Please do not manually edit this file.
    It is rebuilt each time SpiNNakerManchester/IntegrationTests is run

    If it is out of date please edit and run script_builder.py
    Then the new file can be added to github for reference only.
    """
# flake8: noqa

    def test_mcmc_examples_lighthouse_lighthouse_float_model(self):
        raise SkipTest("Not a run script")
        self.check_script("mcmc_examples/lighthouse/lighthouse_float_model.py")

    def test_mcmc_examples_lighthouse_lighthouse(self):
        self.check_script("mcmc_examples/lighthouse/lighthouse.py")

    def test_mcmc_examples_lighthouse_lighthouse_fixed_point_model(self):
        raise SkipTest("Not a run script")
        self.check_script("mcmc_examples/lighthouse/lighthouse_fixed_point_model.py")

    def test_mcmc_examples_lighthouse_lighthouse_model(self):
        raise SkipTest("Not a run script")
        self.check_script("mcmc_examples/lighthouse/lighthouse_model.py")

    def test_mcmc_examples_arma_arma(self):
        self.check_script("mcmc_examples/arma/arma.py")

    def test_mcmc_examples_arma_arma_float_model(self):
        raise SkipTest("Not a run script")
        self.check_script("mcmc_examples/arma/arma_float_model.py")

    def test_mcmc_examples_arma_arma_fixed_point_model(self):
        raise SkipTest("Not a run script")
        self.check_script("mcmc_examples/arma/arma_fixed_point_model.py")

    def test_mcmc_examples_arma_arma_model(self):
        raise SkipTest("Not a run script")
        self.check_script("mcmc_examples/arma/arma_model.py")
