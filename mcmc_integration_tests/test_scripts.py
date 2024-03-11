# Copyright (c) 2019 The University of Manchester
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     https://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from spinnaker_testbase import ScriptChecker


class TestScripts(ScriptChecker):
    """
    This file tests the scripts as configured in script_builder.py

    Please do not manually edit this file.
    It is rebuilt each time SpiNNakerManchester/IntegrationTests is run

    If it is out of date please edit and run script_builder.py
    Then the new file can be added to github for reference only.
    """
# flake8: noqa

    def test_mcmc_examples_arma_arma(self):
        self.check_script("mcmc_examples/arma/arma.py")

    # Not testing file due to: Not a run script
    # mcmc_examples/arma/arma_float_model.py

    # Not testing file due to: Not a run script
    # mcmc_examples/arma/arma_model.py

    # Not testing file due to: Not a run script
    # mcmc_examples/arma/arma_fixed_point_model.py

    # Not testing file due to: Not a run script
    # mcmc_examples/lighthouse/lighthouse_model.py

    # Not testing file due to: Not a run script
    # mcmc_examples/lighthouse/lighthouse_float_model.py

    # Not testing file due to: Not a run script
    # mcmc_examples/lighthouse/lighthouse_fixed_point_model.py

    def test_mcmc_examples_lighthouse_lighthouse(self):
        self.check_script("mcmc_examples/lighthouse/lighthouse.py")
