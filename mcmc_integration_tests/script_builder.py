# Copyright (c) 2017-2019 The University of Manchester
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

from spinnaker_testbase import RootScriptBuilder


class ScriptBuilder(RootScriptBuilder):
    """
    This file will recreate the test_scripts.py file
    """

    def build_mcmc_scripts(self):
        # create_test_scripts supports test that are too long or exceptions
        # These scripts raise a SkipTest with the reasons given
        exceptions = {}
        NOT_SCRIPT = "Not a run script"
        exceptions["arma_fixed_point_model.py"] = NOT_SCRIPT
        exceptions["arma_float_model.py"] = NOT_SCRIPT
        exceptions["arma_model.py"] = NOT_SCRIPT
        exceptions["lighthouse_fixed_point_model.py"] = NOT_SCRIPT
        exceptions["lighthouse_float_model.py"] = NOT_SCRIPT
        exceptions["lighthouse_model.py"] = NOT_SCRIPT
        self.create_test_scripts(["mcmc_examples/"], exceptions=exceptions)


if __name__ == '__main__':
    builder = ScriptBuilder()
    builder.build_mcmc_scripts()
