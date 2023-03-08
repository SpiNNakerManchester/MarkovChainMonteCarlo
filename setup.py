# Copyright (c) 2018 The University of Manchester
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

from setuptools import setup
from collections import defaultdict
import os

__version__ = None
exec(open("mcmc/_version.py").read())
assert __version__

# Build a list of all project modules, as well as supplementary files
main_package = "spinnaker_graph_front_end"
extensions = {".aplx", ".boot", ".cfg", ".json", ".sql", ".template", ".xml",
              ".xsd"}
main_package_dir = os.path.join(os.path.dirname(__file__), main_package)
start = len(main_package_dir)
packages = []
package_data = defaultdict(list)
for dirname, dirnames, filenames in os.walk(main_package_dir):
    if '__init__.py' in filenames:
        package = "{}{}".format(
            main_package, dirname[start:].replace(os.sep, '.'))
        packages.append(package)
    for filename in filenames:
        _, ext = os.path.splitext(filename)
        if ext in extensions:
            package = "{}{}".format(
                main_package, dirname[start:].replace(os.sep, '.'))
            package_data[package].append(filename)

setup(
    name="MarkovChainMonteCarlo",
    version=__version__,
    description="Front end to MCMC methods implemented on SpiNNaker",
    url="https://github.com/SpiNNakerManchester/MarkovChainMonteCarlo",
    packages=packages,
    package_data=package_data,
    install_requires=['SpiNNUtilities == 1!6.0.1',
                      'SpiNNMachine == 1!6.0.1',
                      'SpiNNMan == 1!6.0.1',
                      'SpiNNaker_PACMAN == 1!6.0.1',
                      'SpiNNaker_DataSpecification == 1!6.0.1',
                      'SpiNNFrontEndCommon == 1!6.0.1',
                      'SpiNNakerGraphFrontEnd == 1!6.0.1']
)
