/*
 * Copyright (c) 2016 The University of Manchester
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     https://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */
struct mcmc_params {

    // scaling of t transition distribution for MH jumps
    CALC_TYPE alpha_jump_scale;
    CALC_TYPE beta_jump_scale;

    // Alpha range
    CALC_TYPE alpha_min;
    CALC_TYPE alpha_max;

    // Beta range
    CALC_TYPE beta_min;
    CALC_TYPE beta_max;
};

struct mcmc_state {
    CALC_TYPE alpha;
    CALC_TYPE beta;
};
