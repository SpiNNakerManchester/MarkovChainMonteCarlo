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
#include <mcmc_model.h>

#define PPOLYORDER 9
#define QPOLYORDER 9

struct mcmc_params {

    // scaling of transition distribution for MH jumps
	CALC_TYPE jump_scale[PPOLYORDER+QPOLYORDER+2];

};

struct mcmc_state {
    // parameters has size order_p (p poly) + order_q (q poly) + 2 (mu and sigma)
    CALC_TYPE parameters[PPOLYORDER+QPOLYORDER+2];
};
