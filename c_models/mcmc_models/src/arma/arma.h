/*
 * Copyright (c) 2016 The University of Manchester
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
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
