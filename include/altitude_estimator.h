/*
 * This file is part of Cleanflight and Betaflight.
 *
 * Cleanflight and Betaflight are free software. You can redistribute
 * this software and/or modify this software under the terms of the
 * GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * Cleanflight and Betaflight are distributed in the hope that they
 * will be useful, but WITHOUT ANY WARRANTY; without even the implied
 * warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this software.
 *
 * If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once

#include <stdint.h>

#define POS_PAST_MAXLEN 41 // posDiscreteDelay + 1 <= posPast_MAX

typedef struct altitudeEstimator_s {
    float Ts;   // sampling time in sec
    float x[3]; // state (pos, vel, acc bias)
    float K[3]; // estimator / observer gain
    float wa;   // acc bias dynamics, set this to zero for a pure integrator for the state
    uint8_t posDiscreteDelay;       // discrete delay of pos
    float posPast[POS_PAST_MAXLEN]; // past pos values
    uint8_t posPastIndex;           // past pos index
} altitudeEstimator_t;

void altitudeEstimatorSetup(altitudeEstimator_t *altitudeEstimator, float Ts, float w1, float w2, float D2, float wa, uint8_t posDiscreteDelay);
void altitudeEstimatorInitialize(altitudeEstimator_t *altitudeEstimator, float pos, float vel, float accBias);
void altitudeEstimatorUpdate(altitudeEstimator_t *altitudeEstimator, float accX, float accY, float accZ, float pos, float phi, float theta);