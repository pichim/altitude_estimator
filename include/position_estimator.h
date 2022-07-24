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

#define POSITIONPAST_MAXLEN 41 // positionDiscreteDelay + 1 <= POSITIONPAST_MAXLEN

typedef struct positionEstimator_s {
    float state[3]; // position (pos), velocity (vel) and acc bias (accBias)
    float K[3];     // time continous estimator gain
    float wa;       // acc bias dynamics in rad/sec, set this to zero for a pure integrator for the state
    uint8_t positionDiscreteDelay;           // discrete delay of position
    float positionPast[POSITIONPAST_MAXLEN]; // past position values
    uint8_t positionPastIndex;               // past position index
} positionEstimator_t;

void positionEstimatorSetupPt3(positionEstimator_t *positionEstimator, float f_cut);
void positionEstimatorUpdateCutoff(positionEstimator_t *positionEstimator, float filterFreq1, float filterFreq2, float Q2, float fa);
void positionEstimatorInit(positionEstimator_t *positionEstimator, float pos, float vel, float accBias, uint8_t positionDiscreteDelay);
void positionEstimatorApply(positionEstimator_t *positionEstimator, float acc, float pos, float dT);