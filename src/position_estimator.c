#include "position_estimator.h"

#include <math.h>
#include "maths.h"

void positionEstimatorSetupPt3(positionEstimator_t *positionEstimator, float f_cut)
{
    const float order = 3.0f;
    const float orderCutoffCorrection = 1.0f / sqrtf(powf(2.0f, 1.0f / order) - 1.0f);
    const float Q2 = 0.5f;
    const float fa = 0.0f;
    const uint8_t positionDiscreteDelay = 0;
    float filterFreq = f_cut * orderCutoffCorrection;
    positionEstimatorUpdateCutoff(positionEstimator, filterFreq, filterFreq, Q2, fa);
    positionEstimatorInit(positionEstimator, 0.0f, 0.0f, 0.0f, positionDiscreteDelay);
}

void positionEstimatorUpdateCutoff(positionEstimator_t *positionEstimator, float filterFreq1, float filterFreq2, float Q2, float fa)
{
    float wa = 2.0f * M_PIf * fa;
    positionEstimator->wa = wa;

    // continous time gain
    float w1 = 2.0f * M_PIf * filterFreq1;
    float w2 = 2.0f * M_PIf * filterFreq2;
    float D2 = 1.0f / ( 2.0f * Q2 );
    positionEstimator->K[0] = (w1 + 2.0f*D2*w2) - wa;
    positionEstimator->K[1] = (w2*w2 + 2.0f*D2*w1*w2) - positionEstimator->K[0]*wa;
    positionEstimator->K[2] = positionEstimator->K[1]*wa - w1*w2*w2;
}

void positionEstimatorInit(positionEstimator_t *positionEstimator, float pos, float vel, float accBias, uint8_t positionDiscreteDelay)
{   
    // set state
    positionEstimator->state[0] = pos;
    positionEstimator->state[1] = vel;
    positionEstimator->state[2] = accBias;

    // set delayed position estimate
    positionEstimator->positionDiscreteDelay = positionDiscreteDelay;
    for (uint8_t i = 0; i < positionEstimator->positionDiscreteDelay + 1; i++) {
        positionEstimator->positionPast[i] = pos;
    }
    positionEstimator->positionPastIndex = 0;
}

void positionEstimatorApply(positionEstimator_t *positionEstimator, float acc, float pos, float dT)
{
    /* predict position (pos), velocity (vel) and acc bias (accBias)
       pos = pos + dT * ( vel + 0.5 * acc * Ts );
       vel = vel + dT * acc;
       accBias = accBias - dT*wa*accBias; */
    positionEstimator->state[0] += dT * ( positionEstimator->state[1] + 0.5f * acc * dT );
    positionEstimator->state[1] += dT * acc;
    positionEstimator->state[2] += dT * -positionEstimator->wa*positionEstimator->state[2];

    // update delayed position estimation
    positionEstimator->positionPast[positionEstimator->positionPastIndex] = positionEstimator->state[0];
    positionEstimator->positionPastIndex++;
    if (positionEstimator->positionPastIndex > positionEstimator->positionDiscreteDelay) {
        positionEstimator->positionPastIndex = 0;
    }

    // correction based on delayed position estimate
    float e = pos - positionEstimator->positionPast[positionEstimator->positionPastIndex];
    for (uint8_t i = 0; i < 3; i++) {
        positionEstimator->state[i] += positionEstimator->K[i] * e * dT;
    }
}