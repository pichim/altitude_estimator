#include "altitude_estimator.h"

#include "maths.h"

void altitudeEstimatorSetup(altitudeEstimator_t *altitudeEstimator, float Ts, float w1, float w2, float D2, float wa, uint8_t posDiscreteDelay)
{
    altitudeEstimator->Ts = Ts;
    altitudeEstimator->wa = wa;

    // continous time gain
    altitudeEstimator->K[0] = (w1 + 2.0f*D2*w2) - wa;
    altitudeEstimator->K[1] = (w2*w2 + 2.0f*D2*w1*w2) - altitudeEstimator->K[0]*wa;
    altitudeEstimator->K[2] = altitudeEstimator->K[1]*wa - w1*w2*w2;

    // discrete time gain
    for (uint8_t i = 0; i < 3; i++) {
        altitudeEstimator->K[i] *= Ts;
    }

    // set state and delayed position estimate to zero
    for (uint8_t i = 0; i < 3; i++) {
        altitudeEstimator->x[i] = 0.0f;
    }
    altitudeEstimator->posDiscreteDelay = posDiscreteDelay;
    for (uint8_t i = 0; i < posDiscreteDelay + 1; i++) {
        altitudeEstimator->posPast[i] = 0.0f;
    }
    altitudeEstimator->posPastIndex = 0;
}

void altitudeEstimatorInitialize(altitudeEstimator_t *altitudeEstimator, float pos, float vel, float accBias)
{   
    // set state
    altitudeEstimator->x[0] = pos;
    altitudeEstimator->x[1] = vel;
    altitudeEstimator->x[2] = accBias;

    // set delayed position estimate
    for (uint8_t i = 0; i < altitudeEstimator->posDiscreteDelay + 1; i++) {
        altitudeEstimator->posPast[i] = pos;
    }
    altitudeEstimator->posPastIndex = 0;
}

void altitudeEstimatorUpdate(altitudeEstimator_t *altitudeEstimator, float accX, float accY, float accZ, float pos, float phi, float theta)
{
    /* CEB =
    [cos(psi)*cos(theta), cos(psi)*sin(phi)*sin(theta) - cos(phi)*sin(psi), sin(phi)*sin(psi) + cos(phi)*cos(psi)*sin(theta)]
    [cos(theta)*sin(psi), cos(phi)*cos(psi) + sin(phi)*sin(psi)*sin(theta), cos(phi)*sin(psi)*sin(theta) - cos(psi)*sin(phi)]
    [        -sin(theta),                              cos(theta)*sin(phi),                              cos(phi)*cos(theta)]
    */

    // remove acc bias from accZ in body frame and transform to accZ w.r.t. earth frame
    float u = -sin_approx(theta)*accX + cos_approx(theta)*sin_approx(phi)*accY + cos_approx(phi)*cos_approx(theta)*( accZ - altitudeEstimator->x[2] );
    // remove gravity to get the acceleration that acts in positive z direction
    u -= 9.81f;

    /* predict pos (p), vel (v) and acc bias (a)
       p = p + v*Ts + 0.5*a*Ts^2;
       v = v + a*Ts;
       b = b - Ts*wa*b; */
    altitudeEstimator->x[0] += altitudeEstimator->Ts*(altitudeEstimator->x[1] + 0.5f*u*altitudeEstimator->Ts);
    altitudeEstimator->x[1] += altitudeEstimator->Ts*u;
    altitudeEstimator->x[2] += altitudeEstimator->Ts*-altitudeEstimator->wa*altitudeEstimator->x[2];

    // update delayed position estimation
    altitudeEstimator->posPast[altitudeEstimator->posPastIndex] = altitudeEstimator->x[0];
    altitudeEstimator->posPastIndex++;
    if (altitudeEstimator->posPastIndex > altitudeEstimator->posDiscreteDelay) {
        altitudeEstimator->posPastIndex = 0;
    }

    // correction based on delayed position estimate
    float e = pos - altitudeEstimator->posPast[altitudeEstimator->posPastIndex];
    for (uint8_t i = 0; i < 3; i++) {
        altitudeEstimator->x[i] += altitudeEstimator->K[i] * e;
    }
}