#include <iostream>
#include <fstream>
#include <chrono>
#include <stdint.h>
#include <math.h>

#include "maths.h"
#include "position_estimator.h"

int main(int argc, char *argv[])
{
    const float dT = 0.02f;
    const float f_a = 0.0f;
    const float f_cut = 0.1f;
    const uint8_t nd_pos = 3;
    positionEstimator_t positionEstimator;
    positionEstimatorUpdateGain(&positionEstimator, f_cut, f_a, dT);
    positionEstimatorInit(&positionEstimator, 0.0f, 0.0f, 0.0f, nd_pos);

    std::cout << " --- G ---" << std::endl;
    // positionEstimator->position = position_k + positionEstimator->a12 * velocity_k                                      + positionEstimator->dT * positionEstimator->a12 * u1 + positionEstimator->k1 * u2;
    // positionEstimator->velocity =              positionEstimator->a22 * velocity_k                                      + positionEstimator->dT * positionEstimator->a22 * u1 + positionEstimator->k2 * u2;
    // positionEstimator->accBias =               positionEstimator->a32 * velocity_k + positionEstimator->a33 * accBias_k + positionEstimator->dT * positionEstimator->a32 * u1 + positionEstimator->k3 * u2;
    std::cout << 1 << ", " << positionEstimator.a12 << ", " <<                     0 << ", " << positionEstimator.dT * positionEstimator.a12 << ", " << positionEstimator.k1 << std::endl;
    std::cout << 0 << ", " << positionEstimator.a22 << ", " <<                     0 << ", " << positionEstimator.dT * positionEstimator.a22 << ", " << positionEstimator.k2 << std::endl;
    std::cout << 0 << ", " << positionEstimator.a32 << ", " << positionEstimator.a33 << ", " << positionEstimator.dT * positionEstimator.a32 << ", " << positionEstimator.k3 << std::endl;
    
    float acc[10][3] = {-0.001053720712662,  0.005268692970276,  9.627979278564453,
                        -0.001053720712662, -0.004311382770538,  9.599238395690918,
                        -0.010633796453476, -0.013891458511353,  9.580078125000000,
                        -0.001053720712662, -0.023471534252167,  9.580078125000000,
                        -0.001053720712662, -0.013891458511353,  9.618398666381836,
                        -0.001053720712662, -0.004311382770538,  9.599238395690918,
                         0.008526355028152, -0.004311382770538,  9.608819007873535,
                         0.008526355028152, -0.004311382770538,  9.589658737182617,
                        -0.010633796453476,  0.014848768711090,  9.599238395690918,
                        -0.010633796453476, -0.013891458511353,  9.599238395690918};

    float baro[10] = {-0.004699736042421,
                      -0.032584879266236,
                      -0.081586452659935,
                      -0.109091273965032,
                      -0.120680865284535,
                      -0.135303794889971,
                      -0.175326437075895,
                      -0.173230376064771,
                      -0.087376019283139,
                      -0.033925864942545};

    float est_rpy[10][3] = { 0.041758361476241, -0.058642777730711,  0.060613539972110,
                             0.008951534255175,  0.005222202162258,  0.004050674760947,
                             0.002768863396341,  0.015821724446141, -0.027824047720060,
                            -0.003265461145929, -0.021485197066795,  0.005150544438948,
                             0.017439821021981, -0.058641417126637, -0.028308715627645,
                             0.031977062462829, -0.025539040507283, -0.053970084991306,
                             0.145460930070840, -0.058847843320109, -0.096178060630336,
                             0.218763307202607, -0.156554131535813, -0.158050315803848,
                             0.388234155252576, -0.295557547360659, -0.197181303519756,
                             0.495851389132440, -0.514213752467185, -0.210484111448750};
    for (size_t i = 0; i < 10; i++) {
        for (size_t j = 0; j < 3; j++) {
            est_rpy[i][j] = est_rpy[i][j] * 1.0e-3f;
        }
    }

    uint32_t cntr = 0;

    // open file for writing
    std::ofstream datafile("output/position_estimator_00.txt");

    std::cout << " --- states ---" << std::endl;

    bool main_execute = true;
    while (main_execute) {
        std::chrono::steady_clock::time_point time_begin = std::chrono::steady_clock::now();

        // CEB = [cos(psi)*cos(theta), cos(psi)*sin(phi)*sin(theta) - cos(phi)*sin(psi), sin(phi)*sin(psi) + cos(phi)*cos(psi)*sin(theta)]
        //       [cos(theta)*sin(psi), cos(phi)*cos(psi) + sin(phi)*sin(psi)*sin(theta), cos(phi)*sin(psi)*sin(theta) - cos(psi)*sin(phi)]
        //       [        -sin(theta),                              cos(theta)*sin(phi),                              cos(phi)*cos(theta)]
        // remove acc bias from accZ in body frame and transform to accZ w.r.t. earth frame
        float accZ = -sin_approx(est_rpy[cntr][1]) * acc[cntr][0] 
            + cos_approx(est_rpy[cntr][1]) * sin_approx(est_rpy[cntr][0]) * acc[cntr][1]
            + cos_approx(est_rpy[cntr][0]) * cos_approx(est_rpy[cntr][1]) * (acc[cntr][2] - positionEstimator.a33 * positionEstimator.accBias);
        // remove gravity to get the acceleration that acts in positive z direction w.r.t earth frame
        accZ -= 9.81f;
        positionEstimatorApply(&positionEstimator, accZ, baro[cntr]);

        std::chrono::steady_clock::time_point time_end = std::chrono::steady_clock::now();
        int64_t time_elapsed_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(time_end - time_begin).count();
        
        if (cntr++ < 10)
        {
            std::cout << cntr << ", " << positionEstimator.position << ", " << positionEstimator.velocity << ", " << positionEstimator.accBias << ", " << time_elapsed_ns << std::endl;
            datafile << cntr << ", " << positionEstimator.position << ", " << positionEstimator.velocity << ", " << positionEstimator.accBias << ", " << time_elapsed_ns << std::endl;
        }
        else
        {
            main_execute = false;
        }
    }

    datafile.close();

    return 0;
}
