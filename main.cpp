#include <iostream>
#include <fstream>
#include <chrono>
#include <stdint.h>
#include <math.h>

#include "maths.h"
#include "position_estimator.h"

int main(int argc, char *argv[])
{
    float dT = 0.02f;
    float f_a = 0.01f;
    float f_cut = 0.1f;
    uint8_t nd_pos = 3;
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
    
    float acc[10][3] = {-0.000257140562159,  0.001285724628871, 2.349525803562698,
                        -0.000848143575505,  0.001902953988124, 7.742581336210291,
                        -0.003351659043102, -0.004979671157791, 9.235621230364092,
                        -0.006420034714341, -0.014031958667323, 9.521973231718768,
                        -0.002558909315191, -0.018844881751511, 9.580155420965871,
                        -0.001353664220325, -0.012988037832677, 9.604819002434644,
                         0.001231659968945, -0.006468956699479, 9.603649226166686,
                         0.006648759220139, -0.004724804800184, 9.602717604278396,
                         0.003486611537764,  0.000293380126959, 9.595451925103136,
                        -0.006958540400100,  0.004085750011728, 9.598045605586460};

    float baro[10] = {-0.009782689670254,
                      -0.053625560591952,
                      -0.109018738313538,
                      -0.123976722835630,
                      -0.126801293146205,
                      -0.148636878848532,
                      -0.197803083656959,
                      -0.158299006762922,
                      -0.038591423628320,
                      -0.004736179887408};

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
