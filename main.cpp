#include <iostream>
#include <fstream>
#include <chrono>
#include <stdint.h>

#include "altitude_estimator.h"

int main(int argc, char *argv[])
{
    float Ts = 0.02f;
    float w1 = 1.2566371f;
    float w2 = 1.2566371f;
    float D2 = 1.0000000f;
    float wa = 0.0f;
    float posDiscreteDelay = 5;
    altitudeEstimator_t altitudeEstimator;
    altitudeEstimatorSetup(&altitudeEstimator, Ts, w1, w2, D2, wa, posDiscreteDelay);

    float acc[10][3] = {-0.000670582056046, 0.000957667827606, 9.695039749145508,
                        0.018489569425583, 0.010537743568420, 9.647139549255371,
                        0.008909493684769, -0.008622407913208, 9.637558937072754,
                        0.018489569425583, 0.010537743568420, 9.637558937072754,
                        0.008909493684769, 0.010537743568420, 9.637558937072754,
                        0.018489569425583, 0.010537743568420, 9.637558937072754,
                        0.008909493684769, 0.000957667827606, 9.647139549255371,
                        0.008909493684769, 0.010537743568420, 9.637558937072754,
                        -0.000670582056046, -0.008622407913208, 9.647139549255371,
                        0.008909493684769, 0.000957667827606, 9.627979278564453};

    float baro[10] = {-0.416809082031250,
                      -0.637451171875000,
                      -0.637451171875000,
                      -0.504852294921875,
                      -0.504852294921875,
                      -0.505126953125000,
                      -0.505126953125000,
                      -0.681701660156250,
                      -0.681701660156250,
                      -0.505126953125000};

    float est_rpy[10][3] = {-0.000027377293009, -0.000003089276106, -0.000025992456358,
                            -0.000044101929234, 0.000052410850913, -0.000030681145290,
                            -0.000028863129046, 0.000038669422793, -0.000010648071111,
                            0.000002565495379, 0.000072978378739, 0.000008982236068,
                            -0.000026080982934, 0.000133628986077, -0.000017334623408,
                            0.000006878058684, 0.000116625829833, -0.000005359741408,
                            0.000232226171647, 0.000023925575078, -0.000015440193238,
                            0.000434604473412, -0.000243514848989, -0.000087071479356,
                            0.000716522976290, -0.000686427985784, -0.000106513150968,
                            0.000748570892029, -0.001031362684444, -0.000120513854199};

    uint32_t cntr = 0;

    // open file for writing
    std::ofstream datafile("output/altitude_estimator_00.txt");

    bool main_execute = true;
    while (main_execute) {
        std::chrono::steady_clock::time_point time_begin = std::chrono::steady_clock::now();
        altitudeEstimatorUpdate(&altitudeEstimator, acc[cntr][0], acc[cntr][1], acc[cntr][2], baro[cntr], est_rpy[cntr][0], est_rpy[cntr][1]);
        std::chrono::steady_clock::time_point time_end = std::chrono::steady_clock::now();
        int64_t time_elapsed_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(time_end - time_begin).count();
        //cntr++;
        if (cntr++ < 10)
        {
            std::cout << cntr << ", " << altitudeEstimator.x[0] << ", " << altitudeEstimator.x[1] << ", " << altitudeEstimator.x[2] << ", " << time_elapsed_ns << std::endl;
            datafile << cntr << ", " << altitudeEstimator.x[0] << ", " << altitudeEstimator.x[1] << ", " << altitudeEstimator.x[2] << ", " << time_elapsed_ns << std::endl;
        }
        else
        {
            main_execute = false;
        }
    }

    datafile.close();

    return 0;
}
