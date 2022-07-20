#include <iostream>
#include <fstream>

#include <chrono>

#include "chirp_cpp.h"
#include "chirp_c.h"
#include "filter.h"

int main(int argc, char *argv[])
{
    bool do_plot_in_terminal = false;
    bool do_use_cpp_implementation = false; // 19.02.2022 stopped also taking cpp implementation further
    bool do_use_filter = true;
    bool do_use_exp_chirp = true;

    // set up chirp generator
    float Ts = 1.0f / 8.0e3f;
    float f0 = 0.2f;
    float f1 = 3950.0f;
    uint32_t N = 30 * 8e3;
    float alpha = 7.0;

    CHIRP *chirp_cpp = new CHIRP();
    chirp_cpp->init(f0, f1, N, Ts);

    // chirp_t *chirp_c = new chirp_t;
    // chirpInit(chirp_c, f0, f1, N, Ts);
    chirp_t chirp_c;
    if (do_use_exp_chirp) {
        chirpInitExp(&chirp_c, f0, f1, N, Ts);
    } else {
        chirpInitAlpha(&chirp_c, f0, f1, N, Ts, alpha);
    }    

    // float f_cut = 10.0f;
    // biquadFilter_t *filter_c = new biquadFilter_t;
    // biquadFilterInitLPF(filter_c, f_cut, (uint32_t)(Ts * 1.0e6f));
    // pt1TustinFilter_t *filter_c = new pt1TustinFilter_t;
    // pt1TustinFilterInit(filter_c, f_cut, Ts);
    // pt1TustinFilter_t *filter_c2 = new pt1TustinFilter_t;
    // pt1TustinFilterInit(filter_c2, 100.0f * f_cut, Ts);
    // pt2Filter_t *filter_c = new pt2Filter_t;
    // pt2FilterInit(filter_c, pt2FilterGain(f_cut, Ts));
    leadlag1Filter_t *filter_c = new leadlag1Filter_t;
    leadlag1FilterInit(filter_c, 30.0f, 3.0f, Ts);

    uint32_t cntr = 0;

    // open file for writing
    std::ofstream datafile("output/chirp_meas_01.txt");

    bool main_execute = true;
    while (main_execute) {
        cntr++; 
        if (do_use_cpp_implementation) {
            std::chrono::steady_clock::time_point time_begin = std::chrono::steady_clock::now();
            bool wasChripUpdate =  chirp_cpp->update();
            std::chrono::steady_clock::time_point time_end = std::chrono::steady_clock::now();
            int64_t time_elapsed_ns = std::chrono::duration_cast<std::chrono::nanoseconds> (time_end - time_begin).count();
            if (wasChripUpdate) {
                float exc = chirp_cpp->exc();
                float excf = 0.0f;
                if (do_use_filter) {
                    // excf = biquadFilterApply(filter_c, exc);
                    // excf = pt1TustinFilterApply(filter_c, exc);
                    // excf = pt1TustinFilterApply(filter_c2, excf);
                    // excf = pt2FilterApply(filter_c, exc);
                    excf = leadlag1FilterApply(filter_c, exc);
                }
                if (do_plot_in_terminal) {
                    std::cout << cntr << ", " << exc << ", " << excf << ", " << chirp_cpp->sinarg() << ", " << chirp_cpp->fchirp() << ", " << time_elapsed_ns << std::endl;
                }
                datafile << cntr << ", " << exc << ", " << excf << ", " << chirp_cpp->sinarg() << ", " << chirp_cpp->fchirp() << ", " << time_elapsed_ns << std::endl;
            } else {
                std::cout << "---   chirp finished" << std::endl;
                main_execute = false;
            }
        } else {
            std::chrono::steady_clock::time_point time_begin = std::chrono::steady_clock::now();
            bool wasChripUpdate =  chirpUpdate(&chirp_c);
            std::chrono::steady_clock::time_point time_end = std::chrono::steady_clock::now();
            int64_t time_elapsed_ns = std::chrono::duration_cast<std::chrono::nanoseconds> (time_end - time_begin).count();
            if (wasChripUpdate) {
                float exc = chirp_c.exc;
                float excf = 0.0f;
                if (do_use_filter) {
                    // excf = biquadFilterApply(filter_c, exc);
                    // excf = pt1TustinFilterApply(filter_c, exc);
                    // excf = pt1TustinFilterApply(filter_c2, excf);
                    // excf = pt2FilterApply(filter_c, exc);
                    excf = leadlag1FilterApply(filter_c, exc);
                }
                if (do_plot_in_terminal) {
                    std::cout << cntr << ", " << exc << ", " << excf << ", " << chirp_c.sinarg << ", " << chirp_c.fchirp << ", " << time_elapsed_ns << std::endl;
                }
                datafile << cntr << ", " << exc << ", " << excf << ", " << chirp_c.sinarg << ", " << chirp_c.fchirp << ", " << time_elapsed_ns << std::endl;
            } else {
                std::cout << "---   chirp finished" << std::endl;
                main_execute = false;
            }
        }
    }

    datafile.close();

    return 0;
}
