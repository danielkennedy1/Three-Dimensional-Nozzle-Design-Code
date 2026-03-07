#include <cmath>
#include <iostream>
#include <chrono>

#include "MOC_Grid_BDE/MOC_GridCalc_BDE.h"

using std::chrono::high_resolution_clock;
using std::chrono::duration_cast;
using std::chrono::duration;
using std::chrono::milliseconds;

int main() {
    auto t1 = high_resolution_clock::now();
    legacy::MOC_GridCalc* gridCalc = new legacy::MOC_GridCalc();

    double chamber_pressure_psi = 20.0 * 14.5038;        // CEA t_p
    double chamber_temp_rankine = 2550.7 * 1.8;            // CEA t_t
    double gamma = 1.27;                              // CEA t_gamma
    double exit_pressure_psi = 1.01 * 14.5038;


    gridCalc->SetInitialProperties(
        chamber_pressure_psi,
        chamber_temp_rankine,
        21.65,              // CEA t_mw
        gamma,
        exit_pressure_psi,  // ambient
        141,
        1.5,                // rwtu
        0.8,                // rwtd
        1.0,
        15,
        30,
        30,
        0,
        0,
        213.5               // CEA isp
    );

    gridCalc->SetSolutionParameters(
        legacy::AXI,
        legacy::RAO,
        legacy::EXITPRESSURE,
        exit_pressure_psi,
        15
    );

    gridCalc->SetPrintMode(1);

    int out = gridCalc->CreateMOCGrid();
    std::cout << out << std::endl;

    auto t2 = high_resolution_clock::now();

    duration<double, std::milli> elapsed = t2 - t1;

    std::cout << elapsed.count() << "ms\n";
}
