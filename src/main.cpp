#include <cmath>
#include <iostream>

#include "MOC_Grid_BDE/MOC_GridCalc_BDE.h"

int main() {
    legacy::MOC_GridCalc* gridCalc = new legacy::MOC_GridCalc();

    double chamber_pressure_psi = 20.0 * 14.5038;        // CEA t_p
    double chamber_temp_rankine = 2550.7 * 1.8;            // CEA t_t
    double gamma = 1.27;                              // CEA t_gamma
    double exit_pressure_psi = 1.1 * 14.5038;


    gridCalc->SetInitialProperties(
        chamber_pressure_psi,
        chamber_temp_rankine,
        21.65,              // CEA t_mw
        gamma,
        exit_pressure_psi,  // ambient
        25,
        1.0,                // rwtu
        1.0,                // rwtd
        0.1,
        5,
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
        10
    );

    gridCalc->SetPrintMode(1);

    int out = gridCalc->CreateMOCGrid();
    std::cout << out << std::endl;
    return 0;
}
