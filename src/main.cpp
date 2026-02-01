#include <iostream>
#include <cmath>

#include "MOC_Grid_BDE/MOC_GridCalc_BDE.h"

int main() {
    legacy::MOC_GridCalc* gridCalc = new legacy::MOC_GridCalc();

    double presPSI = 20.0 * 14.5038;        // CEA t_p
    double throatTempRankine = 2550.7 * 1.8;            // CEA t_t
    double gamma = 1.2545;                              // CEA t_gamma
    double ambient = 1.0 * 14.5038;


    gridCalc->SetInitialProperties(
        presPSI,
        throatTempRankine,
        21.65,              // CEA t_mw
        gamma,
        ambient,            // ambient
        141,
        1.0,                // rwtu
        1.0,                // rwtd
        0.05,
        5,
        30,
        30,
        0,
        0,
        213.5               // CEA isp
    );

    gridCalc->SetSolutionParameters(
        legacy::AXI,
        legacy::PERFECT,
        legacy::EXITPRESSURE,
        ambient * 1.1,
        15
    );

    gridCalc->SetPrintMode(1);

    int out = gridCalc->CreateMOCGrid();
    std::cout << out << std::endl;
    return 0;
}
