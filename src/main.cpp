#include <iostream>
#include <cmath>

#include "MOC_Grid_BDE/MOC_GridCalc_BDE.h"

int main() {
    MOC_GridCalc* gridCalc = new MOC_GridCalc();

    double throatPressurePSI = 11.078 * 14.5038;     // CEA t_p  
    double throatTempRankine = 2270.7 * 1.8;         // CEA t_t
    double gamma = 1.2545;                            // CEA t_gamma
    double R = 1545.0 / 21.65;                        // ft-lbf/(lbm-R)
    double throatVelocity = sqrt(gamma * R * throatTempRankine);

    std::cout << "Throat velocity: " << throatVelocity << std::endl;

    gridCalc->SetInitialProperties(
        throatPressurePSI,
        throatTempRankine,
        21.65,               // CEA t_mw
        gamma,
        1.0 * 14.5038,       // ambient
        141,
        1.0,                 // rwtu
        1.0,                 // rwtd
        0.05,
        5,
        30,
        30,
        throatVelocity,
        1,                   // throatFlag = 1
        213.5                // CEA isp
    );

    gridCalc->SetSolutionParameters(
        AXI,
        RAO,
        EPS,
        3.3263,
        15
    );

    gridCalc->SetPrintMode(1);

    int out = gridCalc->CreateMOCGrid();
    std::cout << out << std::endl;
    return 0;
}
