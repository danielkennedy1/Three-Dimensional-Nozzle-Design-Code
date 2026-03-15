#include <gtest/gtest.h>
#include <cmath>
#include "MOC_Grid_BDE/MOC_GridCalc_BDE.h"

const double GASCON = 1545; // Universal Gas Constant ft-lbf/lbm-mol-R
const double GRAV			=	32.174;

struct NozzleTestCase {
    double pTotalPSI;
    double tempRankine;
    double gamma;
    double R;
    double molecularWeight;
    double ambientPressure;
    int nCharacteristics;
    double rwtu;
    double rwtd;
    double dtLimit;
    int nRRCAboveBD;
    int nSLi;
    int nSLj;
    int throatFlag;
    double idealIsp;

    double designParamValue;
    double thetaBi;
};

// Test case from main.cpp - RAO axisymmetric nozzle
NozzleTestCase GetRaoAxiCase() {
    NozzleTestCase tc;
    tc.pTotalPSI = 20 * 14.5038;                    // CEA t_p
    tc.tempRankine = 2550.0 * 1.8;                  // CEA t_t
    tc.gamma = 1.2545;                              // CEA t_gamma
    tc.R = 1545.0 / 21.65;                          // ft-lbf/(lbm-R)
    tc.molecularWeight = 21.65;                     // CEA t_mw
    tc.ambientPressure = 1.0 * 14.5038;             // ambient
    tc.nCharacteristics = 141;
    tc.rwtu = 1.0;
    tc.rwtd = 1.0;
    tc.dtLimit = 0.05;
    tc.nRRCAboveBD = 5;
    tc.nSLi = 30;
    tc.nSLj = 30;
    tc.throatFlag = 0;                              // total conditions
    tc.idealIsp = 213.5;                            // CEA isp

    tc.designParamValue = 3.3263;
    tc.thetaBi = 15;

    return tc;
}

TEST(MOCLegacyTest, ThroatFlag_IsentropicPropertiesConstantAlongInitialLine) {
    NozzleTestCase tc = GetRaoAxiCase();
    double g = tc.gamma;
    double ratio = 1.0 + (g - 1.0) / 2.0;  // M=1
    double staticTemp = tc.tempRankine / ratio;
    double staticPres = tc.pTotalPSI / pow(ratio, g / (g - 1.0));
    double throatVelocity = sqrt(g * GASCON / tc.molecularWeight * GRAV * staticTemp);

    legacy::MOC_GridCalc* leg = new legacy::MOC_GridCalc();
    leg->SetInitialProperties(
        staticPres,
		staticTemp,
		tc.molecularWeight,
		tc.gamma,
        tc.ambientPressure,
		tc.nCharacteristics,
		tc.rwtu,
		tc.rwtd,
        tc.dtLimit,
		tc.nRRCAboveBD,
		tc.nSLi,
		tc.nSLj,
        throatVelocity,
		1,
		tc.idealIsp
    );

    leg->InitializeDataMembers();
    leg->SetPrintMode(1);
    leg->CalcInitialThroatLine(tc.rwtu, leg->nC - 1, tc.gamma, tc.ambientPressure, legacy::nozzleGeom::AXI, 1, 1.0);

    double pRef   = leg->pres[0][0];
    double tRef   = leg->temp[0][0];
    double rhoRef = leg->rho[0][0];

    for (int i = 1; i <= leg->iLast[0]; ++i) {
        EXPECT_DOUBLE_EQ(leg->pres[i][0],  pRef)   << "pres["  << i << "][0]";
        EXPECT_DOUBLE_EQ(leg->temp[i][0],  tRef)   << "temp["  << i << "][0]";
        EXPECT_DOUBLE_EQ(leg->rho[i][0],   rhoRef) << "rho["   << i << "][0]";
    }

    delete leg;
}

int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
