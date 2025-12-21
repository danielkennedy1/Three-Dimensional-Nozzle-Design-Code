#include <gtest/gtest.h>
#include <cmath>
#include "MOC_Grid_BDE/MOC_GridCalc_BDE.h"
#include "MOC_2D/MOC_GridCalc.h"

// Helper function to set up identical test parameters
struct NozzleTestCase {
    double throatPressurePSI;
    double throatTempRankine;
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
    tc.throatPressurePSI = 11.078 * 14.5038;     // CEA t_p
    tc.throatTempRankine = 2270.7 * 1.8;         // CEA t_t
    tc.gamma = 1.2545;                            // CEA t_gamma
    tc.R = 1545.0 / 21.65;                        // ft-lbf/(lbm-R)
    tc.molecularWeight = 21.65;                   // CEA t_mw
    tc.ambientPressure = 1.0 * 14.5038;          // ambient
    tc.nCharacteristics = 141;
    tc.rwtu = 1.0;
    tc.rwtd = 1.0;
    tc.dtLimit = 0.05;
    tc.nRRCAboveBD = 5;
    tc.nSLi = 30;
    tc.nSLj = 30;
    tc.throatFlag = 1;                            // throat conditions
    tc.idealIsp = 213.5;                          // CEA isp

    tc.designParamValue = 3.3263;
    tc.thetaBi = 15;

    return tc;
}

// Regression test: Compare legacy MOC_Grid_BDE vs new MOC_2D
TEST(MOCRegressionTest, RaoAxiNozzle_IdenticalResults) {
    NozzleTestCase tc = GetRaoAxiCase();
    double throatVelocity = sqrt(tc.gamma * tc.R * tc.throatTempRankine);

    // ===== LEGACY VERSION =====
    legacy::MOC_GridCalc* legacy = new legacy::MOC_GridCalc();
    legacy->SetInitialProperties(
        tc.throatPressurePSI,
        tc.throatTempRankine,
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
        tc.throatFlag,
        tc.idealIsp
    );
    legacy->SetSolutionParameters(
        legacy::nozzleGeom::AXI,
        legacy::nozzleType::RAO,
        legacy::param::EPS,
        tc.designParamValue,
        tc.thetaBi
    );
    legacy->SetPrintMode(0);  // Suppress output during testing
    int legacyResult = legacy->CreateMOCGrid();

    // ===== NEW VERSION =====
    moc_2d::MOC_2D* v2 = new moc_2d::MOC_2D();
    v2->SetInitialProperties(
        tc.throatPressurePSI,
        tc.throatTempRankine,
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
        tc.throatFlag,
        tc.idealIsp
    );
    v2->SetSolutionParameters(
        moc_2d::nozzleGeom::AXI,
        moc_2d::nozzleType::RAO,
        moc_2d::param::EPS,
        tc.designParamValue,
        tc.thetaBi
    );
    v2->SetPrintMode(0);  // Suppress output during testing
    int v2Result = v2->CreateMOCGrid();

    // ===== VERIFICATION =====
    // Both should complete successfully
    EXPECT_EQ(legacyResult, v2Result) << "Return codes should match";
    EXPECT_EQ(legacyResult, 1) << "Both should complete successfully";

    // TODO: Add more detailed comparisons once we expose grid data
    // For now, we verify both run to completion with same result code
    // Future additions:
    // - Compare grid points (x, r coordinates)
    // - Compare flow properties (Mach, P, T, rho)
    // - Compare performance metrics (thrust, Isp)

    delete legacy;
    delete v2;
}

// Test that both versions handle errors the same way
TEST(MOCRegressionTest, InvalidInput_SameBehavior) {
    // Test with invalid gamma (should fail)
    legacy::MOC_GridCalc* legacy = new legacy::MOC_GridCalc();
    moc_2d::MOC_2D* v2 = new moc_2d::MOC_2D();

    int legacyResult = legacy->SetInitialProperties(
        100.0,   // pressure
        1000.0,  // temp
        21.65,   // MW
        -1.0,    // INVALID gamma (negative)
        14.5,    // ambient
        141,     // n
        1.0, 1.0, 0.05, 5, 30, 30,
        1000.0,  // velocity
        1,       // throatFlag
        200.0    // Isp
    );

    int v2Result = v2->SetInitialProperties(
        100.0, 1000.0, 21.65, -1.0, 14.5, 141,
        1.0, 1.0, 0.05, 5, 30, 30,
        1000.0, 1, 200.0
    );

    // Both should return 0 (failure) for invalid input
    EXPECT_EQ(legacyResult, v2Result) << "Error handling should match";
    EXPECT_EQ(legacyResult, 0) << "Both should reject invalid gamma";

    delete legacy;
    delete v2;
}

int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
