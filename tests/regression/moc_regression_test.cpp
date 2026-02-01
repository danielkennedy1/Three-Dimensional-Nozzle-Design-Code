#include <gtest/gtest.h>
#include <cmath>
#include "MOC_Grid_BDE/MOC_GridCalc_BDE.h"
#include "MOC_2D/MOC_GridCalc.h"

// Helper function to set up identical test parameters
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
    tc.pTotalPSI = 20 * 14.5038;     // CEA t_p
    tc.tempRankine = 2550.0 * 1.8;         // CEA t_t
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
    tc.throatFlag = 0; // total conditions
    tc.idealIsp = 213.5;                          // CEA isp

    tc.designParamValue = 3.3263;
    tc.thetaBi = 15;

    return tc;
}

// Regression test: Compare legacy MOC_Grid_BDE vs new MOC_2D
TEST(MOCRegressionTest, RaoAxiNozzle_IdenticalResults) {
    NozzleTestCase tc = GetRaoAxiCase();

    double throatVelocity = 0.0;

    // ===== LEGACY VERSION =====
    legacy::MOC_GridCalc* legacy = new legacy::MOC_GridCalc();
    legacy->SetInitialProperties(
        tc.pTotalPSI,
        tc.tempRankine,
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
        legacy::nozzleType::PERFECT,
        legacy::param::EXITPRESSURE,
        tc.ambientPressure * 1.1,
        tc.thetaBi
    );
    legacy->SetPrintMode(0);  // Suppress output during testing
    int legacyResult = legacy->CreateMOCGrid();

    // ===== NEW VERSION =====
    moc_2d::MOC_2D* v2 = new moc_2d::MOC_2D();
    v2->SetInitialProperties(
        tc.pTotalPSI,
        tc.tempRankine,
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
        moc_2d::nozzleType::PERFECT,
        moc_2d::param::EXITPRESSURE,
        tc.ambientPressure * 1.1,
        tc.thetaBi
    );
    v2->SetPrintMode(0);  // Suppress output during testing
    int v2Result = v2->CreateMOCGrid();

    // ===== VERIFICATION =====
    // Both should complete successfully
    EXPECT_EQ(legacyResult, v2Result) << "Return codes should match";
    EXPECT_EQ(legacyResult, 3) << "Both should complete successfully";

    // Only proceed with detailed comparison if both succeeded
    if (legacyResult != 1 || v2Result != 1) {
        delete legacy;
        delete v2;
        return;
    }

    // ===== DETAILED GRID COMPARISON =====
    // Compare grid dimensions
    EXPECT_EQ(legacy->maxLRC, v2->maxLRC) << "maxLRC should match";
    EXPECT_EQ(legacy->maxRRC, v2->maxRRC) << "maxRRC should match";
    EXPECT_EQ(legacy->lastRRC, v2->lastRRC) << "lastRRC should match";
    EXPECT_EQ(legacy->iBD, v2->iBD) << "iBD should match";
    EXPECT_EQ(legacy->jBD, v2->jBD) << "jBD should match";
    EXPECT_EQ(legacy->jDELast, v2->jDELast) << "jDELast should match";

    //// Get grid dimensions for iteration
    int maxLRC = legacy->maxLRC;
    int maxRRC = legacy->maxRRC;
    int jDELast = legacy->jDELast;

    //// Compare iLast array
    const int* iLastLegacy = legacy->iLast;
    const int* iLastV2 = v2->iLast;
    for (int j = 0; j <= maxLRC; ++j) {
        EXPECT_EQ(iLastLegacy[j], iLastV2[j]) << "iLast[" << j << "] should match";
    }

    //// Get 2D grid arrays
    const double* const* machLegacy = legacy->mach;
    const double* const* machV2 = v2->mach;
    const double* const* presLegacy = legacy->pres;
    const double* const* presV2 = v2->pres;
    const double* const* tempLegacy = legacy->temp;
    const double* const* tempV2 = v2->temp;
    const double* const* rhoLegacy = legacy->rho;
    const double* const* rhoV2 = v2->rho;
    const double* const* xLegacy = legacy->x;
    const double* const* xV2 = v2->x;
    const double* const* rLegacy = legacy->r;
    const double* const* rV2 = v2->r;
    const double* const* gammaLegacy = legacy->gamma;
    const double* const* gammaV2 = v2->gamma;
    const double* const* thetaLegacy = legacy->theta;
    const double* const* thetaV2 = v2->theta;
    const double* const* massflowLegacy = legacy->massflow;
    const double* const* massflowV2 = v2->massflow;
    const double* const* thrustLegacy = legacy->thrust;
    const double* const* thrustV2 = v2->thrust;
    const double* const* sthrustLegacy = legacy->Sthrust;
    const double* const* sthrustV2 = v2->Sthrust;

    //// Compare all 2D grid arrays point by point
    for (int j = 0; j <= maxLRC; ++j) {
        for (int i = 0; i <= iLastLegacy[j]; ++i) {
            EXPECT_DOUBLE_EQ(machLegacy[i][j], machV2[i][j])
                << "mach[" << i << "][" << j << "] should match exactly";
            EXPECT_DOUBLE_EQ(presLegacy[i][j], presV2[i][j])
                << "pres[" << i << "][" << j << "] should match exactly";
            EXPECT_DOUBLE_EQ(tempLegacy[i][j], tempV2[i][j])
                << "temp[" << i << "][" << j << "] should match exactly";
            EXPECT_DOUBLE_EQ(rhoLegacy[i][j], rhoV2[i][j])
                << "rho[" << i << "][" << j << "] should match exactly";
            EXPECT_DOUBLE_EQ(xLegacy[i][j], xV2[i][j])
                << "x[" << i << "][" << j << "] should match exactly";
            EXPECT_DOUBLE_EQ(rLegacy[i][j], rV2[i][j])
                << "r[" << i << "][" << j << "] should match exactly";
            EXPECT_DOUBLE_EQ(gammaLegacy[i][j], gammaV2[i][j])
                << "gamma[" << i << "][" << j << "] should match exactly";
            EXPECT_DOUBLE_EQ(thetaLegacy[i][j], thetaV2[i][j])
                << "theta[" << i << "][" << j << "] should match exactly";
            EXPECT_DOUBLE_EQ(massflowLegacy[i][j], massflowV2[i][j])
                << "massflow[" << i << "][" << j << "] should match exactly";
            EXPECT_DOUBLE_EQ(thrustLegacy[i][j], thrustV2[i][j])
                << "thrust[" << i << "][" << j << "] should match exactly";
            EXPECT_DOUBLE_EQ(sthrustLegacy[i][j], sthrustV2[i][j])
                << "Sthrust[" << i << "][" << j << "] should match exactly";
        }
    }

    //// Compare DE (expansion) arrays
    const double* mDELegacy = legacy->mDE;
    const double* mDEV2 = v2->mDE;
    const double* pDELegacy = legacy->pDE;
    const double* pDEV2 = v2->pDE;
    const double* tDELegacy = legacy->tDE;
    const double* tDEV2 = v2->tDE;
    const double* rhoDELegacy = legacy->rhoDE;
    const double* rhoDEV2 = v2->rhoDE;
    const double* xDELegacy = legacy->xDE;
    const double* xDEV2 = v2->xDE;
    const double* rDELegacy = legacy->rDE;
    const double* rDEV2 = v2->rDE;
    const double* gDELegacy = legacy->gDE;
    const double* gDEV2 = v2->gDE;
    const double* thetaDELegacy = legacy->thetaDE;
    const double* thetaDEV2 = v2->thetaDE;
    const double* massDELegacy = legacy->massDE;
    const double* massDEV2 = v2->massDE;

    // Compare all DE arrays point by point
    for (int j = 0; j <= jDELast; ++j) {
        EXPECT_DOUBLE_EQ(mDELegacy[j], mDEV2[j])
            << "mDE[" << j << "] should match exactly";
        EXPECT_DOUBLE_EQ(pDELegacy[j], pDEV2[j])
            << "pDE[" << j << "] should match exactly";
        EXPECT_DOUBLE_EQ(tDELegacy[j], tDEV2[j])
            << "tDE[" << j << "] should match exactly";
        EXPECT_DOUBLE_EQ(rhoDELegacy[j], rhoDEV2[j])
            << "rhoDE[" << j << "] should match exactly";
        EXPECT_DOUBLE_EQ(xDELegacy[j], xDEV2[j])
            << "xDE[" << j << "] should match exactly";
        EXPECT_DOUBLE_EQ(rDELegacy[j], rDEV2[j])
            << "rDE[" << j << "] should match exactly";
        EXPECT_DOUBLE_EQ(gDELegacy[j], gDEV2[j])
            << "gDE[" << j << "] should match exactly";
        EXPECT_DOUBLE_EQ(thetaDELegacy[j], thetaDEV2[j])
            << "thetaDE[" << j << "] should match exactly";
        EXPECT_DOUBLE_EQ(massDELegacy[j], massDEV2[j])
            << "massDE[" << j << "] should match exactly";
    }

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
