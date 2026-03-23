#pragma once
#include <gtest/gtest.h>
#include "src/MOC_Grid_BDE/MOC_GridCalc_BDE.h"


//blatant claude
class Legacy2d : public ::testing::Test {
protected:
    legacy::MOC_GridCalc calc;

    // Member accessors
    int get_nC() { return calc.nC; }
    int get_maxLRC() { return calc.maxLRC; }
    int get_maxRRC() { return calc.maxRRC; }
    int* get_iLast() { return calc.iLast; }
    int get_nozzleType() { return calc.nozzleType; }
    int get_nozzleGeom() { return calc.nozzleGeom; }
    int get_designParam() { return calc.designParam; }
    int get_lastRRC() { return calc.lastRRC; }
    int get_nRRCAboveBD() { return calc.nRRCAboveBD; }
    int get_nSLi() { return calc.nSLi; }
    int get_nSLj() { return calc.nSLj; }
    int get_throatFlag() { return calc.throatFlag; }
    int get_printMode() { return calc.printMode; }
    int get_jDELast() { return calc.jDELast; }
    int get_iBD() { return calc.iBD; }
    int get_jBD() { return calc.jBD; }

    double** get_mach() { return calc.mach; }
    double** get_pres() { return calc.pres; }
    double** get_temp() { return calc.temp; }
    double** get_rho() { return calc.rho; }
    double** get_x() { return calc.x; }
    double** get_r() { return calc.r; }
    double** get_gamma() { return calc.gamma; }
    double** get_theta() { return calc.theta; }
    double** get_massflow() { return calc.massflow; }

    double* get_mDE() { return calc.mDE; }
    double* get_pDE() { return calc.pDE; }
    double* get_tDE() { return calc.tDE; }
    double* get_rhoDE() { return calc.rhoDE; }
    double* get_xDE() { return calc.xDE; }
    double* get_rDE() { return calc.rDE; }
    double* get_gDE() { return calc.gDE; }
    double* get_thetaDE() { return calc.thetaDE; }
    double* get_massDE() { return calc.massDE; }

    double** get_thrust() { return calc.thrust; }
    double** get_Sthrust() { return calc.Sthrust; }

    double get_pTotal() { return calc.pTotal; }
    double get_tTotal() { return calc.tTotal; }
    double get_molWt() { return calc.molWt; }
    double get_gamma_i() { return calc.gamma_i; }
    double get_pAmbient() { return calc.pAmbient; }
    double get_mThroat() { return calc.mThroat; }
    double get_ispIdeal() { return calc.ispIdeal; }
    double get_mdotErrRatio() { return calc.mdotErrRatio; }
    double* get_designValue() { return calc.designValue; }
    double get_RWTD() { return calc.RWTD; }
    double get_RWTU() { return calc.RWTU; }
    double get_DTLIMIT() { return calc.DTLIMIT; }
    double get_thetaBi() { return calc.thetaBi; }
    double get_thetaBMin() { return calc.thetaBMin; }
    double get_thetaBMax() { return calc.thetaBMax; }
    double get_conCrit() { return calc.conCrit; }
    double get_thetaBAns() { return calc.thetaBAns; }
    double get_initialLineAngle() { return calc.initialLineAngle; }

    // Private method passthroughs
    void InitializeDataMembers() { calc.InitializeDataMembers(); }
    void DeleteDataMembers() { calc.DeleteDataMembers(); }

    int CalcMOC_Grid(double g, double pAmb, int tFlag, double mThroat, int nType,
        int dparam, int geom, int nSLi, int nSLj) {
        return calc.CalcMOC_Grid(g, pAmb, tFlag, mThroat, nType, dparam, geom, nSLi, nSLj);
    }
    int CalcContouredNozzle(int paramType, double* designParamValue,
        double g, double pAmb, int nozzleGeom, int nRRCPlus, int nozzleType,
        int nSLi, int nSLj) {
        return calc.CalcContouredNozzle(paramType, designParamValue, g, pAmb, nozzleGeom, nRRCPlus, nozzleType, nSLi, nSLj);
    }
    int CalcConeNozzle(int paramType, double* designParamValue,
        double pAmb, int nozzleGeom, int nSLi, int nSLj) {
        return calc.CalcConeNozzle(paramType, designParamValue, pAmb, nozzleGeom, nSLi, nSLj);
    }

    dummyStruct CalcHallLine(double rUp, double x, double r, double g, int type) {
        return calc.CalcHallLine(rUp, x, r, g, type);
    }
    int CalcInitialThroatLine(double rUp, int n, double g, double pAmb,
        int geom, int tFlag, double Mach_Throat) {
        return calc.CalcInitialThroatLine(rUp, n, g, pAmb, geom, tFlag, Mach_Throat);
    }
    void Sauer(int i, int geom, double RS) { calc.Sauer(i, geom, RS); }
    int KLThroat(int i, int geom, double RS) { return calc.KLThroat(i, geom, RS); }

    int CalcRRCsAlongArc(int j, double rad, double alphaMax, double dALimit,
        double pAmb, int geom) {
        return calc.CalcRRCsAlongArc(j, rad, alphaMax, dALimit, pAmb, geom);
    }
    int CalcArcWallPoint(int j, double rad, double alphaMax, double dALimit, int geom) {
        return calc.CalcArcWallPoint(j, rad, alphaMax, dALimit, geom);
    }
    void CalcContourWallPoint(int j, int iBottom) { calc.CalcContourWallPoint(j, iBottom); }
    void CalcSpecialWallPoint(int j, double RWTD, double alpha) { calc.CalcSpecialWallPoint(j, RWTD, alpha); }
    void CalcAxialMeshPoint(int j, int i) { calc.CalcAxialMeshPoint(j, i); }
    int CalcInteriorMeshPoints(int j, int iStart, int iEnd, int flag, int geom) {
        return calc.CalcInteriorMeshPoints(j, iStart, iEnd, flag, geom);
    }
    int CheckRRCForNegativePoints(int j) { return calc.CheckRRCForNegativePoints(j); }
    dummyStruct CalcRRCsAlongCone(int j, double dr) { return calc.CalcRRCsAlongCone(j, dr); }
    void CalcConeWallPoint(int j, double wall_angle, int geom) { calc.CalcConeWallPoint(j, wall_angle, geom); }

    dummyStruct CalcLRCDE(int j, int nLRC, double pAmb, int geom, int nRRCPlus,
        int nozzleType, double paramMatch, int pointFlag) {
        return calc.CalcLRCDE(j, nLRC, pAmb, geom, nRRCPlus, nozzleType, paramMatch, pointFlag);
    }
    double CalcMdotBD(int j, double xD) { return calc.CalcMdotBD(j, xD); }
    dummyStruct FindPointE(int j, double xD, double mdot, int nozzleGeom, int nType,
        int nRRCPlus, int pointFlag) {
        return calc.FindPointE(j, xD, mdot, nozzleGeom, nType, nRRCPlus, pointFlag);
    }
    dummyStruct RungeKutta(double dr, double r0, double x0, double mach0,
        double theta0, double gamma0) {
        return calc.RungeKutta(dr, r0, x0, mach0, theta0, gamma0);
    }
    dummyStruct RungeKuttaFehlberg(double h, double r0, double x0, double mach0,
        double theta0, double gamma0) {
        return calc.RungeKuttaFehlberg(h, r0, x0, mach0, theta0, gamma0);
    }
    double Deriv(int i, double r0, double mach0, double theta0, double gamma0) {
        return calc.Deriv(i, r0, mach0, theta0, gamma0);
    }
    int CalcDE(int iD, int jD, int jEnd, int nType, int geom) {
        return calc.CalcDE(iD, jD, jEnd, nType, geom);
    }

    void CalcRemainingMesh(int iD, int jD, int jEnd, int geom) { calc.CalcRemainingMesh(iD, jD, jEnd, geom); }
    void CalcBDERegion(int ii, int jStart, int jEnd, int geom) { calc.CalcBDERegion(ii, jStart, jEnd, geom); }

    void CropNozzleToLength(int jEnd) { calc.CropNozzleToLength(jEnd); }
    void CalcWallContour(int iD, int jD, int nRRC, int geom) { calc.CalcWallContour(iD, jD, nRRC, geom); }

    void CalcIsentropicP_T_RHO(int i, int j, double gamma, double mach) {
        calc.CalcIsentropicP_T_RHO(i, j, gamma, mach);
    }
    dummyStruct CalcIsentropicP_T_RHO(double gamma, double mach) {
        return calc.CalcIsentropicP_T_RHO(gamma, mach);
    }

    double InterpolateMassFlowAlongJUsingX_2D(int j, double xP, int geom) {
        return calc.InterpolateMassFlowAlongJUsingX_2D(j, xP, geom);
    }
    void CalcMassFlowAndThrustAlongMesh(int jStart, int jEnd, double pAmb, int geom) {
        calc.CalcMassFlowAndThrustAlongMesh(jStart, jEnd, pAmb, geom);
    }

    void ResetGrid(int iStart, int iEnd, int jStart, int jEnd) {
        calc.ResetGrid(iStart, iEnd, jStart, jEnd);
    }

    void OutputInitialKernel(int jEnd) { calc.OutputInitialKernel(jEnd); }
    void OutputFinalKernel(int iD, int jD, int jEnd) { calc.OutputFinalKernel(iD, jD, jEnd); }
    void OutputUncroppedKernel(int jEnd) { calc.OutputUncroppedKernel(jEnd); }
    void OutputCenterlineData(std::string fileName) { calc.OutputCenterlineData(fileName); }
    void OutputMOC_Grid() { calc.OutputMOC_Grid(); }
    void OutputSummaryFile() { calc.OutputSummaryFile(); }
    void OutputJ(int j, std::string fileName) { calc.OutputJ(j, fileName); }
    void OutputPrimaryChars(int nType) { calc.OutputPrimaryChars(nType); }
    void OutputStreamlines(int iEnd, int jEnd, int nSLi, int nSLj, int geom) {
        calc.OutputStreamlines(iEnd, jEnd, nSLi, nSLj, geom);
    }
    void OutputTDKRAODataFile(int jStart, int jEnd) { calc.OutputTDKRAODataFile(jStart, jEnd); }

    double CalcMu(double mach) { return calc.CalcMu(mach); }
    double CalcA(double mach, double gamma) { return calc.CalcA(mach, gamma); }
    double CalcB(double mach, double theta, double r) { return calc.CalcB(mach, theta, r); }
    double Calcb(double mach, double theta, double r) { return calc.Calcb(mach, theta, r); }
    double lDyDx(double theta, double mu) { return calc.lDyDx(theta, mu); }
    double rDyDx(double theta, double mu) { return calc.rDyDx(theta, mu); }
    double CalcR(double mach, double theta, double r) { return calc.CalcR(mach, theta, r); }
    double CalcRStar(double mach, double theta, double r) { return calc.CalcRStar(mach, theta, r); }
    double TanAvg(double x, double y) { return calc.TanAvg(x, y); }
    double MM(double mach) { return calc.MM(mach); }
    double CalcPMFunction(double mach, double gamma) { return calc.CalcPMFunction(mach, gamma); }
    double SetThetaB(int paramType, double err, double thetaB0) { return calc.SetThetaB(paramType, err, thetaB0); }

    double FindSAFx(double xNew, int geom) { return calc.FindSAFx(xNew, geom); }
    double CalcNozzleSurfaceArea(int lastRRC, int geom) { return calc.CalcNozzleSurfaceArea(lastRRC, geom); }
};
