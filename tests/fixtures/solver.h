#pragma once
#include <gtest/gtest.h>
#include "src/MOC_2D/solver/NozzleSolver.h"

class TestNozzleSolver : public moc_2d::NozzleSolver {
public:
    using moc_2d::NozzleSolver::NozzleSolver;
    using moc_2d::NozzleSolver::calc_initial_throat_line;
    using moc_2d::NozzleSolver::calc_rrcs_along_arc;
    using moc_2d::NozzleSolver::calc_transsonic_velocity;
    
    std::expected<moc_2d::NozzleResult, std::string> solve() override { return {}; }

    moc_2d::NozzleGrid get_grid() {return grid;}
    moc_2d::NozzleProblem get_problem() {return problem;}
};
