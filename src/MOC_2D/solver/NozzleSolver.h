#pragma once

#include <expected>
#include <cmath>
#include <iostream>

#include "src/MOC_2D/common/Grid.h"
#include "src/MOC_2D/common/Result.h"
#include "src/MOC_2D/common/Math.h"
#include "src/MOC_2D/config/Problem.h"

// abstract class with things common between nozzles, then different types
// inherit from this and have their own solve() function or something
//
// HAS: a geometry, which is the solution i guess
//
class MOC2DTest_ThroatInitialLine_Test;

namespace moc_2d {
class NozzleSolver {
friend class ::MOC2DTest_ThroatInitialLine_Test;
public:
    NozzleSolver(
        NozzleProblem problem
    ): 
        problem(problem)
    {};
    virtual std::expected<NozzleResult, std::string> solve() = 0;
private:
    NozzleProblem problem;
    NozzleGrid grid{};
protected:
    std::expected<void, std::string> calc_initial_throat_line();
    std::expected<void, std::string> calc_rrcs_along_arc();
};
}
