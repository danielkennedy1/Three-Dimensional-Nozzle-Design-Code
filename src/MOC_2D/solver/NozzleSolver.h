#pragma once

#include <expected>
#include <cmath>
#include <cassert>


#include "src/MOC_2D/config/Problem.h"

#include "src/MOC_2D/common/Grid.h"
#include "src/MOC_2D/common/Result.h"
#include "src/MOC_2D/common/Model.h"

// NOTE: There's probably scope for neater polymorphism
// with the velocity calculations in the transsonic region given 
// total conditions. The transonic velocity calculation "belongs"
// to the geometry type, i.e. axisymmetric / planar, and whether
// it gets called is decided by the conditions type, so possible
// to offload the complication of that dispatch there. could be
// very deep DI chain so potentially worse to comprehend and debug


namespace moc_2d {
class NozzleSolver {
public:
    NozzleSolver(
        NozzleProblem problem
    ): 
        problem(problem)
    {};
    virtual std::expected<NozzleResult, std::string> solve() = 0;
protected:
    NozzleProblem problem;
    NozzleGrid grid{};
    std::expected<void, std::string> calc_initial_throat_line();
    std::expected<void, std::string> calc_rrcs_along_arc();
    TransonicVelocity calc_transsonic_velocity(GridPoint point, Geometry geometry, mp_units::quantity<mp_units::one> upstream_radius);
};
}
