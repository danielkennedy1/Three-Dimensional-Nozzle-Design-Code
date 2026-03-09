#pragma once

#include "Configuration.h"

namespace moc_2d {
struct NozzleProblem {
    InitialConditions<Total>  conditions;
    ThroatCurve               throat_curve;
    Geometry                  geometry;
    SolverConfiguration       solver_config;
    DesignTarget              target;

    NozzleProblem& with_conditions(InitialConditions<Total> c) {
        conditions = c; return *this;
    }
    NozzleProblem& with_throat(ThroatCurve t) {
        throat_curve = t; return *this;
    }
    NozzleProblem& with_geometry(Geometry g) {
        geometry = g; return *this;
    }
    NozzleProblem& with_config(SolverConfiguration c) {
        solver_config = c; return *this;
    }
    NozzleProblem& with_target(DesignTarget t) {
        target = t; return *this;
    }
};

}
