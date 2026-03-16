#include <gtest/gtest.h>
#include <cmath>
#include "src/MOC_2D/solver/NozzleSolver.h"

using namespace moc_2d;
using namespace mp_units;

class TestSolver : public NozzleSolver {
public:
    using NozzleSolver::NozzleSolver;
    std::expected<NozzleResult, std::string> solve() override {
        return std::unexpected("not implemented");
    }
};

TEST(MOC2DTest, ThroatInitialLine) {
    InitialConditions conditions {
        .type = ConditionType::Throat,
        .pressure = 20 * si::mega<si::pascal>,
        .temperature = 2550.2 * si::kelvin,
        .molecular_weight = 21.623e-3 * si::kilogram / si::mole,
        .specific_heat_ratio = 1.2411 * one,
        .ambient_pressure = 1 * si::mega<si::pascal>,
        .velocity = 1103.2 * si::metre / si::second,
        .ideal_specific_impulse = 213.59 * si::second,
        .initial_theta_b = 20.0 * angular::degree,
    };

    auto config = SolverConfiguration::create(
        101,
        1.0 * angular::degree,
        20
    ).value();

    ThroatCurve throat_curve {
        .upstream_radius = 1.5 * one,
        .downstream_radius = 0.382 * one,
    };

    DesignTarget target = ExitMachTarget{ 2.6 * one };

    auto solver = TestSolver(
        NozzleProblem()
            .with_conditions(conditions)
            .with_throat(throat_curve)
            .with_config(config)
            .with_geometry(Geometry::Axisymmetric)
            .with_target(target)
    );

    auto result = solver.calc_initial_throat_line();

    ASSERT_TRUE(result.has_value()) << result.error();

    // TODO: Proper debug logging
    // These values look vaguely correct which is promising anyway

    //std::println("Axial Positions:");
    //for (const auto& row : solver.grid.points) {
    //    for (const auto& x : row) std::print("{} ", x.axial_position);
    //    std::println("");
    //}

    //std::println("Radial Positions:");
    //for (const auto& row : solver.grid.points) {
    //    for (const auto& x : row) std::print("{} ", x.radial_position);
    //    std::println("");
    //}
}

int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
