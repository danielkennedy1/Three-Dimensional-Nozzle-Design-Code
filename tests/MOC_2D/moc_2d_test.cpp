#include <gtest/gtest.h>
#include <cmath>

#include "spdlog/spdlog.h"

#include "tests/fixtures/solver.h"

using namespace moc_2d;
using namespace mp_units;

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

    auto solver = TestNozzleSolver(
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

    spdlog::debug("Axial Positions:");
    for (const auto& row : solver.get_grid().points) {
        for (const auto& x : row) spdlog::debug("{} ", x.axial_position.numerical_value_in(mp_units::one));
        spdlog::debug("");
    }

    spdlog::debug("Radial Positions:");
    for (const auto& row : solver.get_grid().points) {
        for (const auto& x : row) spdlog::debug("{} ", x.radial_position.numerical_value_in(mp_units::one));
        spdlog::debug("");
    }
}

int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    spdlog::set_level(spdlog::level::debug);
    return RUN_ALL_TESTS();
}
