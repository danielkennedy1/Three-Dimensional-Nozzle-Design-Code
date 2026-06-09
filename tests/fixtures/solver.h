#pragma once

#include <gtest/gtest.h>

#include <mp-units/systems/si.h>
#include <mp-units/systems/angular.h>

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

class NozzleSolverFixture : virtual public ::testing::Test {
protected:
    TestNozzleSolver make_default_throat_solver() {
        moc_2d::InitialConditions conditions {
            .type = moc_2d::ConditionType::Throat,
            .pressure = 11 * mp_units::si::mega<mp_units::si::pascal>,
            .temperature = 2270.7 * mp_units::si::kelvin,
            .molecular_weight = 21.65 * mp_units::si::gram / mp_units::si::mole,
            .specific_heat_ratio = 1.2545 * mp_units::one,
            .ambient_pressure = 1.0 * mp_units::si::mega<mp_units::si::pascal>,
            .velocity = 1045.8 * mp_units::si::metre / mp_units::si::second,
            .ideal_specific_impulse = 213.5 * mp_units::si::second,
            .initial_theta_b = 15 * mp_units::angular::degree,
        };
        auto config = moc_2d::SolverConfiguration::create(101, 1.0 * mp_units::angular::degree, 20).value();
        moc_2d::ThroatCurve throat_curve {
            .upstream_radius = 1.5 * mp_units::one,
            .downstream_radius = 0.382 * mp_units::one,
        };
        return TestNozzleSolver(
            moc_2d::NozzleProblem()
                .with_conditions(conditions)
                .with_throat(throat_curve)
                .with_config(config)
                .with_geometry(moc_2d::Geometry::Axisymmetric)
                .with_target(moc_2d::ExitMachTarget{ 2.6 * mp_units::one })
        );
    }

    TestNozzleSolver make_default_total_solver() {
        moc_2d::InitialConditions conditions {
            .type = moc_2d::ConditionType::Total,
            .pressure = 20 * mp_units::si::mega<mp_units::si::pascal>,
            .temperature = 2550 * mp_units::si::kelvin,
            .molecular_weight = 21.65 * mp_units::si::kilogram / mp_units::si::mole,
            .specific_heat_ratio = 1.2545 * mp_units::one,
            .ambient_pressure = 1.0 * mp_units::si::mega<mp_units::si::pascal>,
            .ideal_specific_impulse = 213.5 * mp_units::si::second,
            .initial_theta_b = 15 * mp_units::angular::degree,
        };
        auto config = moc_2d::SolverConfiguration::create(141, 1.0 * mp_units::angular::degree, 5).value();
        moc_2d::ThroatCurve throat_curve {
            .upstream_radius = 1.5 * mp_units::one,
            .downstream_radius = 0.382 * mp_units::one,
        };
        return TestNozzleSolver(
            moc_2d::NozzleProblem()
                .with_conditions(conditions)
                .with_throat(throat_curve)
                .with_config(config)
                .with_geometry(moc_2d::Geometry::Axisymmetric)
                .with_target(moc_2d::ExitMachTarget{ 2.6 * mp_units::one })
        );
    }
};
