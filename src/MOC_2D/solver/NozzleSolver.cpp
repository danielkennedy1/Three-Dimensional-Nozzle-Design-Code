#include "NozzleSolver.h"

namespace moc_2d {

// MAYBE: Put the throat / total polymorphism in their own classes
// (although I think they're the only two types, and throat should really be "sonic")
std::expected<void, std::string> NozzleSolver::calc_initial_throat_line() {
    spdlog::info("Calculting Initial Throat Line");

    int n = problem.solver_config.n_characteristics;
    grid.add_rrc(n);

    grid[0][0].axial_position = 0.0;

    if (problem.conditions.type == ConditionType::Throat) {
        spdlog::info("Calculating from throat conditions");
        // NOTE:  J = kg * m^2 / s^2, so J / kg = m^2 / s^2. sqrt(J/kg) isn't implemented or something so out of mp and back in is neater
        auto speed_of_sound = sqrt((problem.conditions.specific_heat_ratio * moc_2d::math::R_universal /
                                    problem.conditions.molecular_weight * problem.conditions.temperature)
                                       .numerical_value_in(mp_units::si::joule / mp_units::si::kilogram)) *
                              mp_units::si::metre / mp_units::si::second;

        spdlog::debug("Speed of sound = {} m/s", speed_of_sound.numerical_value_in(mp_units::si::metre / mp_units::si::second));

        auto throat_mach = problem.conditions.velocity / speed_of_sound;

        if (throat_mach < 1.0) {
            return std::unexpected("Calculated throat mach number < 1.0");
        }

        for (int i = 0; i < n; i++) {
            // NOTE: Conditions are constant, only position varies
            grid[0][i].mach = throat_mach;
            grid[0][i].specific_heat_ratio = problem.conditions.specific_heat_ratio;
            grid[0][i].flow_angle = 0.0 * mp_units::angular::radian;
            grid[0][i].pressure = problem.conditions.pressure;
            grid[0][i].temperature = problem.conditions.temperature;
            grid[0][i].density = (problem.conditions.pressure * problem.conditions.molecular_weight) / (moc_2d::math::R_universal * problem.conditions.temperature);

            grid[0][i].radial_position = pow(sin(((std::numbers::pi / 2) * (n - i) / n)), moc_2d::math::RAO_CLUSTERING_EXPONENT);
            if (i > 0) {
                auto slope = moc_2d::math::right_characteristic_slope(grid[0][i - 1].flow_angle, moc_2d::math::mach_angle(grid[0][i - 1].mach));
                grid[0][i].axial_position = grid[0][i - 1].axial_position + (grid[0][i].radial_position - grid[0][i - 1].radial_position) / slope;
            }
        }
    }
    else if (problem.conditions.type == ConditionType::Total) {
        return std::unexpected("Total Conditions inital throat line calculation not implemented");
    }
    else {
        return std::unexpected("Unknown Conditions Type inital throat line calculation not implemented");
    }

    return {};
}

    std::expected<void, std::string> NozzleSolver::calc_rrcs_along_arc() { return {}; }
}
