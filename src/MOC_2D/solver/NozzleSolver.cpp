#include "NozzleSolver.h"

#include "spdlog/spdlog.h"

#include "src/MOC_2D/common/Math.h"

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

        for (int i = 0; i < n; i++) {

            grid[0][i].specific_heat_ratio = problem.conditions.specific_heat_ratio;

            auto transonic_velocity = calc_transsonic_velocity(
                grid[0][i], 
                problem.geometry, 
                problem.throat_curve.upstream_radius
            );

            if (fabs(transonic_velocity.radial.numerical_value_in(mp_units::one)) < 1e-5) transonic_velocity.radial = 0;

            auto flow_angle = mp_units::angular::atan2(transonic_velocity.radial, transonic_velocity.axial);
            if (flow_angle.numerical_value_in(mp_units::angular::radian) < 1e-5) flow_angle = 0 * mp_units::angular::radian;

            grid[0][i].flow_angle = flow_angle;

            auto mach = sqrt(
                ( 
                    transonic_velocity.axial * transonic_velocity.axial 
                    + 
                    transonic_velocity.radial * transonic_velocity.radial 
                ).numerical_value_in(mp_units::one)
            );

            if (mach < 1.0) {
                return std::unexpected("Subsonic flow in initial throat line, increase initial line angle ( ?? )");
            }

            grid[0][i].mach = mach;

            grid[0][i].pressure = moc_2d::math::isentropic_pressure(
                problem.conditions.pressure,
                mach,
                problem.conditions.specific_heat_ratio
            );
            grid[0][i].temperature = moc_2d::math::isentropic_temperature(
                problem.conditions.temperature,
                mach, 
                problem.conditions.specific_heat_ratio
            );
            grid[0][i].density = moc_2d::math::isentropic_density(
                problem.conditions.pressure,
                problem.conditions.temperature,
                problem.conditions.molecular_weight,
                mach,
                problem.conditions.specific_heat_ratio
            );

        }
    }
    else {
        return std::unexpected("Unknown Conditions Type inital throat line calculation not implemented");
    }

    return {};
}

// NOTE: IDK if I'll ever be happy with this, many magic numbers which are just plucked out of
// J. R. Kliegel and J. N. Levine, “Transonic flow in small throat radius of curvature nozzles.,” AIAA Journal, vol. 7, no. 7, pp. 1375–1378, Jul. 1969, doi: 10.2514/3.5355
// which references and corrects
// I. M. HALL, "Transonic Flow in Two-Dimensional and Axially-Symmetric Nozzles.," The Quarterly Journal of Mechanics and Applied Mathematics, Volume 15, Issue 4, November 1962, Pages 487–508, https://doi.org/10.1093/qjmam/15.4.487
// It looks like there's another "Sauer" method which probably implements another paper's approach
//
// MAYBE: have another swing at cleaning up esp. wrt. operator ordering because there's a lot of implicit ordering in here

TransonicVelocity NozzleSolver::calc_transsonic_velocity(
    GridPoint point,
    Geometry geometry,
    mp_units::quantity<mp_units::one> upstream_radius
) {
    assert(point.axial_position != 0.0);
    assert(point.radial_position != 0.0);
    assert(point.specific_heat_ratio != 0.0);

    const auto x = point.axial_position;
    const auto r = point.radial_position;
    const auto gamma = point.specific_heat_ratio;

    mp_units::quantity<one> z, u_1, u_2, u_3, v_1, v_2, v_3, U, V;

    if (geometry == Geometry::Axisymmetric) {
		z = x * sqrt(  (  2 * upstream_radius / (  gamma + 1 ) ).numerical_value_in( mp_units::one) ); // Hall Eq. 12
        u_1 = r * r / 2 - 0.25 + z;
        u_2 = ( 2 * gamma + 9 ) * r * r * r * r / 24 - ( 4 * gamma + 15 ) * r * r / 24 + ( 10 * gamma + 57 ) / 288 + z * ( r * r - 5 / 8 ) - ( 2 * gamma - 3 ) * z * z / 6;
        u_3 = ( 556 * gamma * gamma + 1737 * gamma + 3069 ) * r * r * r * r * r * r / 10368 - ( 388 * gamma * gamma + 1161 * gamma + 1881 ) * r * r * r * r / 2304 + 
        ( 304 * gamma * gamma + 831 * gamma + 1242 ) * r * r / 1728 - ( 2708 * gamma * gamma + 7839 * gamma + 14211 ) / 82944 + 
        z * ( ( 52 * gamma * gamma + 51 * gamma + 327 ) * r * r * r * r / 34 - ( 52 * gamma * gamma + 75 * gamma + 279 ) * r * r / 192 + ( 92 * gamma * gamma + 180 * gamma + 639 ) / 1152 ) + 
        z * z * ( - ( 7 * gamma - 3 ) * r * r / 8 + ( 13 * gamma - 27 ) / 48 ) + ( 4 * gamma * gamma - 57 * gamma + 27 ) * z * z * z / 144;

        v_1 = r * r * r / 4 - r / 4 + r * z;
        v_2 = ( gamma + 3 ) * r * r * r * r * r / 9 - ( 20 * gamma + 63 ) * r * r * r / 96 + ( 28 * gamma + 93 ) * r / 288 + z * ( ( 2 * gamma + 9 ) * r * r * r / 6 - ( 4 * gamma + 15 ) * r / 12 ) + r * z * z;
        v_3 = ( 6836 * gamma * gamma + 23031 * gamma + 30627 ) * r * r * r * r * r * r * r / 82944 - ( 3380 * gamma * gamma + 11391 * gamma + 15291 ) * r * r * r * r * r / 13824 + 
        ( 3424 * gamma * gamma + 11271 * gamma + 15228 ) * r * r * r / 13824 - ( 7100 * gamma * gamma + 22311 * gamma + 30249 ) * r / 82944 + 
        z * ( ( 556 * gamma * gamma + 1737 * gamma + 3069 ) * r * r * r * r * r / 1728 * ( 388 * gamma * gamma + 1161 * gamma + 1181 ) * r * r / 576 + 
        ( 304 * gamma * gamma + 831 * gamma + 1242 ) * r / 864 ) + z * z * ( ( 52 * gamma * gamma + 51 * gamma + 327 ) * r * r * r / 192 - ( 52 * gamma * gamma + 75 * gamma + 279 ) * r / 192 ) - 
        z * z * z * ( 7 * gamma - 3 ) * r / 12;


        U = 1 
        + u_1 / ( upstream_radius + 1 ) 
        + ( u_1 + u_2 ) / ( ( upstream_radius + 1 ) * ( upstream_radius + 1 ) ) 
        + ( u_1 + 2 * u_2 + u_3 ) / ( ( upstream_radius + 1 ) * ( upstream_radius + 1 ) * ( upstream_radius + 1 ) );
        V = sqrt( ( ( gamma + 1 ) / ( 2 * ( upstream_radius + 1 ) ) ).numerical_value_in(mp_units::one) ) * ( v_1 / ( upstream_radius + 1 ) 
        + ( 1.5 * v_1 + v_2 ) / ( ( upstream_radius + 1 ) * ( upstream_radius + 1 ) ) 
        + ( 15. / 8. * v_1 + 2.5 * v_2 + v_3 ) / ( ( upstream_radius + 1 ) * ( upstream_radius + 1 ) * ( upstream_radius + 1 ) ) );

        return TransonicVelocity(U, V);
    } else if (geometry == Geometry::Planar) {

		z = x * sqrt( ( upstream_radius/(gamma+1) ).numerical_value_in(mp_units::one));

		u_1 = 0.5 * r * r - 1 / 6 + z;
        u_2 = ( r+6 ) * r * r * r * r / 18 - ( 2 * gamma+9 ) * r * r / 18 + ( gamma+30 ) / 270 + z * ( r * r-0.5 ) - ( 2 * gamma-3 ) * z * z / 6;
		u_3 = ( 362 * gamma * gamma+1449 * gamma+3177 ) * r * r * r * r * r * r / 12960 - ( 194 * gamma * gamma + 837 * gamma + 1665 ) * r * r * r * r / 2592 + 
			( 854 * gamma * gamma + 3687 * gamma + 6759 ) * r * r / 12960 - ( 782 * gamma * gamma + 5523 + 2 * gamma * 2887 ) / 272160 + 
			z * ( ( 26 * gamma * gamma + 27 * gamma + 237 ) * r * r * r * r / 288 - ( 26 * gamma * gamma + 51 * gamma + 189 ) * r * r / 144 +
			( 134 * gamma * gamma + 429 * gamma + 1743 ) / 4320 ) + z * z * ( -5 * gamma * r * r / 4 + ( 7 * gamma - 18 ) / 36 ) + 
			z * z * z * ( 2 * gamma * gamma - 33 * gamma + 9 ) / 72;

		v_1 = r * r * r / 6 - r / 6 + r * z;
		v_2 = ( 22 * gamma+75 ) * r * r * r * r * r / 360 - ( 5 * gamma+21 ) * r * r * r / 54 + ( 34 * gamma+195 ) * r / 1080 + z / 9 * ( ( 2 * gamma+12 ) * r * r * r - ( 2 * gamma+9 ) * r ) + r * z * z;
		v_3 = ( 6574 * gamma * gamma + 26481 * gamma + 40059 ) * r * r * r * r * r * r * r / 181440 - ( 2254 * gamma * gamma + 10113 * gamma + 16479 ) * r * r * r * r * r / 25920 + 
			( 5026 * gamma * gamma + 25551 * gamma + 46377 ) * r * r * r / 77760 - ( 7570 * gamma * gamma + 45927 * gamma + 98757 ) * r / 544320 + 
			z * ( ( 362 * gamma * gamma + 1449 * gamma + 3177 ) * r * r * r * r * r / 2160  *  ( 194 * gamma * gamma + 837 * gamma + 1665 ) * r * r * r / 648 + 
			( 854 * gamma * gamma + 3687 * gamma + 6759 ) * r / 6480 ) + z * z * ( ( 26 * gamma * gamma + 27 * gamma + 237 ) * r * r * r / 144 - ( 26 * gamma * gamma + 51 * gamma + 189 ) / 144 ) + 
			z * z * z * ( -5 * gamma * r / 6 );
		
		U = 1 
            + u_1 / upstream_radius 
            + u_2 / upstream_radius / upstream_radius 
            + u_3 / upstream_radius / upstream_radius / upstream_radius;
		V = sqrt(
            (
                ( 
                    ( gamma + 1 ) / upstream_radius 
                ) 
                * 
                ( 
                      v_1 / upstream_radius 
                    + v_2 / upstream_radius / upstream_radius 
                    + v_3 / upstream_radius / upstream_radius / upstream_radius 
                ) 
            ).numerical_value_in(mp_units::one) 
        );
    } else {
        throw std::logic_error("Transonic velocity for unknown geometry type not implemented");
    }

    return moc_2d::TransonicVelocity{
        .axial = U,
        .radial = V
    };

}

std::expected<void, std::string> NozzleSolver::calc_rrcs_along_arc() { return {}; }
}
