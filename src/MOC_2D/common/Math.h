#pragma once

#include <cmath>

#include "mp-units/systems/si.h"
#include "mp-units/systems/angular.h"

namespace moc_2d {
namespace math {
    constexpr auto R_universal = 8.314 * mp_units::si::joule / (mp_units::si::mole * mp_units::si::kelvin);
    constexpr auto RAO_CLUSTERING_EXPONENT = 1.5;

    mp_units::quantity<mp_units::angular::radian> mach_angle(mp_units::quantity<mp_units::one> mach);

    mp_units::quantity<mp_units::one> right_characteristic_slope(
        mp_units::quantity<mp_units::angular::radian> flow_angle,
        mp_units::quantity<mp_units::angular::radian> mach_angle
    );

    mp_units::quantity<mp_units::si::pascal> isentropic_pressure(
        mp_units::quantity<mp_units::si::pascal> total_pressure,
        mp_units::quantity<mp_units::one> mach,
        mp_units::quantity<mp_units::one> gamma
    );

    mp_units::quantity<mp_units::si::kelvin> isentropic_temperature(
        mp_units::quantity<mp_units::si::kelvin> total_temperature,
        mp_units::quantity<mp_units::one> mach,
        mp_units::quantity<mp_units::one> gamma
    );

    mp_units::quantity<mp_units::si::kilogram / mp_units::cubic(mp_units::si::metre)> isentropic_density(
        mp_units::quantity<mp_units::si::pascal> total_pressure,
        mp_units::quantity<mp_units::si::kelvin> total_temperature,
        mp_units::quantity<mp_units::si::kilogram / mp_units::si::mole> molecular_weight,
        mp_units::quantity<mp_units::one> mach,
        mp_units::quantity<mp_units::one> gamma
    );
}
}
