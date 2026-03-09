#pragma once

#include <mp-units/systems/si.h>
#include <mp-units/systems/angular.h>
#include <mp-units/systems/imperial.h>

namespace moc_2d {

using mp_units::quantity;
using mp_units::one;

struct GridPoint {
    quantity<one>                                                   axial_position{};      // x
    quantity<one>                                                   radial_position{};     // r
    quantity<one>                                                   mach{};                // M
    quantity<mp_units::angular::radian>                             flow_angle{};          // theta
    quantity<one>                                                   specific_heat_ratio{}; // gamma
    quantity<mp_units::si::pascal>                                  pressure{};            // p
    quantity<mp_units::si::kelvin>                                  temperature{};         // T
    quantity<mp_units::si::kilogram / cubic(mp_units::si::metre)>   density{};             // rho
};
}
