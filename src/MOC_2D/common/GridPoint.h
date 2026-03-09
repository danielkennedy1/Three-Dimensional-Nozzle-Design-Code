#pragma once

#include <mp-units/systems/si.h>
#include <mp-units/systems/angular.h>
#include <mp-units/systems/imperial.h>

using namespace mp_units;

namespace moc_2d {
struct GridPoint {
    quantity<one>                                axial_position{};      // x
    quantity<one>                                radial_position{};     // r
    quantity<one>                                mach{};                // M
    quantity<angular::radian>                    flow_angle{};          // theta
    quantity<one>                                specific_heat_ratio{}; // gamma
    quantity<si::pascal>                         pressure{};            // p
    quantity<si::kelvin>                         temperature{};         // T
    quantity<si::kilogram / cubic(si::metre)>    density{};             // rho
};
}
