#pragma once

#include <variant>
#include <expected>

#include <mp-units/systems/si.h>
#include <mp-units/systems/angular.h>

namespace moc_2d {

using mp_units::quantity;
using mp_units::one;

enum class ConditionType { Total, Throat };

struct InitialConditions {
    ConditionType                                           type;               // Whether at t (total) or * (throat, sonic)
    quantity<mp_units::si::pascal>                          pressure;           // p_0
    quantity<mp_units::si::kelvin>                          temperature;        // T_0
    quantity<mp_units::si::kilogram / mp_units::si::mole>   molecular_weight;   // Mwt
    quantity<one>                                           specific_heat_ratio;// gamma
    quantity<mp_units::si::pascal>                          ambient_pressure;   // p_a
    quantity<mp_units::si::metre / mp_units::si::second>    velocity;           // v
    quantity<mp_units::si::second>                          ideal_specific_impulse; // Isp
    quantity<mp_units::angular::radian>                     initial_theta_b;    // theta_B
};

struct ThroatCurve {
    quantity<one> upstream_radius;      // Ru / R*
    quantity<one> downstream_radius;    // Ru / R*
};

enum class Geometry { Planar, Axisymmetric };

struct SolverConfiguration {
    size_t n_characteristics;
    quantity<mp_units::angular::radian> delta_angle_limit;
    size_t n_rrc_above_bd;

    static std::expected<SolverConfiguration, std::string> create(
        size_t n_characteristics,
        quantity<mp_units::angular::radian> delta_angle_limit,
        size_t n_rrc_above_bd
    ) {
        if (n_characteristics < 5)
            return std::unexpected("n_characteristics must be >= 5");
        return SolverConfiguration{n_characteristics, delta_angle_limit, n_rrc_above_bd};
    }
};

struct ExitMachTarget       { quantity<one> mach; };
struct ExitPressTarget      { quantity<mp_units::si::pascal> pressure; };
struct AreaRatioTarget      { quantity<one> eps; };
struct NozzleLengthTarget   { quantity<one> length; };  // normalised to R*
struct EndpointTarget       { quantity<one> x; quantity<one> r; };

using DesignTarget = std::variant<ExitMachTarget, ExitPressTarget, AreaRatioTarget, NozzleLengthTarget, EndpointTarget>;

}
