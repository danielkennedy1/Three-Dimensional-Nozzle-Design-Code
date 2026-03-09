#pragma once

#include <variant>
#include <expected>

#include <mp-units/systems/si.h>
#include <mp-units/systems/angular.h>

using namespace mp_units;

namespace moc_2d {

struct Total {};
struct Throat {};

template<typename T>
concept Conditions = std::same_as<T, Total> || std::same_as<T, Throat>;

template<Conditions c>
struct InitialConditions {
    quantity<si::pascal>                          pressure;           // p₀
    quantity<si::kelvin>                          temperature;        // T₀
    quantity<si::kilogram / si::mole>             molecular_weight;   // Mwt
    quantity<one>                                 gamma;              // γ
    quantity<si::pascal>                          ambient_pressure;   // pₐ
    quantity<si::metre / si::second>              velocity;           // v
    quantity<si::second>                          ideal_specific_impulse; // Isp
    quantity<angular::radian>                     initial_theta_b;    // θ_B
};

using InitialTotalConditions  = InitialConditions<Total>;
using InitialThroatConditions = InitialConditions<Throat>;

struct ThroatCurve {
    quantity<one> upstream_radius;      // Ru / R*
    quantity<one> downstream_radius;    // Ru / R*
};

enum class Geometry { Planar, Axisymmetric };

struct SolverConfiguration {
    size_t n_characteristics;
    quantity<angular::radian> delta_angle_limit;
    size_t n_rrc_above_bd;

    static std::expected<SolverConfiguration, std::string> create(
        size_t n_characteristics,
        quantity<angular::radian> delta_angle_limit,
        size_t n_rrc_above_bd
    ) {
        if (n_characteristics < 5)
            return std::unexpected("n_characteristics must be >= 5");
        return SolverConfiguration{n_characteristics, delta_angle_limit, n_rrc_above_bd};
    }
};

struct ExitMachTarget       { quantity<one> mach; };
struct ExitPressTarget      { quantity<si::pascal> pressure; };
struct AreaRatioTarget      { quantity<one> eps; };
struct NozzleLengthTarget   { quantity<one> length; };  // normalised to R*
struct EndpointTarget       { quantity<one> x; quantity<one> r; };

using DesignTarget = std::variant<ExitMachTarget, ExitPressTarget, AreaRatioTarget, NozzleLengthTarget, EndpointTarget>;

}
