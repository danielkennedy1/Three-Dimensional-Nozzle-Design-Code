#pragma once

#include <mp-units/systems/si.h>
#include <mp-units/systems/angular.h>

namespace moc_2d {
struct NozzleResult {
    mp_units::quantity<mp_units::si::second>                            specific_impulse;        // Isp
    mp_units::quantity<mp_units::one>                                   thrust_coefficient;      // Cfg
    mp_units::quantity<mp_units::si::metre / mp_units::si::second>      characteristic_velocity; // C*
    mp_units::quantity<mp_units::one>                                   discharge_coefficient;   // CD
    mp_units::quantity<mp_units::si::kilogram / mp_units::si::second>   mass_flow;               // mdot
    mp_units::quantity<mp_units::one>                                   expansion_ratio;         // epsilon
    mp_units::quantity<mp_units::one>                                   nozzle_length;           // L/R*
    mp_units::quantity<mp_units::si::metre * mp_units::si::metre>       surface_area;            // A_s
    mp_units::quantity<mp_units::one>                                   exit_mach_wall;          // M_e,wall
    mp_units::quantity<mp_units::one>                                   exit_mach_centreline;    // M_e,cl
    mp_units::quantity<mp_units::angular::radian>                       theta_b;                 // theta_B (converged)
};
}
