#pragma once

#include <mp-units/systems/si.h>
#include <mp-units/systems/angular.h>

namespace moc_2d {
// NOTE: Not sure about the units for this, the calcs come out with 
// dimensionless but might be missing quantity information
struct TransonicVelocity {
    mp_units::quantity<mp_units::one> axial;  // U
    mp_units::quantity<mp_units::one> radial; // V
};
}
