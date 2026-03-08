#pragma once

#include "mp-units/systems/si.h"
#include "mp-units/systems/angular.h"

using namespace mp_units;

namespace moc_2d::math {
    auto mach_angle(QuantityOf<dimensionless> auto mach);
}
