#include "Math.h"

using namespace mp_units;

namespace moc_2d::math {
  auto mach_angle(QuantityOf<dimensionless> auto mach)
  {
      return angular::asin(1 / mach);
  }
}
