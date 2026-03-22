#include "Math.h"

using namespace mp_units;

// NOTE: For the internals here, I'm bringing the mp_units types back
// into standard double values, to make the calculations clearer

namespace moc_2d::math {
    

    quantity<angular::radian> mach_angle(quantity<one> mach)
    {
        return angular::asin(1.0 / mach);
    }

    // rDyDx
    mp_units::quantity<mp_units::one> right_characteristic_slope(
        mp_units::quantity<mp_units::angular::radian> flow_angle, // theta
        mp_units::quantity<mp_units::angular::radian> mach_angle // mu
    ) {
        return mp_units::angular::tan(flow_angle - mach_angle);
    }

    quantity<si::pascal> isentropic_pressure(
        quantity<si::pascal> total_pressure,
        quantity<one> mach,
        quantity<one> gamma
    ) {
        double p_t = total_pressure.numerical_value_in(si::pascal);
        double g = gamma.numerical_value_in(one);
        double M = mach.numerical_value_in(one);

        double ratio = std::pow(1 + ( g - 1 ) / 2 * M * M, g / (g - 1));
        
        return (p_t / ratio) * si::pascal;
    }

    quantity<si::kelvin> isentropic_temperature(
        quantity<si::kelvin> total_temperature,
        quantity<one> mach,
        quantity<one> gamma
    ) {
        double t_t = total_temperature.numerical_value_in(si::kelvin);
        double g = gamma.numerical_value_in(one);
        double M = mach.numerical_value_in(one);
        
        double ratio = 1 + ( g - 1 ) / 2 * M * M;
        
        return (t_t / ratio) * si::kelvin;
    }

    quantity<si::kilogram / cubic(si::metre)> isentropic_density(
        quantity<si::pascal> total_pressure,
        quantity<si::kelvin> total_temperature,
        quantity<si::kilogram / si::mole> molecular_weight,
        quantity<one> mach,
        quantity<one> gamma
    ) {
        auto pressure = isentropic_pressure(total_pressure, mach, gamma);
        auto temperature = isentropic_temperature(total_temperature, mach, gamma);

        return (pressure * molecular_weight) / (R_universal * temperature);
    }

}
