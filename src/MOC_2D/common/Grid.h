#pragma once

#include <mp-units/systems/si.h>
#include <mp-units/systems/angular.h>

namespace moc_2d {

using mp_units::quantity;
using mp_units::one;

// NOTE: Some attributes could be made optional here, as all information about a point
// is not determined simultaneously
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

struct NozzleGrid {
    std::vector<std::vector<GridPoint>> points;

    static constexpr size_t max_rrc = 1000;

    void add_rrc(size_t radial_count) {
        if (points.size() > max_rrc) {
            throw std::runtime_error("Right running characteristic count exceeded maximum");
        }
        points.emplace_back(radial_count);
    }

    std::vector<GridPoint> wall_contour() const {
        std::vector<GridPoint> wall;
        wall.reserve(points.size());
        for (const auto& rrc : points)
            wall.push_back(rrc[0]);
        return wall;
    }

    std::span<const GridPoint> rrc(size_t j) const {
        return points[j];
    }

    const GridPoint& centreline_at(size_t j) const {
        return points[j].back();
    }

    std::vector<GridPoint>& operator[](size_t j) {
        return points[j];
    }

    const std::vector<GridPoint>& operator[](size_t j) const {
        return points[j];
    }

};
}
