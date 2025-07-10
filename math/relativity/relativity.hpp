#pragma once

#include "../core/vec4.hpp"
#include "../core/tensor.hpp"

#include <cmath>

namespace constants {
constexpr double SUN_MASS = 1.989e30;
constexpr double G        = 6.6743e-11;
constexpr double C        = 299792458;
constexpr double PI       = M_PI;
constexpr double mass_sm  = 1000;
constexpr double mass     = mass_sm * SUN_MASS;
constexpr double rs       = (2 * G * mass) / (C * C);
constexpr double r_0      = 10 * rs;
constexpr double theta_0  = PI / 4;
constexpr double phi_0    = 10;
}   // namespace constants

using namespace constants;

namespace parameters {
double x = r_0 * std::sin(theta_0) * std::cos(phi_0);
double y = r_0 * std::sin(theta_0) * std::sin(phi_0);
double z = r_0 * std::cos(theta_0);

vec4<double> c_0(0.0, x, y, z);
}   // namespace parameters

using namespace parameters;

// REQUIRES: cartesian is a cartesian coordinate
// EFFECTS: converts given cartesian coordinates in to schwarzchild coordinates
inline vec4<double> convert_to_schwarzchild(const vec4<double> cartesian) {
    vec4<double> schwarzchild;

    auto t = cartesian[0];
    auto x = cartesian[1];
    auto y = cartesian[2];
    auto z = cartesian[3];

    schwarzchild[0] = cartesian[0];
    schwarzchild[1] = sqrt(x * x + y * y + z * z);
    schwarzchild[2] = acos(z / r_0);
    schwarzchild[3] = atan2(y, x);

    return schwarzchild;
}

// REQUIRES: c_pos is a cartesian position
// EFFECTS: given c_pos, outputs the schwarzchild metric there
inline tensor<double> get_schwarzchild_metric(const vec4<double> c_pos) {
    tensor<double> schwarzchild_metric(4, 4);
    vec4<double> s_pos = convert_to_schwarzchild(c_pos);

    schwarzchild_metric(0, 0) = -(1 - (rs / r_0));
    schwarzchild_metric(1, 1) = 1 / (1 - (rs / r_0));
    schwarzchild_metric(2, 2) = r_0 * r_0;
    schwarzchild_metric(3, 3) = r_0 * r_0 * sin(s_pos.theta()) * sin(s_pos.theta());

    return schwarzchild_metric;
}

// REQUIRES: c_pos is a cartesian coordinate
// EFFECTS: produces a tetrad for the cartesian coordinate
inline std::vector<vec4<double>> get_tetrad(const vec4<double> c_pos) {
    auto s_pos   = convert_to_schwarzchild(c_pos);
    double r     = s_pos.r();
    double theta = s_pos.theta();

    auto e_0 = vec4<double>(1 / std::sqrt(1 - rs / r), 0, 0, 0);
    auto e_1 = vec4<double>(0, std::sqrt(1 - rs / r), 0, 0);
    auto e_2 = vec4<double>(0, 0, 1 / r, 0);
    auto e_3 = vec4<double>(0, 0, 0, 1 / (r * std::sin(theta)));

    return {e_0, e_1, e_2, e_3};
}