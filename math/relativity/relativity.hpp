#pragma once

#include "../core/tensor.hpp"
#include "../core/vec4.hpp"
#include "geodesic.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>

namespace constants {
constexpr double G  = 1.0;
constexpr double C  = 1.0;
constexpr double PI = M_PI;

constexpr double M = 1.0;

constexpr double R_S      = 2.0;
constexpr double R_0      = 8.0;
constexpr double R_ESCAPE = 100.0;

constexpr double THETA_0 = PI / 2 + 1e-4;
constexpr double PHI_0   = 0;

constexpr double H_MAX = 0.1;
constexpr double H_MIN = 1e-5;
}   // namespace constants

using namespace constants;

namespace parameters {
double x = R_0 * std::sin(THETA_0) * std::cos(PHI_0);
double y = R_0 * std::sin(THETA_0) * std::sin(PHI_0);
double z = R_0 * std::cos(THETA_0);

vec4<double> c_0(0.0, x, y, z, false);
}   // namespace parameters

using namespace parameters;

namespace relativity {

// EFFECTS: produces a tetrad for the initial camera position
inline std::vector<vec4<double>> get_tetrad() {
    vec4<double> e_0(1.0 / std::sqrt(1 - (R_S / R_0)), 0.0, 0.0, 0.0, true);
    vec4<double> e_1(0.0, std::sqrt(1.0 - (R_S / R_0)), 0.0, 0.0, true);
    vec4<double> e_2(0.0, 0.0, (1.0 / R_0), 0.0, true);
    vec4<double> e_3(0.0, 0.0, 0.0, 1.0 / (R_0 * std::sin(THETA_0)), true);

    return {e_0, e_1, e_2, e_3};
}

// REQUIRES: cartesian position is a cartesian coordinate
// EFFECTS: converts given cartesian position coordinates in to schwarzschild position coordinates
inline vec4<double> convert_to_schwarzschild_pos(const vec4<double>& cartesian) {
    if (cartesian.schwarzschild) { return cartesian; }
    double t = cartesian[0];
    double x = cartesian[1];
    double y = cartesian[2];
    double z = cartesian[3];

    double r     = std::sqrt(x * x + y * y + z * z);
    double theta = std::acos(z / r);
    double phi   = std::atan2(y, x);

    vec4<double> v(t, r, theta, phi, true);

    return v;
}

// REQUIRES: schwarzschild position is a schwarzschild coordinate
// EFFECTS: converts given schwarzschild position coordinates in to cartesian position coordinates
inline vec4<double> convert_to_cartesian_pos(const vec4<double>& schwarzschild) {
    if (!schwarzschild.schwarzschild) { return schwarzschild; }
    double t     = schwarzschild[0];
    double r     = schwarzschild[1];
    double theta = schwarzschild[2];
    double phi   = schwarzschild[3];

    if (r == 0) { r = 0.00001; }

    double x = r * std::sin(theta) * std::cos(phi);
    double y = r * std::sin(theta) * std::sin(phi);
    double z = r * std::cos(theta);

    vec4<double> v(t, x, y, z, false);

    return v;
}

// REQUIRES: cartesian is from the camera position
// EFFECTS: converts the cartesian velocity to schwarzschild
inline vec4<double> convert_to_schwarzschild_vel(const vec4<double>& cartesian) {
    if (cartesian.schwarzschild) { return cartesian; }

    auto tetrad = get_tetrad();

    vec4<double> v =
      tetrad[0] * cartesian[0] + tetrad[1] * cartesian[1] + tetrad[2] * cartesian[2] + tetrad[3] * cartesian[3];
    v.schwarzschild = true;

    return v;
}

// EFFECTS: given c_pos, outputs the schwarzschild metric there
inline tensor<double> get_schwarzschild_metric(const vec4<double> pos) {
    tensor<double> schwarzschild_metric(4, 4);
    vec4<double> s_pos = convert_to_schwarzschild_pos(pos);

    double r = s_pos.r();

    schwarzschild_metric(0, 0) = -(1 - (R_S / r));
    schwarzschild_metric(1, 1) = 1 / (1 - (R_S / r));
    schwarzschild_metric(2, 2) = r * r;
    schwarzschild_metric(3, 3) = r * r * sin(s_pos.theta()) * sin(s_pos.theta());

    return schwarzschild_metric;
}

// REQUIRES: u is a schwarzschild position vector
// EFFECTS: returns the christoffel symbol corresponding to the correct mu, alpha, beta, and vector
inline double get_christoffel(int mu, int alpha, int beta, const vec4<double> u) {
    auto r     = u.r();
    auto theta = u.theta();

    switch (mu) {
        case 0:
            if ((alpha == 1 && beta == 0) || (alpha == 0 && beta == 1)) { return M / (r * (r - 2 * M)); }
            break;
        case 1:
            if (alpha == 0 && beta == 0) {
                return (M * (r - 2 * M)) / (r * r * r);
            } else if (alpha == 1 && beta == 1) {
                return -(M / (r * (r - 2 * M)));
            } else if (alpha == 2 && beta == 2) {
                return -(r * (1 - (2 * M) / r));
            } else if (alpha == 3 && beta == 3) {
                return -(r * (1 - (2 * M) / r)) * (std::sin(theta) * std::sin(theta));
            }
            break;
        case 2:
            if ((alpha == 1 && beta == 2) || (alpha == 2 && beta == 1)) {
                return 1 / r;
            } else if (alpha == 3 && beta == 3) {
                return -(std::sin(theta) * std::cos(theta));
            }
            break;
        case 3:
            if ((alpha == 3 && beta == 1) || (alpha == 1 && beta == 3)) {
                return 1 / r;
            } else if ((alpha == 3 && beta == 2) || (alpha == 2 && beta == 3)) {
                double s = std::sin(theta);
                if (std::abs(s) < 1e-6) s = (s >= 0 ? 1e-6 : -1e-6);
                return std::cos(theta) / s;
            }
            break;
    }

    return 0.0;
}

// REQUIRES: g is a light-like geodesic
// EFFECTS: returns the new state geodesic based on g
inline geodesic rk4_step(const geodesic& g) {
    vec4<double> acc;
    auto pos = g.position;
    auto vel = g.velocity;

    for (int mu = 0; mu < 4; mu++) {
        double sum = 0.0;

        for (int alpha = 0; alpha < 4; alpha++) {
            for (int beta = 0; beta < 4; beta++) {
                sum -= get_christoffel(mu, alpha, beta, pos) * vel[alpha] * vel[beta];
            }
        }

        acc[mu] = sum;
    }

    acc.schwarzschild = true;

    geodesic f_g = {vel, acc};
    return f_g;
}

// REQUIRES: g's components are in schwarzschild
// EFFECTS: returns a geodesic with a fixed time velocity component so it is light-like
inline geodesic enforce_null_condition(const geodesic& g) {
    auto v = g.velocity;
    auto p = g.position;

    auto r     = g.position[1];
    auto theta = g.position[2];

    auto u_r  = g.velocity[1];
    auto u_th = g.velocity[2];
    auto u_p  = g.velocity[3];

    auto f = 1 - (2 * M) / g.position.r();

    double u_t =
      std::sqrt((((u_r * u_r) / f) + r * r * u_th * u_th + r * r * std::sin(theta) * std::sin(theta) * u_p * u_p) / f);

    geodesic g_null = {g.position, vec4(u_t, u_r, u_th, u_p, true), g.lambda};

    return g_null;
}

// EFFECTS: integrates the geodesic one full step
inline geodesic integrate(const geodesic& g) {
    double r = g.position.r();
    double h = (H_MAX * r) / (r + R_0);

    if (r <= R_S * 1.5) h *= .25;

    h = std::clamp(h, H_MIN, H_MAX);

    geodesic g_1 = g;
    auto k_1     = rk4_step(g_1);
    geodesic g_2 = {g.position + (h / 2) * k_1.position, g.velocity + (h / 2) * k_1.velocity};
    auto k_2     = rk4_step(g_2);
    geodesic g_3 = {g.position + (h / 2) * k_2.position, g.velocity + (h / 2) * k_2.velocity};
    auto k_3     = rk4_step(g_3);
    geodesic g_4 = {g.position + h * k_3.position, g.velocity + h * k_3.velocity};
    auto k_4     = rk4_step(g_4);

    geodesic g_final = {g.position + (h / 6) * (k_1.position + 2 * k_2.position + 2 * k_3.position + k_4.position),
      g.velocity + (h / 6) * (k_1.velocity + 2 * k_2.velocity + 2 * k_3.velocity + k_4.velocity)};

    return g_final;
}

// REQUIRES: g is a light-like geodesic
// EFFECTS: returns whether the geodesic should stop being integrated
inline bool check_termination(const geodesic& g) {
    auto r = g.position.r();

    return (r > R_ESCAPE) || (r <= R_S);
}
}   // namespace relativity

using namespace relativity;