#pragma once

#include "../core/vec4.hpp"

struct geodesic {
    vec4<double> position;
    vec4<double> velocity;
    double lambda;
};