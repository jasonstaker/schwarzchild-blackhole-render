
#include "libs/stb/stb_image.h"
#include "libs/stb/stb_image_write.h"
#include "math/tensor.hpp"
#include "math/vec4.hpp"

#include <cmath>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

namespace constants {
constexpr double SUN_MASS = 1.989e30;
constexpr double G        = 6.6743e-11;
constexpr double C        = 299792458;
constexpr double mass_sm  = 1000;
constexpr double mass     = mass_sm * SUN_MASS;
}   // namespace constants

using namespace constants;

namespace parameters {
double rs  = (2 * G * mass) / (C * C);
double r_0 = 10 * rs;
vec4<double> c_0(0.0, 0.0, 0.0, r_0);
}   // namespace parameters

using namespace parameters;

vec4<double> convert_to_schwarzchild(const vec4<double> cartesian) {
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

tensor<double> get_schwarzchild_metric(const vec4<double> c_pos) {
    tensor<double> schwarzchild_metric(4, 4);
    vec4<double> s_pos = convert_to_schwarzchild(c_pos);

    schwarzchild_metric(0, 0) = -(1 - (rs / r_0));
    schwarzchild_metric(1, 1) = 1 / (1 - (rs / r_0));
    schwarzchild_metric(2, 2) = r_0 * r_0;
    schwarzchild_metric(3, 3) = r_0 * r_0 * sin(s_pos.theta()) * sin(s_pos.theta());

    return schwarzchild_metric;
}

int main() {
    auto sm = get_schwarzchild_metric(c_0);

    int width, height, channels;
    unsigned char* img = stbi_load("assets/background.jpg", &width, &height, &channels, 0);

    size_t img_size      = width * height * channels;
    int gray_channels    = channels == 4 ? 2 : 1;
    size_t gray_img_size = width * height * gray_channels;

    unsigned char* gray_img = static_cast<unsigned char*>(malloc(gray_img_size));

    for (unsigned char *p = img, *pg = gray_img; p != img + img_size; p += channels, pg += gray_channels) {
        *pg = (uint8_t)((*p + *(p + 1) + *(p + 2)) / 3.0);
        if (channels == 4) { *(pg + 1) = *(p + 3); }
    }

    stbi_write_jpg("output.png", width, height, gray_channels, gray_img, 100);
    stbi_image_free(img);
}