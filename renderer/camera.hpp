#pragma once

#include "../math/core/vec4.hpp"
#include "../math/relativity/geodesic.hpp"
#include "../math/relativity/relativity.hpp"

#include <array>
#include <iostream>
#include <omp.h>
#include <vector>

struct color {
    int r, g, b;
};

class camera {
  private:
    double pixel_samples_scale;   // Color scale factor for a sum of pixel samples
    vec4<double> center;          // Camera center
    vec4<double> pixel00_loc;     // Location of pixel 0, 0
    vec4<double> pixel_delta_u;   // Offset to pixel to the right
    vec4<double> pixel_delta_v;   // Offset to pixel below
    vec4<double> u, v, w;         // Camera frame basis vectors

    void initialize() {
        image_height = int(image_width / aspect_ratio);
        image_height = (image_height < 1) ? 1 : image_height;

        pixel_samples_scale = 1.0 / samples_per_pixel;
        center              = lookfrom;

        auto theta           = degrees_to_radians(vfov);
        auto h               = std::tan(theta / 2);
        auto viewport_height = 2 * h;
        auto viewport_width  = viewport_height * (double(image_width) / image_height);

        w = unit_vector(lookfrom - lookat);
        u = unit_vector(cross(vup, w));
        v = cross(w, u);

        auto viewport_u = viewport_width * u;
        auto viewport_v = viewport_height * -v;

        pixel_delta_u = viewport_u / image_width;
        pixel_delta_v = viewport_v / image_height;

        auto viewport_upper_left = center - w - viewport_u / 2 - viewport_v / 2;
        pixel00_loc              = viewport_upper_left + 0.5 * (pixel_delta_u + pixel_delta_v);
    }

    const geodesic get_geodesic(int i, int j) const {
        geodesic g = {vec4<double>(), vec4<double>()};
        return g;
    }

    vec4<float> sample_square() const {
        return vec4(0, random_double() - 0.5, random_double - 0.5, 0);
    }

    color geodesic_pixel(const geodesic& g) const {
        color c = {0, 0, 0};
        return c;
    }

    // unsigned char* map_to_background(const vec4<float>& direction) {}

  public:
    double aspect_ratio = 1.0;    // Ratio of image width over height
    int image_width     = 100;    // Rendered image width in pixel count
    int image_height;             // Rendered image height
    int samples_per_pixel = 10;   // Count of random samples for each pixel
    int max_depth         = 10;   // Maximum number of ray bounces into scene

    double vfov           = 90;                         // Vertical view angle (field of view)
    vec4<double> lookfrom = c_0;                        // Point camera is looking from
    vec4<double> lookat   = vec4(0.0, 0.0, 0.0, 0.0);   // Point camera is looking at
    vec4<double> vup      = vec4(0.0, 0.0, 0.0, 1.0);   // Camera-relative "up" direction

    unsigned char* render(int width, int height, int channels, unsigned char* img) {
        initialize();

        size_t img_size = (image_width * image_height * channels);
        std::vector<unsigned char> image_buffer(img_size);

        unsigned char* new_img = static_cast<unsigned char*>(malloc(img_size));

        unsigned char* p = new_img;

#pragma omp parallel for schedule(dynamic)
        for (int j = 0; j < image_height; j++) {
            for (int i = 0; i < image_width; i++) {
                color pixel_color = {0, 0, 0};

                for (int sample = 0; sample < samples_per_pixel; sample++) {
                    geodesic g  = get_geodesic(i, j);
                    color pixel = geodesic_pixel(g);

                    pixel_color.r += pixel.r;
                    pixel_color.g += pixel.g;
                    pixel_color.b += pixel.b;
                }

                image_buffer[j * image_height + 3 * image_width + 0] = pixel_color.r * pixel_samples_scale;
                image_buffer[j * image_height + 3 * image_width + 1] = pixel_color.g * pixel_samples_scale;
                image_buffer[j * image_height + 3 * image_width + 2] = pixel_color.b * pixel_samples_scale;
            }
        }

        for (int i = 0; i < img_size; i += 3) {
            p[0] = image_buffer[i];
            p[1] = image_buffer[i + 1];
            p[2] = image_buffer[i + 2];

            p += 3;
        }

        return new_img;
    }
};