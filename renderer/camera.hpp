#pragma once

#include "../math/core/vec4.hpp"
#include "../math/relativity/geodesic.hpp"
#include "../math/relativity/relativity.hpp"

#include <algorithm>
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
    unsigned char* background;
    int background_width;
    int background_height;

    // MODIFIES: this 
    // EFFECTS: initializes all the necessary variables for rendering
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

    // REQUIRES: i < image_width && j < image_height
    // EFFECTS: returns a random geodesic 0.5 units away from the pixel center
    geodesic get_geodesic(int i, int j) const {
        auto v_i = pixel_delta_u * i;
        auto v_j = pixel_delta_v * j;

        auto pixel_sample = pixel00_loc + v_i + v_j + (pixel_delta_u * (random_double() - 0.5)) +
                            (pixel_delta_v * (random_double() - 0.5));

        geodesic g = {center, unit_vector(pixel_sample - center) * constants::C, 0.0};
        return g;
    }

    // EFFECTS: returns the pixel mapped from spherical coordinates to texture coordinates
    color geodesic_pixel(const geodesic& g) const {
        auto velocity = unit_vector(g.velocity);

        double theta = std::acos(velocity[3]);
        double phi   = std::atan2(velocity[2], velocity[1]);

        double u = std::fmod((phi + PI)/(2*PI) + 0.5, 1.0);
        double v = (theta / PI);
        int x    = static_cast<int>(u * (background_width - 1));
        int y    = static_cast<int>(v * (background_height - 1));

        int idx   = (y * background_width + x) * 3;
        int red   = background[idx];
        int green = background[idx + 1];
        int blue  = background[idx + 2];

        color c = {red, green, blue};

        return c;
    }

  public:
    double aspect_ratio = 1.0;    // Ratio of image width over height
    int image_width     = 100;    // Rendered image width in pixel count
    int image_height;             // Rendered image height
    int samples_per_pixel = 10;   // Count of random samples for each pixel

    double vfov           = 90;                         // Vertical view angle (field of view)
    vec4<double> lookfrom = c_0;                        // Point camera is looking from
    vec4<double> lookat   = vec4(0.0, 0.0, 0.0, 0.0);   // Point camera is looking at
    vec4<double> vup      = vec4(0.0, 0.0, 0.0, 1.0);   // Camera-relative "up" direction

    // REQUIRES: width, height, channels, and img all correspond to the background image
    // MODIFIES: this
    // EFFECTS: renders and integrates all the geodesics
    unsigned char* render(int width, int height, int channels, unsigned char* img) {
        this->background  = img;
        background_width  = width;
        background_height = height;
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

                int r = std::clamp(int(pixel_color.r * pixel_samples_scale + 0.5), 0, 255);
                image_buffer[(j * image_width + i) * 3 + 0] = static_cast<unsigned char>(r);
                int g = std::clamp(int(pixel_color.g * pixel_samples_scale + 0.5), 0, 255);
                image_buffer[(j * image_width + i) * 3 + 1] = static_cast<unsigned char>(g);
                int b = std::clamp(int(pixel_color.b * pixel_samples_scale + 0.5), 0, 255);
                image_buffer[(j * image_width + i) * 3 + 2] = static_cast<unsigned char>(b);
            }
        }

        memcpy(new_img, image_buffer.data(), img_size);

        return new_img;
    }
};