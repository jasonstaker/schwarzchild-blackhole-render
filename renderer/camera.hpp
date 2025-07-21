#pragma once

#include "../math/core/vec4.hpp"
#include "../math/relativity/geodesic.hpp"
#include "../math/relativity/relativity.hpp"

#include <algorithm>
#include <array>
#include <atomic>
#include <iomanip>
#include <iostream>
#include <omp.h>
#include <vector>

struct color {
    int r, g, b;
};

// camera class specialized for black hole rendering
class camera {
  private:
    double pixel_samples_scale;   // color scale factor for a sum of pixel samples
    vec4<double> center;          // camera center
    vec4<double> pixel00_loc;     // location of pixel 0, 0
    vec4<double> pixel_delta_u;   // offset to pixel to the right
    vec4<double> pixel_delta_v;   // offset to pixel below
    vec4<double> u, v, w;         // camera frame basis vectors

    unsigned char* background;   // background image location
    int background_width;        // background image width
    int background_height;       // background image height
    int background_channels;     // background image channels

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

        geodesic g = {center, unit_vector(pixel_sample - center), 0.0};
        return g;
    }

    // REQUIRES: g's position vector is in schwarzcschild coordinates
    // EFFECTS: returns the pixel mapped from spherical coordinates to texture coordinates
    color geodesic_pixel(const geodesic& g) {
        auto s_position = g.position;

        double mag   = s_position.r();
        double theta = s_position.theta();
        double phi   = s_position.phi();

        if (s_position.r() <= R_S) { return {0, 0, 0}; }

        double raw_u = phi / (2 * PI);
        double raw_v = theta / PI;

        double u = raw_u - std::floor(raw_u);
        double v = std::clamp(raw_v, 0.0, 1.0);

        int x = int(u * (background_width - 1) + 0.5);
        int y = int(v * (background_height - 1) + 0.5);
        x     = std::clamp(x, 0, background_width - 1);
        y     = std::clamp(y, 0, background_height - 1);

        int idx = (y * background_width + x) * 3;
        return {background[idx + 0], background[idx + 1], background[idx + 2]};
    }

    // EFFECTS: integrates a geodesic by max_steps, enforcing the null condition each time
    geodesic run_geodesic(const geodesic& g_init) const {
        int max_steps = 5000;
        geodesic g    = {
          convert_to_schwarzschild_pos(g_init.position), convert_to_schwarzschild_vel(g_init.velocity), 0.0};

        g = enforce_null_condition(g);

        for (int step = 1; step <= max_steps; step++) {
            g = integrate(g);
            g = enforce_null_condition(g);

            if (check_termination(g)) { break; }
        }

        return g;
    }

  public:
    double aspect_ratio = 1.0;    // Ratio of image width over height
    int image_width     = 100;    // Rendered image width in pixel count
    int image_height;             // Rendered image height
    int samples_per_pixel = 10;   // Count of random samples for each pixel

    double vfov           = 90;                                // Vertical view angle (field of view)
    vec4<double> lookfrom = c_0;                               // Point camera is looking from
    vec4<double> lookat   = vec4(0.0, 0.0, 0.0, 0.0, false);   // Point camera is looking at
    vec4<double> vup      = vec4(0.0, 0.0, 0.0, 1.0, false);   // Camera-relative "up" direction

    // REQUIRES: width, height, channels, and img all correspond to the background image
    // MODIFIES: this
    // EFFECTS: renders and integrates all the geodesics
    unsigned char* render(int width, int height, int channels, unsigned char* img) {
        this->background    = img;
        background_width    = width;
        background_height   = height;
        background_channels = channels;
        initialize();

        size_t img_size = (image_width * image_height * background_channels);
        std::vector<unsigned char> image_buffer(img_size);

        unsigned char* new_img = static_cast<unsigned char*>(malloc(img_size));

        unsigned char* p = new_img;

        std::atomic<int> lines_done = 0;

#pragma omp parallel for schedule(dynamic)
        for (int j = 0; j < image_height; j++) {
            for (int i = 0; i < image_width; i++) {
                bool allEscaped = true;
                color sumColor  = {0, 0, 0};

                for (int sample = 0; sample < samples_per_pixel; ++sample) {
                    geodesic g_init = get_geodesic(i, j);
                    geodesic g      = run_geodesic(g_init);

                    if (g.position.r() <= R_S) {
                        allEscaped = false;
                        break;
                    }

                    color pix = geodesic_pixel(g);
                    sumColor.r += pix.r;
                    sumColor.g += pix.g;
                    sumColor.b += pix.b;
                }

                color finalPix = {0, 0, 0};
                if (allEscaped) {
                    double scale = 1.0 / samples_per_pixel;
                    finalPix.r   = std::clamp(int(sumColor.r * scale + 0.5), 0, 255);
                    finalPix.g   = std::clamp(int(sumColor.g * scale + 0.5), 0, 255);
                    finalPix.b   = std::clamp(int(sumColor.b * scale + 0.5), 0, 255);
                }

                image_buffer[(j * image_width + i) * background_channels + 0] = static_cast<unsigned char>(finalPix.r);
                image_buffer[(j * image_width + i) * background_channels + 1] = static_cast<unsigned char>(finalPix.g);
                image_buffer[(j * image_width + i) * background_channels + 2] = static_cast<unsigned char>(finalPix.b);
            }

            int done = ++lines_done;
            if (omp_get_thread_num() == 0)
                std::cout << "\rScanlines remaining: " << (image_height - done) << ' ' << std::flush;
        }

        memcpy(new_img, image_buffer.data(), img_size);

        return new_img;
    }
};