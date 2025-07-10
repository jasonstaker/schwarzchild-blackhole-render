#include "math/relativity/relativity.hpp"
#include "renderer/camera.hpp"

#include <iostream>
#include <stdlib.h>

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"


#include <filesystem>

int main() {

    // Load Image
    int width, height, channels;
    unsigned char* img = stbi_load("assets/background.jpg", &width, &height, &channels, 0);

    // Camera
    camera cam;

    cam.aspect_ratio = 16.0 / 9.0;
    cam.image_width = 400;
    cam.samples_per_pixel = 500;
    cam.max_depth = 50;
    cam.vfov = 20;

    auto new_img = cam.render(width, height, channels, img);
    // Spacetime values
    auto sm            = get_schwarzchild_metric(c_0);
    auto camera_tetrad = get_tetrad(c_0);

    // Render

    // Load Image
    stbi_write_jpg("output.jpg", cam.image_width, cam.image_height, 3, new_img, 100);
    stbi_image_free(img);
}