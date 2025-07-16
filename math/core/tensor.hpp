#pragma once

#include <vector>

template<typename T> class tensor {
  private:
    std::vector<T> data;
    int width, height;

  public:
    tensor(int width, int height) : width(width), height(height), data(width * height) {}

    const T& operator()(int x, int y) const { return data[y * width + x]; }
    T& operator()(int x, int y) { return data[y * width + x]; }
};