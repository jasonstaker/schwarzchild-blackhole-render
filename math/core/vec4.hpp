#pragma once

#include <vector>

template<typename T> class vec4 {
  private:
    T e[4];

  public:
    vec4() : e {0, 0, 0, 0} {}
    vec4(T t, T x, T y, T z) : e {t, x, y, z} {}

    T lambda() { return e[0]; }
    T x() { return e[1]; }
    T r() { return e[1]; }
    T y() { return e[2]; }
    T theta() { return e[2]; }
    T z() { return e[3]; }
    T phi() { return e[3]; }

    const T& operator[](int i) const { return e[i]; }
    T& operator[](int i) { return e[i]; }
};