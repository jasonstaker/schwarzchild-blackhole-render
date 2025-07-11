#pragma once

#include <cmath>
#include <math.h>
#include <vector>
#include <random>

template<typename T> class vec4 {
  private:
    T e[4];

  public:
    vec4() : e {0, 0, 0, 0} {}
    vec4(T t, T x, T y, T z) : e {t, x, y, z} {}

    T x() { return e[1]; }
    T r() { return e[1]; }
    T y() { return e[2]; }
    T theta() { return e[2]; }
    T z() { return e[3]; }
    T phi() { return e[3]; }

    vec4<T> operator-() const { return vec4(-e[0], -e[1], -e[2], -e[3]); }
    const T& operator[](int i) const { return e[i]; }
    T& operator[](int i) { return e[i]; }

    double length_squared() const { return e[0] * e[0] + e[1] * e[1] + e[2] * e[2] + e[3] * e[3]; }

    double length() const { return std::sqrt(length_squared()); }
};

template<typename T> vec4<T> operator+(const vec4<T>& u, const vec4<T>& v) {
    return vec4(u[0] + v[0], u[1] + v[1], u[2] + v[2], u[3] + v[3]);
}

template<typename T> vec4<T> operator-(const vec4<T>& u, const vec4<T>& v) {
    return vec4(u[0] - v[0], u[1] - v[1], u[2] - v[2], u[3] - v[3]);
}

template<typename T> vec4<T> operator*(const vec4<T>& u, double t) {
    return vec4(u[0] * t, u[1] * t, u[2] * t, u[3] * t);
}

template<typename T> vec4<T> operator*(double t, const vec4<T>& u) {
    return vec4(u[0] * t, u[1] * t, u[2] * t, u[3] * t);
}

template<typename T> vec4<T> operator/(const vec4<T>& u, double t) { return u * (1 / t); }

template<typename T> vec4<T> unit_vector(const vec4<T>& u) { return u / u.length(); }

template<typename T> vec4<T> cross(const vec4<T>& u, const vec4<T>& v) {
    T c1 = u[2] * v[3] - u[3] * v[2];
    T c2 = u[3] * v[1] - u[1] * v[3];
    T c3 = u[1] * v[2] - u[2] * v[1];

    return vec4<T>(0, c1, c2, c3);
}

template<typename T> double dot(const vec4<T>& u, const vec4<T>& v) {
    return v[0] * u[0] + v[1] * u[1] + v[2] * u[2] + v[3] * u[3];
}

inline double degrees_to_radians(double degrees) { return degrees * (M_PI / 180); }

inline double random_double() {
    return std::rand() / (RAND_MAX + 1.0);
}