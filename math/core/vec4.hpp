#pragma once

#include <cmath>
#include <math.h>
#include <random>
#include <stdexcept>
#include <vector>

template<typename T> class vec4 {
  private:
    T e[4];

  public:
    bool schwarzschild;

    vec4() : e {0, 0, 0, 0} {}
    vec4(T t, T x, T y, T z, bool is_schwarzschild) : e {t, x, y, z}, schwarzschild(is_schwarzschild) {}

    T t() { return e[0]; }
    T t() const { return e[0]; }
    T x() { return e[1]; }
    T x() const { return e[1]; }
    T r() { return e[1]; }
    T r() const { return e[1]; }
    T y() { return e[2]; }
    T y() const { return e[2]; }
    T theta() { return e[2]; }
    T theta() const { return e[2]; }
    T z() { return e[3]; }
    T z() const { return e[3]; }
    T phi() { return e[3]; }
    T phi() const { return e[3]; }

    vec4<T> operator-() const { return vec4(-e[0], -e[1], -e[2], -e[3], schwarzschild); }
    const T& operator[](int i) const { return e[i]; }
    T& operator[](int i) { return e[i]; }

    double length_squared() const { return e[1] * e[1] + e[2] * e[2] + e[3] * e[3]; }

    double length() const { return std::sqrt(length_squared()); }
};

template<typename T> vec4<T> operator+(const vec4<T>& u, const vec4<T>& v) {
    if (u.schwarzschild != v.schwarzschild) {
        throw std::runtime_error("Error: Attempt to add vec4 with differing flags");
    }

    bool is_schwarzschild = u.schwarzschild;

    return vec4<T>(u[0] + v[0], u[1] + v[1], u[2] + v[2], u[3] + v[3], is_schwarzschild);
}

template<typename T> vec4<T> operator-(const vec4<T>& u, const vec4<T>& v) {
    if (u.schwarzschild != v.schwarzschild) {
        throw std::runtime_error("Error: Attempt to subtract vec4 with differing flags");
    }

    bool is_schwarzschild = u.schwarzschild;

    return vec4(u[0] - v[0], u[1] - v[1], u[2] - v[2], u[3] - v[3], is_schwarzschild);
}

template<typename T> vec4<T> operator*(const vec4<T>& u, double t) {
    bool is_schwarzschild = u.schwarzschild;

    return vec4(u[0] * t, u[1] * t, u[2] * t, u[3] * t, is_schwarzschild);
}

template<typename T> vec4<T> operator*(double t, const vec4<T>& u) {
    bool is_schwarzschild = u.schwarzschild;

    return vec4(u[0] * t, u[1] * t, u[2] * t, u[3] * t, is_schwarzschild);
}

template<typename T> vec4<T> operator/(const vec4<T>& u, double t) { return u * (1 / t); }

template<typename T> vec4<T> unit_vector(const vec4<T>& u) { return u / u.length(); }

template<typename T> vec4<T> cross(const vec4<T>& u, const vec4<T>& v) {
    if (u.schwarzschild != v.schwarzschild) {
        throw std::runtime_error("Error: Attempt to cross vec4 with differing flags");
    }

    bool is_schwarzschild = u.schwarzschild;

    T c1 = u[2] * v[3] - u[3] * v[2];
    T c2 = u[3] * v[1] - u[1] * v[3];
    T c3 = u[1] * v[2] - u[2] * v[1];

    return vec4<T>(0, c1, c2, c3, is_schwarzschild);
}

template<typename T> double dot(const vec4<T>& u, const vec4<T>& v) {
    if (u.schwarzschild != v.schwarzschild) { throw std::runtime_error("Schwarzschild flags differ in dot product"); }

    return v[1] * u[1] + v[2] * u[2] + v[3] * u[3];
}

inline double degrees_to_radians(double degrees) { return degrees * (M_PI / 180); }

inline double random_double() { return std::rand() / (RAND_MAX + 1.0); }