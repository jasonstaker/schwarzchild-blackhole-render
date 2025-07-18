#include "../math/core/vec4.hpp"

#include <cmath>
#include <gtest/gtest.h>
#include <stdexcept>

using Vec                   = vec4<double>;
static constexpr double eps = 1e-9;

TEST(Vec4, DefaultConstructor) {
    Vec v;
    EXPECT_FALSE(v.schwarzschild);
    for (int i = 0; i < 4; ++i) EXPECT_EQ(v[i], 0.0);
}

TEST(Vec4, ValueConstructorAndAccessors) {
    Vec v(1.0, 2.0, 3.0, 4.0, true);
    EXPECT_TRUE(v.schwarzschild);
    EXPECT_DOUBLE_EQ(v.t(), 1.0);
    EXPECT_DOUBLE_EQ(v.x(), 2.0);
    EXPECT_DOUBLE_EQ(v.y(), 3.0);
    EXPECT_DOUBLE_EQ(v.z(), 4.0);

    EXPECT_DOUBLE_EQ(v[0], 1.0);
    EXPECT_DOUBLE_EQ(v[1], 2.0);
    EXPECT_DOUBLE_EQ(v[2], 3.0);
    EXPECT_DOUBLE_EQ(v[3], 4.0);
}

TEST(Vec4, UnaryMinusPreservesFlag) {
    Vec v(1, -2, 3, -4, true);
    Vec u = -v;
    EXPECT_TRUE(u.schwarzschild);
    EXPECT_DOUBLE_EQ(u[0], -1.0);
    EXPECT_DOUBLE_EQ(u[1], 2.0);
    EXPECT_DOUBLE_EQ(u[2], -3.0);
    EXPECT_DOUBLE_EQ(u[3], 4.0);
}

TEST(Vec4, AdditionAndSubtraction) {
    Vec a(1, 2, 3, 4, false), b(0.5, 1, 1, 1, false);
    Vec c = a + b;
    EXPECT_DOUBLE_EQ(c[0], 1.5);
    EXPECT_DOUBLE_EQ(c[1], 3.0);
    EXPECT_DOUBLE_EQ(c[2], 4.0);
    EXPECT_DOUBLE_EQ(c[3], 5.0);
    EXPECT_FALSE(c.schwarzschild);

    Vec d = a - b;
    EXPECT_DOUBLE_EQ(d[0], 0.5);
    EXPECT_DOUBLE_EQ(d[1], 1.0);
    EXPECT_DOUBLE_EQ(d[2], 2.0);
    EXPECT_DOUBLE_EQ(d[3], 3.0);
    EXPECT_FALSE(d.schwarzschild);
}

TEST(Vec4, AdditionThrowsOnMixedFlag) {
    Vec a(0, 0, 0, 0, false), b(0, 0, 0, 0, true);
    EXPECT_THROW(a + b, std::runtime_error);
    EXPECT_THROW(a - b, std::runtime_error);
}

TEST(Vec4, ScalarMultiplyDividePreserveFlag) {
    Vec v(1, 2, 3, 4, true);
    Vec a = v * 2.0;
    Vec b = 2.0 * v;
    Vec c = v / 2.0;
    EXPECT_TRUE(a.schwarzschild && b.schwarzschild && c.schwarzschild);
    EXPECT_DOUBLE_EQ(a[1], 4.0);
    EXPECT_DOUBLE_EQ(b[2], 6.0);
    EXPECT_DOUBLE_EQ(c[3], 2.0);
}

TEST(Vec4, LengthAndSquared) {
    Vec v(0, 3, 4, 0, false);
    EXPECT_DOUBLE_EQ(v.length_squared(), 25.0);
    EXPECT_DOUBLE_EQ(v.length(), 5.0);
}

TEST(Vec4, UnitVector) {
    Vec v(0, 3, 4, 0, false);
    Vec u = unit_vector(v);
    EXPECT_DOUBLE_EQ(u.length(), 1.0);
    EXPECT_NEAR(u[1], 3.0 / 5.0, eps);
    EXPECT_NEAR(u[2], 4.0 / 5.0, eps);
}

TEST(Vec4, DotAndCross) {
    Vec u(0, 1, 0, 0, false), v(0, 0, 1, 0, false);
    double d = dot(u, v);
    EXPECT_DOUBLE_EQ(d, 0.0);

    Vec w = cross(u, v);

    EXPECT_DOUBLE_EQ(w[1], 0.0);
    EXPECT_DOUBLE_EQ(w[2], 0.0);
    EXPECT_DOUBLE_EQ(w[3], 1.0);
    EXPECT_FALSE(w.schwarzschild);

    Vec a(0, 1, 0, 0, false), b(0, 1, 0, 0, true);
    EXPECT_THROW(dot(a, b), std::runtime_error);
    EXPECT_THROW(cross(a, b), std::runtime_error);
}
