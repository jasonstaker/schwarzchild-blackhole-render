#include "../math/core/tensor.hpp"
#include "../math/core/vec4.hpp"
#include "../math/relativity/relativity.hpp"

#include <cmath>
#include <gtest/gtest.h>
#include <numeric>

using namespace relativity;
using namespace constants;
using namespace parameters;

static constexpr double eps = 1e-6;

double lightlike_invariant(const geodesic& g) {
    auto M     = get_schwarzschild_metric(g.position);
    double sum = 0;
    for (int μ = 0; μ < 4; ++μ)
        for (int ν = 0; ν < 4; ++ν) sum += M(μ, ν) * g.velocity[μ] * g.velocity[ν];
    return sum;
}

TEST(Relativity, TetradOrthonormality) {
    auto tet = get_tetrad();
    ASSERT_EQ(tet.size(), 4u);

    for (int i = 1; i < 4; ++i) {
        for (int j = 1; j < 4; ++j) {
            double d = dot(tet[i], tet[j]);
            if (i == j) {
                ASSERT_GT(std::abs(d), 0.0);
            } else {
                ASSERT_NEAR(d, 0.0, eps);
            }
        }
    }
}

TEST(Relativity, ConvertCartesianToSchwarzschildAndBack) {
    vec4<double> c0(0.5, 3.0, -4.0, 2.0, false);
    auto s = convert_to_schwarzschild_pos(c0);
    EXPECT_TRUE(s.schwarzschild);

    double r = std::sqrt(3 * 3 + 4 * 4 + 2 * 2);
    EXPECT_NEAR(s.r(), r, eps);
    EXPECT_NEAR(s.theta(), std::acos(2.0 / r), eps);
    EXPECT_NEAR(s.phi(), std::atan2(-4.0, 3.0), eps);

    auto c1 = convert_to_cartesian_pos(s);
    EXPECT_FALSE(c1.schwarzschild);
    EXPECT_NEAR(c1.x(), c0.x(), eps);
    EXPECT_NEAR(c1.y(), c0.y(), eps);
    EXPECT_NEAR(c1.z(), c0.z(), eps);
}

TEST(Relativity, ConvertVelocityToSchwarzschild) {
    vec4<double> v_cart(1.0, 0.1, 0.2, 0.3, false);
    auto v_sch = convert_to_schwarzschild_vel(v_cart);
    EXPECT_TRUE(v_sch.schwarzschild);

    EXPECT_NE(v_sch[0], v_cart[0]);
}

TEST(Relativity, MetricComponentsAtVariousR) {
    for (auto rr : {R_0, R_0 * 2, R_0 / 2}) {
        vec4<double> pos(0, rr, M_PI / 4, M_PI / 3, true);
        auto g = get_schwarzschild_metric(pos);
        EXPECT_NEAR(g(0, 0), -(1 - R_S / rr), eps);
        EXPECT_NEAR(g(1, 1), 1.0 / (1 - R_S / rr), eps);
        EXPECT_NEAR(g(2, 2), rr * rr, eps);
        EXPECT_NEAR(g(3, 3), rr * rr * std::sin(pos.theta()) * std::sin(pos.theta()), eps);

        EXPECT_NEAR(g(0, 1), 0.0, eps);
        EXPECT_NEAR(g(2, 3), 0.0, eps);
    }
}

TEST(Relativity, ChristoffelSpecialAndGeneral) {
    vec4<double> u(0, 10.0, M_PI / 3, M_PI / 6, true);

    double g01 = get_christoffel(0, 1, 0, u);
    double g10 = get_christoffel(0, 0, 1, u);
    EXPECT_NEAR(g01, g10, eps);
    EXPECT_NEAR(g01, M / (10.0 * (10.0 - 2.0 * M)), eps);

    EXPECT_NEAR(get_christoffel(1, 1, 1, u), -(M / (10.0 * (10.0 - 2 * M))), eps);

    EXPECT_NEAR(get_christoffel(2, 0, 3, u), 0.0, eps);
    EXPECT_NEAR(get_christoffel(3, 2, 0, u), 0.0, eps);
}

TEST(Relativity, RK4StepConservesNullConditionApprox) {

    geodesic g0 {{1.0, R_0, M_PI / 2, 0, true}, {1.0, 0.1, 0.1, 0.1, true}, 0.0};
    auto g1 = rk4_step(g0);

    EXPECT_NE(g1.velocity[1], g0.velocity[1]);
}

TEST(Relativity, EnforceNullConditionProducesLightlike) {
    geodesic g0 {{1.0, R_0, M_PI / 2, 0, true}, {0.0, 1.0, 2.0, 3.0, true}, 0.0};
    auto g1 = enforce_null_condition(g0);
    EXPECT_NEAR(lightlike_invariant(g1), 0.0, 1e-5);
}

TEST(Relativity, IntegrateMovesRay) {
    geodesic g {{0, 20.0, M_PI / 2, 0, true}, {1.0, -0.5, 0, 0, true}, 0.0};
    auto g2 = integrate(g);
    EXPECT_LT(g2.position.r(), g.position.r());
}

TEST(Relativity, CheckTerminationBehavior) {
    geodesic inside {{0, 1.0, 0, 0, true}, {1, 0, 0, 0, true}, 0.0};
    geodesic escape {{0, R_ESCAPE + 1, 0, 0, true}, {1, 0, 0, 0, true}, 0.0};
    geodesic normal {{0, (R_S + R_ESCAPE) / 2, 0, 0, true}, {1, 0, 0, 0, true}, 0.0};
    EXPECT_TRUE(check_termination(inside));
    EXPECT_TRUE(check_termination(escape));
    EXPECT_FALSE(check_termination(normal));
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
