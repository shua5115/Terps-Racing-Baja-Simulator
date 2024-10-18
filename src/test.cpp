#include <iostream>
#include "trb.hpp"

// bs to define an assert macro
#define MACRO_STR(x) #x
#define MACRO_STR2(x) MACRO_STR(x)
#define ASSERT(cond) test_assert((cond) , "ASSERT FAILED: " #cond " in \"" __FILE__ "\" at line " MACRO_STR2(__LINE__))
void test_assert(bool condition, const char *what) {
    if (condition) return;
    throw std::runtime_error(what);
}

#define TESTFN(name) {#name, name}
typedef void (*test_fn)();

struct TestCase {const char *name; test_fn test;};

void engine_rpm_lookup() {
    for(double w = -4000; w <= 4000; w += 500) {
        double torque = matrix_linear_lookup(CH440_TORQUE_CURVE, w);
        printf("%f\n", torque);
        if (w < 1800) ASSERT(torque == 25.4);
        if (w > 3800) ASSERT(torque == 23.9);
    }
}

void primary_roller_solver() {
    double theta1_guess = 0;
    double theta2_guess = 0;
    double d_r_guess = 0;
    auto ramp = [](double x){ return remap(x, 0, 1, 0.25, 0); };
    double L1 = 1;
    double L2 = 0.125;
    double d_p = 0.9;
    double x_ramp = 1;
    double r_cage = 0.75;
    double r_shoulder = 0.25;
    
    auto S = solve_flyweight_position(
        theta1_guess, theta2_guess, d_r_guess,
        ramp,
        L1, L2, d_p, x_ramp, r_cage, r_shoulder
    );

    printf("Solver %s after %llu iterations\n", S.converged ? "converged" : "did not converge", S.iterations);
    printf("Error: %e = %e^2\n", S.f_of_x, sqrt(S.f_of_x));
    printf("Solution: theta1=%f deg, theta2=%f deg, d_r=%f m]\n", S.x(0)*RAD2DEG, S.x(1)*RAD2DEG, S.x(2));
    ASSERT(S.f_of_x <= 1e-6);
}

void cvt_shift_solver() {

}

TestCase tests[] = {
    TESTFN(engine_rpm_lookup),
    TESTFN(primary_roller_solver),
    TESTFN(cvt_shift_solver),
};

int main() {
    for(const auto &testcase : tests) {
        try {
            std::cout << "\nRUNNING TEST: \"" << testcase.name << "\":\n";
            testcase.test();
            std::cout << "TEST SUCCEEDED: \"" << testcase.name << "\"\n";
        } catch(const std::exception &e) {
            std::cerr << "TEST FAILED: \"" << testcase.name << "\":\n" << e.what();
        }
    }
}