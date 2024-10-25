#include <iostream>
#include <chrono>
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

void belt_length() {
    BajaState state = TR24_GAGED_GX9;
    double tan_phi = tan(state.phi);
    double inv_tan_phi = 1.0/tan_phi;
    double r_s_min = state.r_s_min();
    double r_s_max = r_s_min + state.d_s_max*inv_tan_phi;
    double r_p_min = state.r_p_min();
    double r_p_max = r_p_min + state.d_p_max*inv_tan_phi;
    double L = state.L;
    double L0 = state.L_b0;

    size_t N = 50;
    for(size_t i = 0; i < N; i++) {
        double r_p = remap(i, 0, N-1, r_p_min, r_p_max);
        double r_s = solve_r_s(r_p, r_s_min, r_s_max, L, L0, 8);
        double alpha = 2*acos(clamp((r_s-r_p)/L, -1, 1));
        double beta = 2*PI-alpha;
        double e = r_p*alpha + r_s*beta + 2*sqrt(L*L - (r_p - r_s)*(r_p - r_s)) - L0;
        printf("r_p=%f, r_s=%f, error=%e\n", r_p, r_s, e);
    }
}

// void primary_roller_solver() {
//     auto start = std::chrono::high_resolution_clock::now();
//     double theta1_guess = 0;
//     double theta2_guess = 0;
//     double d_r_guess = 0;
//     auto ramp = [](double x){ return remap(x, 0, 1, 0.25, 0); };
//     double L1 = 1;
//     double L2 = 0.125;
//     double d_p = 0.9;
//     double x_ramp = 1;
//     double r_cage = 0.75;
//     double r_shoulder = 0.25;
    
//     auto S = solve_flyweight_position(
//         theta1_guess, theta2_guess, d_r_guess,
//         ramp,
//         L1, L2, d_p, x_ramp, r_cage, r_shoulder
//     );
//     auto finish = std::chrono::high_resolution_clock::now();
//     auto dur = finish - start;
//     printf("Finished in %f us\n", dur.count()*0.001);
//     printf("Solver %s after %llu iterations\n", S.converged ? "converged" : "did not converge", S.iterations);
//     printf("Error: %e = %e^2\n", S.f_of_x, sqrt(S.f_of_x));
//     printf("Solution: theta1=%f deg, theta2=%f deg, d_r=%f m]\n", S.x(0)*RAD2DEG, S.x(1)*RAD2DEG, S.x(2));
//     ASSERT(S.f_of_x <= 1e-6);
// }

void cvt_shift_solver() {
    BajaState state = TR24_GAGED_GX9;
    state.omega_p = 3000*RPM2RADPS;
    state.tau_s = 18.5*LBF2N/FT2M; // 18.5 lb-ft, max output from engine
    state.controls.throttle = 1;

    auto start = std::chrono::high_resolution_clock::now();
    auto S = solve_cvt_shift_no_slip(state);
    auto finish = std::chrono::high_resolution_clock::now();
    auto dur = finish - start;
    printf("Shift solver finished in %f us\n", dur.count()*0.001);
    printf("Solver %s after %llu iterations\n", S.converged ? "converged" : "did not converge", S.iterations);
    printf("Error: %e = %e^2\n", S.f_of_x, sqrt(S.f_of_x));
    double d_p = S.x(0);
    double r_p = d_p/tan(state.phi) + state.r_p_min();
    double r_s = solve_r_s(r_p, state.r_s_min(), state.r_s_min() + state.d_s_max/tan(state.phi), state.L, state.L_b0, 8);
    double d_s = state.d_s_max - (r_s - state.r_s_min())*tan(state.phi);
    printf("Solution: d_p=%f m, d_s=%f m, r_p=%f, r_s=%f, ratio=%f\n", d_p, d_s, r_p, r_s, r_s/r_p);
    ASSERT(S.converged && (S.f_of_x < 1e-3));
}

const TestCase tests[] = {
    TESTFN(engine_rpm_lookup),
    // TESTFN(primary_roller_solver),
    TESTFN(belt_length),
    TESTFN(cvt_shift_solver),
};

int main() {
    constexpr size_t N_tests = sizeof(tests)/sizeof(TestCase);
    size_t N_success = 0;
    for(const auto &testcase : tests) {
        try {
            std::cout << "\nRUNNING TEST: \"" << testcase.name << "\":\n";
            testcase.test();
            std::cout << "TEST SUCCEEDED: \"" << testcase.name << "\"\n";
            N_success++;
        } catch(const std::exception &e) {
            std::cout << e.what() << "\n";
            std::cout << "TEST FAILED: \"" << testcase.name << "\"\n";
        }
    }
    std::cout << "\nRESULTS: " << N_success << "/" << N_tests <<" tests succeeded, " << (N_tests - N_success) << "/" << N_tests << " tests failed.\n";
}