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
    
}

TestCase tests[] = {
    TESTFN(engine_rpm_lookup),
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