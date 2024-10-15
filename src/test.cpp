#include <iostream>
#include "trb.hpp"

#define TESTFN(name) {#name, name}

typedef void (*test_fn)();


void primary_ramp() {

}

struct {const char *name; test_fn test;} tests[] = {
    TESTFN(primary_ramp),
};

int main() {
    for(auto testcase : tests) {
        try {
            std::cout << "\nRUNNING TEST \"" << testcase.name << "\":\n";
            testcase.test();
        } catch(std::exception &e) {
            std::cerr << e.what();
        }
    }
}