#include <thread>
#include "trb.hpp"

typedef void (*test_fn)();

void primary_ramp();

test_fn tests[] = {
    primary_ramp,
};

int main() {
    std::vector<std::thread> threads;
    for(auto test : tests) {
        auto t = std::thread([test](){
            try {
                test();
            } catch(std::exception &e) {
                e.what();
            }
        });
    }
    for(auto &t : threads) {
        if (t.joinable())
            t.join();
    }
}

void primary_ramp() {
    

}