#include <iostream>
#include "csv.hpp"

int main() {
    auto data = read_csv(std::cin);
    for (const auto &row : data) {
        for (const auto &elem : row) {
            std::cout << '"' << elem << "\",";
        }
        std::cout << '\n';
    }
}