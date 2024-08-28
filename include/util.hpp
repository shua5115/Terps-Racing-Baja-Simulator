#pragma once

#include "Eigen/Dense"

#define DEG2RAD   (0.0174532925199432957692369)
#define RAD2DEG   (57.295779513082320876798155)
#define RPM2RADPS (0.1047197551196597746154214)
#define IN2M (0.0254)

constexpr inline double lerp(double a, double b, double t) {
    return a*(1.0-t) + b*t;
}

constexpr inline double invlerp(double a, double b, double v) {
    return (v-a)/(b-a);
}

constexpr inline double remap(double v, double a1, double b1, double a2, double b2) {
    return lerp(a2, b2, invlerp(a1, b1, v));
}

namespace Eigen {
    template<class T>
    void swap(T&& a, T&& b){
        a.swap(b);
    }
}

inline void matrix_sort_rows(Eigen::MatrixXd mat, Eigen::Index col) {
    std::sort(mat.rowwise().begin(), mat.rowwise().end(),
      [col](auto const& r1, auto const& r2){return r1(col)<r2(col);});
}

inline void matrix_sort_cols(Eigen::MatrixXd mat, Eigen::Index row) {
    std::sort(mat.colwise().begin(), mat.colwise().end(),
      [row](auto const& c1, auto const& c2){return c1(row)<c2(row);});
}

// 1-d linear interpolation lookup
// Matrix should have first column sorted in increasing order
double matrix_linear_lookup(Eigen::MatrixX2d mat, double val, bool extrapolate = false);