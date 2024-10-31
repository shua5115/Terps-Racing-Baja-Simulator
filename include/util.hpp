#pragma once

#include "Eigen/Dense"

#define PI        (3.1415926535897932384626434) // precision? yes.
#define DEG2RAD   (0.0174532925199432957692369) // rad/deg
#define RAD2DEG   (57.295779513082320876798155) // deg/rad
#define RPM2RADPS (0.1047197551196597746154214) // rpm-s/rad
#define RADPS2RPM (9.5492965855137201461330300) // rad/rpm-s
#define IN2M      (0.0254)                      // m/in
#define FT2M      (0.3048)                      // m/ft
#define LBF2N     (4.4482216153)                // lbf/N
#define SLUG2KG   (14.593902937)                // kg/slug
#define LBF2KG    (0.45359236844386)            // lbm/lbf * kg/lbm

// Positive modulo
template<typename T>
constexpr T imod(T i, T n) {
    return (n + (i % n)) % n;
}

constexpr double lerp(double a, double b, double t) {
    return a*(1.0-t) + b*t;
}

constexpr double invlerp(double a, double b, double v) {
    return (v-a)/(b-a);
}

constexpr double remap(double v, double a1, double b1, double a2, double b2) {
    return lerp(a2, b2, invlerp(a1, b1, v));
}

constexpr double clamp(double v, double lo, double hi) {
    return std::min(std::max(v, lo), hi);
}

constexpr double sign(double x) {
    if (x < 0) return -1;
    if (x > 0) return 1;
    return 0;
}

// hack to make std::sort work
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

// Finds the slope of an arbitrary 1-D function at input x numerically using two-sided finite difference.
template<typename F>
double diff_central(F f, double x, double step) {
    return (f(x+step) - f(x-step))*0.5/(step);
}

template<typename F>
double root_bisection(F f, double x_lo, double x_hi, size_t N) {
    double x_mid, y_lo, y_hi;
    y_lo = f(x_lo);
    y_hi = f(x_hi);
    for (size_t i = 0; i < N; i++) {
        x_mid = 0.5*(x_lo + x_hi);
        double y_mid = f(x_mid);
        double sign_left = sign(y_lo)*sign(y_mid);
        double sign_right = sign(y_hi)*sign(y_mid);
        if (sign_left < 0.0) {
            x_hi = x_mid;
        } else if (sign_right < 0.0) {
            x_lo = x_mid;
        } else if (sign_left != 0 && sign_right != 0){
            return x_mid;
        } else break;
    }
    x_mid = 0.5*(x_lo + x_hi);
    return x_mid;
}

template<typename F>
double root_secant(F f, double x0, double x1, unsigned int N) {
    double y0 = f(x0);
    double y1 = f(x1);
    for (size_t i = 0; i < N; i++) {
        if (abs(y1-y0) == 0) break;
        double x_next = x1 - y1*(x1-x0)/(y1-y0);
        x0 = x1;
        x1 = x_next;
        y0 = y1;
        y1 = f(x1);
    }
    return (x0+x1)*0.5;
}

template<typename F>
double root_newton(F f, double x, double diff_step, size_t N) {
    for (size_t i = 0; i < N; i++) {
        double dfdx = diff_central(f, x, diff_step);
        if (dfdx == 0) break;
        x = x - f(x)/dfdx;
    }
    return x;
}