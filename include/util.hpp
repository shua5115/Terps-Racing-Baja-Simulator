#pragma once

#include <tuple>
#include "Eigen/Dense"

#ifndef PI
#define PI        (3.1415926535897932384626434) // precision? yes.
#endif
#ifndef DEG2RAD
#define DEG2RAD   (0.0174532925199432957692369) // rad/deg
#endif
#ifndef RAD2DEG
#define RAD2DEG   (57.295779513082320876798155) // deg/rad
#endif
#ifndef RPM2RADPS
#define RPM2RADPS (0.1047197551196597746154214) // rpm-s/rad
#endif
#ifndef RADPS2RPM
#define RADPS2RPM (9.5492965855137201461330300) // rad/rpm-s
#endif
#ifndef IN2M
#define IN2M      (0.0254)                      // m/in
#endif
#ifndef FT2M
#define FT2M      (0.3048)                      // m/ft
#endif
#ifndef LBF2N
#define LBF2N     (4.4482216153)                // lbf/N
#endif
#ifndef SLUG2KG
#define SLUG2KG   (14.593902937)                // kg/slug
#endif
#ifndef LBF2KG
#define LBF2KG    (0.45359236844386)            // lbm/lbf * kg/lbm
#endif
#ifndef METERPS2MPH
#define METERPS2MPH (2.2369362921)              // mi-s/m-h
#endif

// Positive modulo for types with operator% defined
template<typename T>
constexpr T imod(T i, T n) {
    return (n + (i % n)) % n;
}

// Positive modulo for double
inline double mod_euclid(double i, double n) {
    return fmod(n + fmod(i, n), n);
}

constexpr double Lerp(double a, double b, double t) {
    return a*(1.0-t) + b*t;
}

constexpr double invLerp(double a, double b, double v) {
    return (v-a)/(b-a);
}

constexpr double remap(double v, double a1, double b1, double a2, double b2) {
    return Lerp(a2, b2, invLerp(a1, b1, v));
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

// Sorts rows of a matrix in increasing order
inline void matrix_sort_rows(Eigen::MatrixXd mat, Eigen::Index col) {
    std::sort(mat.rowwise().begin(), mat.rowwise().end(),
      [col](auto const& r1, auto const& r2){return r1(col)<r2(col);});
}

// Sorts columns of a matrix in increasing order
inline void matrix_sort_cols(Eigen::MatrixXd mat, Eigen::Index row) {
    std::sort(mat.colwise().begin(), mat.colwise().end(),
      [row](auto const& c1, auto const& c2){return c1(row)<c2(row);});
}

// 1-d linear interpolation lookup
// Matrix should have first column sorted in increasing order
double matrix_linear_lookup(Eigen::MatrixX2d mat, double val, bool extrapolate = false);

// Finds the minimum and maximum of every entry in the dynamic matrix.
std::tuple<double, double> minmax(const Eigen::MatrixXd &m);

// Finds the slope of an arbitrary 1-D function at input x numerically using two-sided finite difference.
template<typename F>
double diff_central(F f, double x, double step) {
    return (f(x+step) - f(x-step))*0.5/(step);
}

// Finds the root of a function f in between bounds x_left and x_right.
// The bounds must contain a sign change for this to work.
// Error = (x_right-x_left)*2^(-N)
template<typename F>
double root_bisection(F f, double x_left, double x_right, size_t N) {
    double x_mid, y_left, y_right;
    y_left = f(x_left);
    y_right = f(x_right);
    if (y_left == 0) return x_left;
    if (y_right == 0) return x_right;
    for (size_t i = 0; i < N; i++) {
        x_mid = 0.5*(x_left + x_right);
        double y_mid = f(x_mid);
        if (y_mid == 0) return x_mid;
        double sign_left = sign(y_left)*sign(y_mid);
        double sign_right = sign(y_right)*sign(y_mid);
        if (sign_left < 0.0) {
            x_right = x_mid;
        } else if (sign_right < 0.0) {
            x_left = x_mid;
        } else {
            return x_mid;
        }
    }
    x_mid = 0.5*(x_left + x_right);
    return x_mid;
}

// Finds the root of a function f using the secant method.
// This uses a window of two sample points to approximate the slope.
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

// Finds the root of a function f using Newton's method.
// This uses 2nd-order central finite difference to approximate the slope.
template<typename F>
double root_newton(F f, double x, double diff_step, size_t N) {
    for (size_t i = 0; i < N; i++) {
        double dfdx = diff_central(f, x, diff_step);
        if (dfdx == 0) break;
        x = x - f(x)/dfdx;
    }
    return x;
}

// Integrates function f using 4th order explicit Runge Kutta integration for a single time step.
// Function f is of the form f(x, t)
// The integration solves dx/dt = f(x, t)
template<typename T, int N>
Eigen::Vector<T, N> runge_kutta_4_step(std::function<Eigen::Vector<T, N>(Eigen::Vector<T, N>, double)> f, Eigen::Vector<T, N> x, double t, double dt) {
    using V = Eigen::Vector<T, N>;
    using Index = Eigen::Index;
    V f1 = f(x, t);
    V f2 = f(x + 0.5*dt*f1, t + 0.5*dt);
    V f3 = f(x + 0.5*dt*f2, t + 0.5*dt);
    V f4 = f(x + dt*f3, t + dt);
    return x + dt*(f1/6.0 + f2/3.0 + f3/3.0 + f4/6.0);
}

// Performs integration of a multivariable function f using 4th order explicit Runge Kutta method over a range of time using fixed time steps.
// Each intermediate value of x is stored in a table.
// Each row of the table contains the time value in the first column, then x as a row vector for the remaining columns.
// Function f is of the form f(x, t)
// The integration solves dx/dt = f(x, t)
template<typename T, int N>
Eigen::Matrix<T, -1, N+1> runge_kutta_4(std::function<Eigen::Vector<T, N>(Eigen::Vector<T, N>, double)> f, Eigen::Vector<T, N> x0, double t0, double tf, double dt) {
    using V = Eigen::Vector<T, N>;
    using Index = Eigen::Index;
    V x = x0;
    double t = t0;
    Index M = (Index) ((tf - t0)/dt) + 1;
    Eigen::Matrix<T, -1, N+1> res;
    if (M <= 0) return res;
    res.resize(M, N+1);
    for(Index i = 0; i < M; i++) {
        res(i, 0) = t;
        for(Index j = 0; j < N; j++) {
            res(i, j+1) = x(j);
        }
        double actual_dt = std::min(t+dt, tf) - t;
        x = runge_kutta_4_step(f, x, t, actual_dt);
        t = std::min(t+dt, tf);
    }
    return res;
}