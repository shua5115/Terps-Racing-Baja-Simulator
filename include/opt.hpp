#pragma once
// Contains functions related to optimization
#include <iostream>
#include <functional>
#include "Eigen/Dense"

#define INV_GOLDEN (0.6180339887498948482)

// Function which takes vector of double and returns double.
template<int N>
using ScalarFn = std::function<double(Eigen::Vector<double, N>)>;

template<int N>
struct OptResults {
    Eigen::Vector<double, N> x; // x_optimum
    double f_of_x; // value of f(x_optimum)
    size_t iterations;
    bool converged;
};

template<int N>
static double opt_golden_section(ScalarFn<N> &f, Eigen::Vector<double, N> x, Eigen::Vector<double, N> d, double end_tol, double minval, double maxval) {
    double a = minval;
    double b = maxval;
    double alpha1 = a + (1.0-INV_GOLDEN)*(b-a);
    double alpha2 = a + INV_GOLDEN*(b-a);
    Eigen::Vector<double, N> x1, x2;
    for (size_t i = 0; i < N; i++) {
        x1(i) = x(i) + alpha1*d(i);
        x2(i) = x(i) + alpha2*d(i);
    }
    double f1 = f(x1);
    double f2 = f(x2);
    // Reusing function results to reduce calls to 1 per iteration
    while ((b-a) > end_tol) {
        // std::cout << "GOLDEN: x1=" << Eigen::Transpose(x1) << ", x2=" << Eigen::Transpose(x2) << ", section=[" << a << ", " << b << "]\n";
        if (f1 > f2) {
            // b stays the same
            a = alpha1;
            alpha1 = alpha2;
            x1 = x2;
            f1 = f2;
            alpha2 = a + INV_GOLDEN*(b-a);
            x2 = x + alpha2*d;
            f2 = f(x2);
        } else { // f1 < f2
            // a stays the same
            b = alpha2;
            alpha2 = alpha1;
            x2 = x1;
            f2 = f1;
            alpha1 = a + (1-INV_GOLDEN)*(b-a);
            x1 = x + alpha1*d;
            f1 = f(x1);
        }
    }
    double t1 = f2/(f1+f2);
    double t2 = f1/(f1+f2);
    return a*t1 + b*t2; // Return weighted average
}

// Returns the gradient of the function f at x using second-order central finite difference.
// Evaluates the function 2N times, in case performance is a concern.
template<int N>
Eigen::Vector<double, N> gradient(ScalarFn<N> f, Eigen::Vector<double, N> x, double dx) {
    Eigen::Vector<double, N> grad;
    Eigen::Vector<double, N> x_adj = x;
    for (size_t i = 0; i < N; i++) {
        double f1, f2;
        x_adj(i) = x(i) - dx;
        f1 = f(x_adj);
        x_adj(i) = x(i) + dx;
        f2 = f(x_adj);
        x_adj(i) = x(i); // reset this index
        grad(i) = (f2 - f1)*0.5/dx; // 2nd-order central difference
    }
    return grad;
}

// Minimizes the objective function "f" using gradient descent with golden section.
// Starts with initial guess x0.
// Gradient is computed with second-order central finite difference, using grad_step step size.
// Golden section halts when the interval size is less than tol. Initial interval size is called "section size". A larger section size can converge faster, but is less stable.
// Gradient descent halts when the size of the optimal step size is below tol, convergence depends on behavior of objective function.
// An conservative estimate of objective function evaluations for max_iter = 1000 is 1000*(2N+50).
template<int N>
OptResults<N> minimize_gradient_golden(ScalarFn<N> f, Eigen::Vector<double, N> x, const Eigen::Vector<double, N> x_lb, const Eigen::Vector<double, N> x_ub, double tol = 1e-6, double grad_step = 1e-6, double section_size = 1, size_t max_iter = 1000, bool scaled_step = false) {
    // Initially constrain x
    for(Eigen::Index i = 0; i < N; i++) {
        x(i) = clamp(x(i), x_lb(i), x_ub(i));
    }
    double f_val = f(x);
    Eigen::Vector<double, N> grad = gradient<N>(f, x, grad_step);
    double alpha;
    size_t iters = 0;
    bool converged = true;
    while (true) {
        // std::cout << "Iteration " << iters << ": x=[" << Eigen::Transpose(x) << "], f=" << f_val << "\n";
        bool isNaN = false;
        for (Eigen::Index i = 0; i < N; i++) {
            if (std::isnan(x(i)) || std::isnan(grad(i))) {
                isNaN = true;
                break;
            }
        }
        if (iters >= max_iter || std::isnan(f_val) || isNaN) {
            converged = false;
            break;
        }
        iters += 1;
        if (!scaled_step) grad.normalize();
        alpha = opt_golden_section<N>(f, x, grad, tol, -section_size, 0); // because gradient is kept positive, alpha must be <= 0
        grad *= alpha;
        // double step_len = grad.norm();
        Eigen::Vector<double, N> new_x = x + grad;
        for(Eigen::Index i = 0; i < N; i++) {
            new_x(i) = clamp(new_x(i), x_lb(i), x_ub(i));
        }
        double step_len = (new_x - x).norm();
        
        x = new_x;
        double new_f_val = f(x);
        if (step_len < tol) { // If the step size is small enough, then we have reached a local minimum!
            break;
        }
        f_val = new_f_val;
        grad = gradient<N>(f, x, grad_step);
    }
    // std::cout << "Result after iter " << iters << ": x=[" << Eigen::Transpose(x) << "], f=" << f_val << "\n";
    OptResults<N> results;
    results.x = x;
    results.f_of_x = f_val;
    results.iterations = iters;
    results.converged = converged;
    return results;
}

template<int N>
OptResults<N> minimize_gradient_simple(ScalarFn<N> f, Eigen::Vector<double, N> x, const Eigen::Vector<double, N> x_lb, const Eigen::Vector<double, N> x_ub, double tol = 1e-6, double grad_step = 1e-6, double step_scale = 1, size_t max_iter = 1000) {
    // Initially constrain x
    for(Eigen::Index i = 0; i < N; i++) {
        x(i) = clamp(x(i), x_lb(i), x_ub(i));
    }
    double f_val = f(x);
    Eigen::Vector<double, N> grad = gradient(f, x, grad_step);
    double alpha;
    size_t iters = 0;
    bool converged = true;
    while (true) {
        bool isNaN = false;
        for (Eigen::Index i = 0; i < N; i++) {
            if (std::isnan(x(i)) || std::isnan(grad(i))) {
                isNaN = true;
                break;
            }
        }
        if (iters >= max_iter || std::isnan(f_val) || isNaN) {
            converged = false;
            break;
        }
        iters += 1;
        grad *= -step_scale;
        double step_len = grad.norm();
        x += grad;
        for(Eigen::Index i = 0; i < N; i++) {
            x(i) = clamp(x(i), x_lb(i), x_ub(i));
        }
        double new_f_val = f(x);
        // printf("Iteration %llu: x=[%f, %f, %f], f=%f\n", iters, x(0), x(1), x(2), f_val);
        if (step_len < tol) { // If the step size is small enough, then we have reached a local minimum!
        // if (abs(new_f_val - f_val)/step_len < tol) { // if df/dx is very close to zero, then it is a local minimum!
            break;
        }
        f_val = new_f_val;
        grad = gradient(f, x, grad_step);
    }
    OptResults<N> results;
    results.x = x;
    results.f_of_x = f_val;
    results.iterations = iters;
    results.converged = converged;
    return results;
}