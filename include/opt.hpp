// Contains functions related to optimization

#include <functional>
#include "Eigen/Dense"

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

// Returns the gradient of the function f at x using second-order central finite difference.
// Evaluates the function 2N times, in case performance is a concern.
template<int N>
Eigen::Vector<double, N> gradient(ScalarFn<N> f, Eigen::Vector<double, N> x, double dx);

// Minimizes the objective function using gradient descent with golden section.
// Starts with initial guess x0.
// Gradient is computed with second-order central finite difference, using grad_step step size.
// Golden section halts when the interval size is less than tol. Golden section takes about 25-50 function calls per iteration of gradient descent for reasonable values of tol.
// Gradient descent halts when the size of the optimal step size is below tol, convergence depends on behavior of objective function.
// An conservative estimate of objective function evaluations for max_iter = 1000 is 1000*(2N+50).
template<int N>
OptResults<N> minimize_gradient_descent(ScalarFn<N> objective, Eigen::Vector<double, N> x0, double grad_step = 1e-4, double tol = 1e-6, size_t max_iter = 1000);