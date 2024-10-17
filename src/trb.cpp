#include "trb.hpp"

OptResults<3> solve_flyweight_position(
    double theta1_guess, double theta2_guess, double d_r_guess,
    double (*ramp)(double x),
    double L1, double L2,
    double d_p, double x_ramp,
    double r_cage, double r_shoulder
) {
    Eigen::Vector3d x0(theta1_guess, theta2_guess, d_r_guess);

    ScalarFn<3> objective = [=](Eigen::Vector3d x){
        double theta1 = clamp(x(0), -PI*0.5, PI*0.5);
        double theta2 = clamp(x(1), -PI*0.5, PI*0.5);
        double d_r = std::max(x(2), 0.0);

        double ramp_eval = d_r - L2*cos(theta2);

        double eq1 = L1*cos(theta1) - x_ramp - d_p + d_r;
        double eq2 = L1*sin(theta1) + L2*sin(theta2) - r_cage + r_shoulder + ramp(ramp_eval);
        double eq3 = atan(diff_central(ramp, ramp_eval, 1e-4)) + PI*0.5 - theta2;

        return eq1*eq1 + eq2*eq2 + eq3*eq3;
    };

    return minimize_gradient_descent(objective, x0, 1e-5, 1e-6, 0.5, 50);
}