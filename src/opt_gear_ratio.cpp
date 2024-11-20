// The gear ratio optimizer does not use CVT dynamics.

#include <iostream>
#include "trb.hpp"

#define sim_dt (1e-4)

void sim_step(BajaState &baja, double dt, double fixed_engine_rpm, double hill_deg, Eigen::Vector2d ratio_limits) {
    BajaDynamicsResult res = {0};
    double fixed_omega_p = fixed_engine_rpm*RPM2RADPS;
    baja.theta_hill = hill_deg*DEG2RAD;
    // To integrate velocity and position, we can use this definition: u = [x, x'], where u' = [x', x'']
    auto integrand = [&](Eigen::Vector2d u, double t){
        auto du = Eigen::Vector2d();
        double v = u(1);
        // Calculate the ratio required to keep engine at fixed rpm
        double omega_wheel = v/baja.r_wheel;
        double omega_s = omega_wheel*baja.N_g;
        double ratio = clamp(fixed_omega_p/omega_s, ratio_limits(0), ratio_limits(1));
        if (!isfinite(ratio)) ratio = ratio_limits(1);
        double tau_e = matrix_linear_lookup(baja.engine_torque_curve, fixed_engine_rpm);
        double tau_s = tau_e*ratio;

        // dx/dt = x'
        du(0) = u(1);
        // dx'/dt = acceleration of the car
        // Resistance forces on the car always act opposite to direction of v. This is captured with -sign(u(1)).
        du(1) = (tau_s*baja.N_g/(baja.r_wheel) - sign(u(1))*baja.F_total_resist() - baja.m_car*baja.g*sin(baja.theta_hill))/baja.M_effective(false);
        return du;
    };
    Eigen::Vector2d u = runge_kutta_4_step<double, 2>(integrand, Eigen::Vector2d(baja.x, baja.v), 0, dt);
    res.x = u(0);
    res.v = u(1);
    
    baja.omega_p = fixed_omega_p;
    baja.x = res.x;
    baja.v = res.v;
}

double objective(Eigen::Vector2d u) {
    BajaState baja = TR24_GAGED_GX9;
    BajaState baja_hill = TR24_GAGED_GX9;
    baja.N_g = u(0);
    baja_hill.N_g = u(0);
    double peak_power_rpm = u(1)*RADPS2RPM;
    
    const Eigen::Vector2d RATIO_LIMITS(0.9, 3.9);

    double t = 0;
    while(baja.x < 150*FT2M) {
        sim_step(baja, sim_dt, peak_power_rpm, 0, RATIO_LIMITS);
        t += sim_dt;
    }
    double t_hill = 0;
    while(baja_hill.x < 150*FT2M) {
        sim_step(baja_hill, sim_dt, peak_power_rpm, 30, RATIO_LIMITS);
        t_hill += sim_dt;
    }
    printf("N_g = %.12f,\tPeak Power RPM = %.2f,\taccel time = %.5f s,\thill time = %.5f s\n", u(0), u(1)*RADPS2RPM, t, t_hill);

    constexpr double accel_importance = 1.0;
    constexpr double hill_importance = 1.0;
    return t*accel_importance + t_hill*hill_importance;
}

int main() {
    Eigen::Vector2d lb(1, 2000*RPM2RADPS);
    Eigen::Vector2d ub(20, 3800*RPM2RADPS);

    Eigen::Vector2d u0(8.32, 3000*RPM2RADPS);
    
    auto S = minimize_gradient_golden<2>(objective, u0, lb, ub, 1e-6, 1e-2, 1e-1, 1000);

    printf("Solver %s after %llu iterations\n", S.converged ? "converged" : "did not converge", S.iterations);
    printf("Solution: N_g = %f, Peak Power RPM = %f\n", S.x(0), S.x(1)*RADPS2RPM);
}