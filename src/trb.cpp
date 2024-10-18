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

OptResults<CVT_INDEX_COUNT> solve_cvt_shift(const BajaState &baja) {
    // Step 1. solve assuming no slip
    // Step 2. if no slip condition is invalid, then solve with slipping condition
    
    // Initial guess
    CVTState x0;
    x0(R_P) = baja.r_p;
    x0(R_S) = baja.r_s;
    x0(F_F) = baja.F_f;

    // Precalculate derived constants
    double tau_e = matrix_linear_lookup(baja.engine_torque_curve, std::max(RADPS2RPM*baja.omega_p, baja.rpm_idle));
    tau_e = baja.throttle_scale(tau_e, baja.controls.throttle);
    double phi = baja.phi;
    double r_p_min = baja.r_p_min();
    double r_s_min = baja.r_s_min();
    double r_fly = baja.r_shoulder + baja.L_arm*sin(baja.theta1);
    double rho_b = baja.rho_b();
    double sin_phi = sin(phi);
    double cos_phi = cos(phi);
    double tan_phi = tan(phi);
    double tan_helix = tan(baja.cvt_tune.theta_helix);

    ScalarFn<8> objective = [&](CVTState x){
        double r_p = x(R_P);
        double r_s = x(R_S);
        double F_f = x(F_F);

        // CVT system values
        double d_p = (r_p-r_p_min)*tan_phi;
        double d_s = baja.d_s_max - (r_s-r_s_min)*tan_phi;
        double alpha = 2*acos((r_s-r_p)/baja.L);
        double beta = 2*acos((r_p-r_s)/baja.L);
        double theta_s = d_s/(baja.r_helix*tan(baja.cvt_tune.theta_helix));
        double L_b = r_p*alpha + r_s*beta + 2*sqrt(baja.L*baja.L - (r_p - r_s)*(r_p - r_s));
        double T0 = baja.E_b*baja.A_b*(L_b/baja.L_b0 - 1) - rho_b*baja.A_b*(r_s*r_s*baja.omega_s*baja.omega_s);
        double T1 = F_f + T0;
        double N_p = (F_f/log(F_f + 1) - (T0 - 1))*alpha/sin_phi;

        // Primary subsystem forces
        double F_sp = baja.cvt_tune.k_p*(baja.d_p_0 + d_p);
        double F_bp = N_p*cos_phi; // equivalent to full equation
        // whooo boy
        double F_flyarm = (0.25*baja.cvt_tune.m_fly*(baja.r_shoulder + baja.L_arm*sin(baja.theta1))*baja.omega_p*baja.omega_p*baja.L_arm*cos(baja.theta1)*cos(baja.theta2))
            /(baja.L_arm*sin(baja.theta1 + baja.theta2) + baja.r_roller*sin(2*baja.theta2));
        double F1 = std::abs(baja.F1);
        double F2 = std::abs(baja.F2);

        // Secondary subsystem forces
        double F_ss = baja.cvt_tune.k_s*(baja.d_s_0 + d_s);
        double F_bs = (F_f/log(F_f + 1) - (T0 - 1))*beta/tan_phi;
        double F_helix = baja.cvt_tune.kappa_s*(baja.cvt_tune.theta_s_0 + d_s/(baja.r_helix*tan_helix))/(baja.r_helix*tan_helix);

        // Check if slipping
        bool no_slip = tau_e/r_p <= N_p*baja.mu_b;
        double tau_p = no_slip ? (tau_e/r_p) : (N_p*baja.mu_b);

        // Three variables, so three equations are needed:
        double eq1 = tau_p/r_p + baja.tau_s/r_s - F_f; // Sum of forces between sheaves
        double eq2 = F_sp + F_bp - 4*F_flyarm; // Sum of forces in primary
        // F1 always resists shifting, bringing forces closer to equilibrium, or reaching equilibrium if close enough
        if (eq2 > F1) {
            eq2 -= F1;
        } else if(eq2 < -F1) {
            eq2 += F1;
        } else {
            eq2 = 0;
        }
        double eq3 = F_ss - F_bs + F_helix - F2; // Sum of forces in secondary
        if (eq3 > F2) {
            eq3 -= F2;
        } else if(eq3 < -F2) {
            eq3 += F2;
        } else {
            eq3 = 0;
        }

        return eq1*eq1 + eq2*eq2 + eq3*eq3;
    };
}