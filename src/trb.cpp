#include "trb.hpp"

const BajaState TR24_GAGED_GX9 = {
    .engine_torque_curve=CH440_TORQUE_CURVE,
    .controls={
        .steering=0.0,
        .throttle=1.0,
        .brake_pedal=0.0,
        .awd_lever=0.0,
    },
    // Using best accel tune from 2024 season
    .cvt_tune = {
        .p_ramp_fn = [](double x){
            // measured from scan, best represented as piecewise function
            // so the function is continuous
            // with each piecewise section having a continous derivative
            if (x < 0) return 17.5;
            if (x < 5) return remap(x, 0, 5, 17.5, 16.475);
            if (x < 30.5) return 0.01*x*x-0.805*x+20.25; // intersects (5, 16.475) and (30.5, 5)
            if (x < 31) return remap(x, 30.5, 31, 5, 0);
            return 0.0;
        },
        .k_p = 60*LBF2N/IN2M,
        .m_fly = 0.536,
        .k_s = 20*LBF2N/IN2M,
        .kappa_s = 0.4644*LBF2N*IN2M*RAD2DEG, // lbf-in/deg * N/lbf * m/in * deg/rad
        .theta_s_0 = (0.5)*24*DEG2RAD, // each pretension hole is 24 degrees away from the previous
        .theta_helix = DEG2RAD*33,
    },
    .rpm_gov = 3800,
    .rpm_idle = 1800,
    .N_g = 8.32,
    .m_car = 500*LBF2KG,
    .r_wheel = 12.5*IN2M,
    .wheelbase = 1.5003,
    .com_x = 0.5, // estimate
    .com_height = 0.2, // estimate
    .theta_hill = 0,
    .phi = 12.5*DEG2RAD,
    .L = 9.1*IN2M,
    .L_b0 = 28*IN2M,
    .b_min = 0.4*IN2M,
    .b_max = 0.7*IN2M,
    .h_v = IN2M*0.675*7/8,
    .h = IN2M*7/8,
    .A_b = (IN2M*0.675*7/8)*(0.4*IN2M),
    .m_b = LBF2KG*0.6,
    .E_b = 13.8e3*LBF2N/(IN2M*IN2M), // Pa = psi * N/lbf * in^2/m^2, around 95 MPa
    .mu_b = 0.13, // from WVU paper
    .I_e = 1,
    .I_p = 1,
    .I_s = 1,
    .I_w = 1,
    .N_fly = 4,
    .r_p_inner = 0.75*IN2M,
    .d_p_max = 0.75*IN2M,
    .d_p_0 = (3 - 1.875)*IN2M, // 3in is the spring uncompressed length
    .r_cage = 65.23137756e-3,
    .r_shoulder = 41.87742519e-3,
    .L_arm = 33.02e-3,
    .r_roller = 6.35e-3,
    .x_ramp = 37.76893878e-3,
    .F1 = 0,
    .r_s_inner = 47.625e-3,
    .d_s_max = 18e-3,
    .d_s_0 = (3-1.825)*IN2M, // 3in is the spring uncompressed length
    .r_helix = 46.0375e-3,
    .F2 = 0,
    .tau_s = 0,
    .omega_p = 0,
    .omega_s = 0,
    .F_f = 0,
    .r_p = IN2M*1.1875,
    .d_p = 0,
    .d_r = 0,
    .theta1 = 0, // this is wrong, but the solver will correct this
    .theta2 = 0, // this is wrong, but the solver will correct this
    .r_s = IN2M*3.3125,
    .d_s = 0,

};

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

    ScalarFn<CVT_INDEX_COUNT> objective = [&](CVTState x){
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

    return minimize_gradient_descent(objective, x0);
}