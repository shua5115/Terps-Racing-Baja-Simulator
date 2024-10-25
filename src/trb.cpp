#include "trb.hpp"

const BajaState TR24_GAGED_GX9 = {
    .engine_torque_curve=CH440_TORQUE_CURVE,
    .controls={
        .steering=0.0,
        .throttle=0.0,
        .brake_pedal=0.0,
        .awd_lever=0.0,
    },
    // Using best accel tune from 2024 season
    .cvt_tune = {
        .p_ramp_fn = [](double x){
            // measured from scan, best represented as piecewise function
            // so the function is continuous
            // with each piecewise section having a continous derivative
            // the ramp must never have a derivative of zero
            x *= 1000; // convert from m to mm
            // if (x < 0) return 17.5;
            if (x < 5) return remap(x, 0, 5, 17.5, 16.475);
            if (x < 30.) return (0.01*x*x-0.805*x+20.25); // intersects (5, 16.475) and (30.5, 5)
            if (x < 31) return remap(x, 30.5, 31, 5, 0);
            return remap(x, 31, 32, 0, -1);
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
    .L_b0 = 34*IN2M,
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
    .x_ramp = 37.8775e-3,
    .F1 = 0,
    .r_s_inner = 47.625e-3,
    .d_s_max = 18e-3,
    .d_s_0 = (3-1.825)*IN2M, // 3in is the spring uncompressed length
    .r_helix = 46.0375e-3,
    .F2 = 0,
    .tau_s = 0,
    .omega_p = 1800*RPM2RADPS,
    .omega_s = 0,
    .F_f = 2000, // Initial guess is about 4000 N
    .r_p = IN2M*1.1875,
    .d_p = 0.00,
    .d_r = 0,
    .theta1 = 0, // this is wrong, but the solver will correct this
    .theta2 = PI*0.5, // this is wrong, but the solver will correct this
    .r_s = IN2M*3.3125,
    .d_s = 0.00,
};

OptResults<3> solve_flyweight_position(
    double theta1_guess, double theta2_guess, double d_r_guess,
    double (*ramp)(double x),
    double L1, double L2,
    double d_p, double x_ramp,
    double r_cage, double r_shoulder
) {
    Eigen::Vector3d x0(theta1_guess, theta2_guess, d_r_guess);
    Eigen::Vector3d x_lb(-PI*0.5, -PI*0.5, 0);
    Eigen::Vector3d x_ub(PI*0.5, PI*0.5, INFINITY);

    ScalarFn<3> objective = [=](Eigen::Vector3d x){
        double theta1 = clamp(x(0), x_lb(0), x_ub(0));
        double theta2 = clamp(x(1), x_lb(1), x_ub(1));
        double d_r = clamp(x(2), x_lb(2), x_ub(2));
        // double theta1 = x(0);
        // double theta2 = x(1);
        // double d_r = x(2);

        double ramp_eval = d_r - L2*cos(theta2);

        double eq1 = L1*cos(theta1) - x_ramp - d_p + d_r;
        double eq2 = L1*sin(theta1) + L2*sin(theta2) - r_cage + r_shoulder + ramp(ramp_eval);
        double eq3 = atan(diff_central(ramp, ramp_eval, 1e-6)) + PI*0.5 - theta2;

        return 0.1*eq1*eq1 + 0.1*eq2*eq2 + 100*eq3*eq3;
    };

    return minimize_gradient_golden(
        objective, x0, x_lb, x_ub,
        1e-5, 1e-6, 0.05, 1000
    );
}

OptResults<CVT_INDEX_COUNT> solve_cvt_shift(const BajaState &baja) {
    // Step 1. solve assuming no slip
    // Step 2. if no slip condition is invalid, then solve with slipping condition
    constexpr double F_f_scale = 1e-4;
    constexpr double F_f_unscale = 1e4;
    // Initial guess
    CVTState x0;
    x0(R_P) = baja.r_p;
    x0(R_S) = baja.r_s;
    x0(F_F) = baja.F_f*F_f_scale; // convert to similar order of magnitude to radius

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

    // 5 MN is about 1100 lbf
    Eigen::Vector3d x_lb(r_p_min, r_s_min, 0);
    Eigen::Vector3d x_ub(r_p_min + baja.d_p_max/tan_phi, r_s_min + baja.d_s_max/tan_phi, INFINITY);

    ScalarFn<CVT_INDEX_COUNT> objective = [&](CVTState x){
        double r_p = x(R_P); //clamp(x(R_P), x_lb(R_P), x_ub(R_P));
        double r_s = x(R_S); //clamp(x(R_S), x_lb(R_S), x_ub(R_S));
        double F_f = x(F_F)*F_f_unscale; //clamp(x(F_F)*F_f_unscale, x_lb(F_F), x_ub(F_F)); // convert from mega newtons to newtons

        // CVT system values
        double d_p = (r_p-r_p_min)*tan_phi;
        double d_s = baja.d_s_max - (r_s-r_s_min)*tan_phi;
        double alpha = 2*acos(clamp((r_s-r_p)/baja.L, -1, 1));
        double beta = 2*acos(clamp((r_p-r_s)/baja.L, -1, 1));
        double theta_s = d_s/(baja.r_helix*tan(baja.cvt_tune.theta_helix));
        double L_b = r_p*alpha + r_s*beta + 2*sqrt(std::max(0.0, baja.L*baja.L - (r_p - r_s)*(r_p - r_s)));
        double T0 = baja.E_b*baja.A_b*(L_b/baja.L_b0 - 1) - rho_b*baja.A_b*(r_s*r_s*baja.omega_s*baja.omega_s);
        double N_b = ((F_f == 0) ? T0 : (abs(F_f)/log(abs(F_f) + 1) - (T0 - 1))); // normal force from belt on either sheave
        double N_p = N_b*alpha/sin_phi;

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
        double F_bs = N_b*beta/tan_phi;
        double F_helix = baja.cvt_tune.kappa_s*(baja.cvt_tune.theta_s_0 + d_s/(baja.r_helix*tan_helix))/(baja.r_helix*tan_helix);

        // Check if slipping
        bool no_slip = tau_e/r_p <= N_p*baja.mu_b;
        double tau_p = no_slip ? (tau_e) : (N_p*baja.mu_b*r_p);

        // Three variables, so three equations are needed:
        double eq1 = tau_p/r_p + baja.tau_s/r_s - F_f; // Sum of forces in belt section
        double eq2 = F_sp + F_bp - 4*F_flyarm; // Sum of forces in primary
        double eq3 = F_ss - F_bs + F_helix; // Sum of forces in secondary
        // F1 and F2 always resist shifting, bringing forces closer to equilibrium, or reaching equilibrium if close enough
        if (eq2 > F1) {
            eq2 -= F1;
        } else if(eq2 < -F1) {
            eq2 += F1;
        } else {
            eq2 = 0;
        }
        if (eq3 > F2) {
            eq3 -= F2;
        } else if(eq3 < -F2) {
            eq3 += F2;
        } else {
            eq3 = 0;
        }

        // printf("\tx=[%f, %f, %f], eq=[%f, %f, %f]\n", r_p, r_s, F_f, eq1, eq2, eq3);

        return abs(eq1) + abs(eq2) + abs(eq3);
    };
    auto S = minimize_gradient_golden(objective, x0,
        x_lb, x_ub,
        1e-8, 1e-8, 1e-3, 1000
    );
    S.x(F_F) *= F_f_unscale; // convert from mega newtons to newtons
    return S;
}

double solve_r_s(double r_p, double r_s_min, double r_s_max, double L, double L0, unsigned int N) {
    auto belt_err = [r_p, L, L0](double r_s){
        double alpha = 2*acos(clamp((r_s-r_p)/L, -1, 1));
        double beta = 2*PI-alpha;
        double e = r_p*alpha + r_s*beta + 2*sqrt(L*L - (r_p - r_s)*(r_p - r_s)) - L0; // Belt length constraint
        return e;
    };
    return root_secant(belt_err, r_s_min, r_s_max, N);
}

OptResults<4> solve_cvt_shift_no_slip(const BajaState &baja) {
    Eigen::Vector<double, 4> x0(
        baja.d_p,
        baja.theta1,
        baja.theta2,
        baja.d_r
    );

    // Precalculate derived constants
    double tau_e = matrix_linear_lookup(baja.engine_torque_curve, std::max(RADPS2RPM*baja.omega_p, baja.rpm_idle));
    tau_e = baja.throttle_scale(tau_e, baja.controls.throttle);
    double phi = baja.phi;
    double tan_phi = tan(phi);
    double inv_tan_phi = 1.0/tan_phi;
    double r_p_min = baja.r_p_min();
    double r_s_min = baja.r_s_min();
    double d_p_max = baja.d_p_max;
    double d_s_max = baja.d_s_max;
    double r_s_max = d_s_max*inv_tan_phi + r_s_min;
    // double rho_b = baja.rho_b();
    // double sin_phi = sin(phi);
    // double cos_phi = cos(phi);
    double tan_helix = tan(baja.cvt_tune.theta_helix);
    double L = baja.L;
    double L0 = baja.L_b0;
    double L1 = baja.L_arm;
    double L2 = baja.r_roller;
    auto ramp = baja.cvt_tune.p_ramp_fn;

    Eigen::Vector<double, 4> x_lb(0, -PI*0.5, -PI*0.5, 0);
    Eigen::Vector<double, 4> x_ub(d_p_max, PI*0.5, PI*0.5, INFINITY);

    // Eigen::Vector2d x_lb(0, 0);
    // Eigen::Vector2d x_ub(INFINITY, baja.L);

    std::cout << "Lower bound: " << Eigen::Transpose(x_lb) << ", Upper bound: " << Eigen::Transpose(x_ub) << "\n";

    bool print_final = false;
    ScalarFn<4> objective = [&](Eigen::Vector<double, 4> x){
        double d_p = clamp(x(0), x_lb(0), x_ub(0));
        double theta1 = clamp(x(1), x_lb(1), x_ub(1));
        double theta2 = clamp(x(2), x_lb(2), x_ub(2));
        double d_r = clamp(x(3), x_lb(3), x_ub(3));
        double r_p = d_p*inv_tan_phi + r_p_min;

        // Constrain r_s with belt length relation independently of this numerical solution.
        // This makes the solution more stable.
        double r_s = solve_r_s(r_p, r_s_min, r_s_max, L, L0, 30);

        double d_s = d_s_max - (r_s - r_s_min)*tan_phi;

        double r_fly = baja.r_shoulder + baja.L_arm*sin(theta1);
        // CVT system values
        double alpha = 2*acos(clamp((r_s-r_p)/baja.L, -1, 1));
        double beta = 2*PI-alpha; //2*acos(clamp((r_p-r_s)/baja.L, -1, 1));
        
        double theta_s = d_s/(baja.r_helix*tan_helix);

        double ramp_eval = d_r - L2*cos(theta2);

        // Primary subsystem forces
        double F_sp = baja.cvt_tune.k_p*(baja.d_p_0 + d_p);
        // whooo boy
        double F_flyarm = (0.25*baja.cvt_tune.m_fly*(baja.r_shoulder + baja.L_arm*sin(theta1))*baja.omega_p*baja.omega_p*baja.L_arm*cos(theta1)*cos(theta2))
            /(baja.L_arm*sin(theta1 + theta2) + baja.r_roller*sin(2*theta2));
        double F1 = std::abs(baja.F1);

        // Secondary subsystem forces
        double F_ss = baja.cvt_tune.k_s*(baja.d_s_0 + d_s);
        //N = N-m/rad * rad * 1/m
        double F_helix = (baja.cvt_tune.kappa_s*(baja.cvt_tune.theta_s_0 + theta_s) + baja.tau_s*0.5)*baja.r_helix/(tan_helix);

        // Assume no slip
        double tau_p = tau_e;

        // double eq1 = r_p*alpha + r_s*beta + 2*sqrt(baja.L*baja.L - (r_p - r_s)*(r_p - r_s)) - baja.L_b0; // Belt length constraint
        double eq2 = (F_ss + F_helix)*beta/alpha - 4*F_flyarm + F_sp; // Sum of forces between primary and secondary
        double eq3 = L1*cos(theta1) - baja.x_ramp - d_p + d_r;
        double eq4 = L1*sin(theta1) + L2*sin(theta2) - baja.r_cage + baja.r_shoulder + ramp(ramp_eval);
        double eq5 = atan(diff_central(ramp, ramp_eval, 1e-6)) + PI*0.5 - theta2;

        if (print_final) {
            printf("eqs: [%f, %f, %f, %f]\n", eq2, eq3, eq4, eq5);
        }

        // F1 always resists shifting, bringing forces closer to equilibrium, or reaching equilibrium if close enough
        if (eq2 > F1) {
            eq2 -= F1;
        } else if(eq2 < -F1) {
            eq2 += F1;
        } else {
            eq2 = 0;
        }

        return eq2*eq2 + eq3*eq3 + eq4*eq4 + eq5*eq5;
    };
    auto S = minimize_gradient_golden(objective, x0,
        x_lb, x_ub,
        1e-9, 1e-9, 0.005, 100000
    );

    print_final = true;
    objective(S.x);

    return S;
}