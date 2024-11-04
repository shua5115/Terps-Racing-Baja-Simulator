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
            if (x < 5) return 0.001*remap(x, 0, 5, 17.5, 16.475);
            // if (x < 30.) return 0.001*(0.01*x*x-0.805*x+20.25); // intersects (5, 16.475) and (30.5, 5)
            // if (x < 31) return 0.001*remap(x, 30.5, 31, 5, 0);
            // return 0.001*remap(x, 31, 32, 0, -1);
            return 0.001*(0.01*x*x-0.805*x+20.25);
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
    .C_d = 0.7, // from CFD
    .A_front = 0.78, // from CAD
    .x = 0,
    .v = 0,
    .theta_hill = 0,
    .g = 9.81,
    .phi = 12.5*DEG2RAD,
    .L = 9.1*IN2M,
    .L_b0 = 34*IN2M,
    .b_min = 0.4*IN2M,
    .b_max = 0.7*IN2M,
    .h_v = 14.5e-3,
    .h = 14.5e-3,
    .A_b = (IN2M*0.675*7/8)*(0.4*IN2M),
    .m_b = LBF2KG*0.6,
    .E_b = 13.8e3*LBF2N/(IN2M*IN2M), // Pa = psi * N/lbf * in^2/m^2, around 95 MPa
    .mu_b = 100000.0, //0.13, // from WVU paper
    .mu_b_k = 0.5, // assumed to be half
    .I_e = 1,
    .I_p = 1,
    .I_s = 1,
    .I_w = 1,
    .F_resist = 10*LBF2N,
    .F1 = 0,
    .N_fly = 4,
    .r_p_inner = 0.75*IN2M,
    .d_p_max = 10e-3,
    .d_p_0 = (3 - 1.875)*IN2M, // 3in is the spring uncompressed length
    .r_cage = 65.23137756e-3,
    .r_shoulder = 41.87742519e-3,
    .L_arm = 33.02e-3,
    .r_roller = 6.35e-3,
    .x_ramp = 37.8775e-3,
    .r_s_inner = 47.625e-3,
    .d_s_max = 18e-3,
    .d_s_0 = (3-1.825)*IN2M, // 3in is the spring uncompressed length
    .r_helix = 46.0375e-3,
    .tau_s = 0,
    .omega_p = 1800*RPM2RADPS,
    .r_p = IN2M*1.1875,
    .d_p = 0.00,
    .theta1 = 0,
    .theta2 = PI*0.4,
    .r_s = IN2M*3.3125,
    .d_s = 0.00,
};

void BajaState::set_ratio_from_d_p(double d) {
    double tan_phi = tan(phi);
    double inv_tan_phi = 1.0/tan_phi;
    d_p = d;
    r_p = d_p*inv_tan_phi + r_p_min();
    r_s = solve_r_s(r_p, r_s_min(), r_s_min() + d_s_max*inv_tan_phi, L, L_b0, 8);
    d_s = d_s_max - (r_s - r_s_min())*tan_phi;
}

double BajaState::calc_tau_s() const {
    constexpr double rho_air = 1.225; // kg/m^3
    double F_external = m_car*g*sin(theta_hill) + 0.5*rho_air*C_d*A_front*v*v;
    double s_resist = sign(F_external);
    return (r_wheel/N_g)*(s_resist*F_resist + F_external);
}

OptResults<2> solve_flyweight_position(
    double theta1_guess, double theta2_guess,
    double (*ramp)(double x),
    double L1, double L2,
    double d_p, double x_ramp,
    double r_cage, double r_shoulder
) {
    Eigen::Vector2d x0(theta1_guess, theta2_guess);
    Eigen::Vector2d x_lb(-PI*0.5, 0);
    Eigen::Vector2d x_ub(PI*0.5, PI);

    ScalarFn<2> objective = [=](Eigen::Vector2d x){
        double theta1 = clamp(x(0), x_lb(0), x_ub(0));
        double theta2 = clamp(x(1), x_lb(1), x_ub(1));
        
        double d_r = x_ramp + d_p - L1*cos(theta1);
        double ramp_eval = d_r - L2*cos(theta2);

        double eq1 = L1*sin(theta1) + L2*sin(theta2) - r_cage + r_shoulder + ramp(ramp_eval);
        double eq2 = atan(diff_central(ramp, ramp_eval, 1e-6)) + PI*0.5 - theta2;

        return eq1*eq1 + eq2*eq2;
    };

    return minimize_gradient_golden(
        objective, x0, x_lb, x_ub,
        1e-8, 1e-6, 0.05, 500
    );
}

double solve_r_s(double r_p, double r_s_min, double r_s_max, double L, double L0, unsigned int N) {
    auto belt_err = [r_p, L, L0](double r_s){
        double alpha = 2*acos(clamp((r_s-r_p)/L, -1, 1));
        double beta = 2*PI-alpha;
        double e = r_p*alpha + r_s*beta + 2*sqrt(std::max(0.0, L*L - (r_p - r_s)*(r_p - r_s))) - L0; // Belt length constraint
        return e;
    };
    // return clamp(root_secant(belt_err, r_s_min, r_s_max, N), r_s_min, r_s_max);
    return root_secant(belt_err, r_s_min, r_s_max, N);
}

double solve_cvt_shift(const BajaState &baja, int debug) {
    // Precalculate derived constants
    double tau_e = baja.tau_e();
    double phi = baja.phi;
    double sin_phi = sin(phi);
    double cos_phi = cos(phi);
    double tan_phi = tan(phi);
    double inv_tan_phi = 1.0/tan_phi;
    double r_p_min = baja.r_p_min();
    double r_s_min = baja.r_s_min();
    double d_p_max = baja.d_p_max;
    double d_s_max = baja.d_s_max;
    double r_p_max = r_p_min + d_p_max*inv_tan_phi;
    double r_s_max = r_s_min + d_s_max*inv_tan_phi;
    double L = baja.L;
    double L0 = baja.L_b0;
    double L1 = baja.L_arm;
    double L2 = baja.r_roller;
    double tan_helix = tan(baja.cvt_tune.theta_helix);
    auto ramp = baja.cvt_tune.p_ramp_fn;
    double theta1 = baja.theta1;
    double theta2 = baja.theta2;

    auto objective = [&](double d_p){
        d_p = clamp(d_p, 0, d_p_max);
        double r_p = d_p*inv_tan_phi + r_p_min;
        // Constrain r_s with belt length relation independently of this numerical solution.
        double r_s = solve_r_s(r_p, r_s_min, r_s_max, L, L0, 8);

        double d_s = d_s_max - (r_s - r_s_min)*tan_phi;

        double r_fly = baja.r_shoulder + baja.L_arm*sin(theta1);
        // CVT system values
        double alpha = 2*acos(clamp((r_s-r_p)/baja.L, -1, 1));
        double beta = 2*PI-alpha;
        
        double theta_s = d_s/(baja.r_helix*tan_helix);

        // Primary subsystem forces
        double F_sp = baja.cvt_tune.k_p*(baja.d_p_0 + d_p);
        double F_flyarm = (baja.cvt_tune.m_fly*(r_fly)*baja.omega_p*baja.omega_p*L1*cos(theta1)*cos(theta2))
            /(L1*sin(theta1 + theta2) + L2*sin(2*theta2));
        double F1 = std::abs(baja.F1);

        // Secondary subsystem forces
        double F_ss = baja.cvt_tune.k_s*(baja.d_s_0 + d_s);
        double tau_ss = baja.cvt_tune.kappa_s*(baja.cvt_tune.theta_s_0 + theta_s);
        double F_helix = (baja.tau_s*0.5 + tau_ss)/(baja.r_helix*tan_helix);
        
        double eq1 = (F_sp - F_flyarm) + (alpha/beta)*(F_ss + F_helix);

        // F1 always resists shifting, bringing forces closer to equilibrium, or reaching equilibrium if close enough
        if (eq1 > F1) {
            eq1 -= F1;
        } else if(eq1 < -F1) {
            eq1 += F1;
        } else {
            eq1 = 0;
        }

        if (debug > 1){
            // printf("d_p=%f, d_s=%f, ratio=%f, f=%f\t wrap=%f, F_sp=%f, F_fly=%f, F_ss=%f, F_helix=%f\n",
            printf("%f, %f, %f, %f, %f, %f, %f, %f\n",
                d_p, d_s, r_s/r_p, eq1, F_sp, F_flyarm, F_ss, F_helix);
        }

        return eq1;
    };

    // Sample range of objective function looking for a sign change
    //printf("Searching for sign change:\n");

    int N = 100; // precision to 0.1mm
    double eq_prev = objective(0);
    double d_p_prev = 0;
    double d_p_cur = 0;
    double min_d_p = 0;
    double min_eq = eq_prev;
    bool sign_change = false;
    for(int i = 0; i < N; i++) {
        double d_p = remap(i, 0, N-1, 0, d_p_max);
        double eq = objective(d_p);
        if (abs(eq) <= abs(min_eq)) {
            min_eq = eq;
            min_d_p = d_p;
        }
        d_p_cur = d_p;
        if (sign(eq) != sign(eq_prev)) {
            sign_change = true;
            break;
        }
        eq_prev = eq;
        d_p_prev = d_p_cur;
    }
    double S = min_d_p;
    // printf("Sign change %s: [%f, %f], min_d_p=%f\n", sign_change ? "found!" : "not found...", d_p_prev, d_p_cur, min_d_p);
    if (sign_change) {
        S = root_bisection(objective, d_p_prev, d_p_cur, 10);
    }

    if(debug > 0) {
        double d_p = S;
        double r_p = d_p*inv_tan_phi + r_p_min;
        double r_s = solve_r_s(r_p, r_s_min, r_s_max, L, baja.L_b0, 8);
        double d_s = d_s_max - (r_s - r_s_min)*tan_phi;
        double F_sp = baja.cvt_tune.k_p*(baja.d_p_0 + d_p);
        double F_flyarm = (baja.cvt_tune.m_fly*(baja.r_shoulder + baja.L_arm*sin(theta1))*baja.omega_p*baja.omega_p*L1*cos(theta1)*cos(theta2))
            /(L1*sin(theta1 + theta2) + L2*sin(2*theta2));
        // Secondary subsystem forces
        double F_ss = baja.cvt_tune.k_s*(baja.d_s_0 + d_s);
        double tau_ss = baja.cvt_tune.kappa_s*(baja.cvt_tune.theta_s_0 + d_s/(baja.r_helix*tan_helix));
        double F_helix = (baja.tau_s*0.5 + tau_ss)/(baja.r_helix*tan_helix);
        double alpha = 2*acos(clamp((r_s-r_p)/baja.L, -1, 1));
        double beta = 2*PI-alpha;
        double eq1 = beta*(F_sp - F_flyarm) + alpha*(F_ss + F_helix);
        printf("%f, %f, %f, %f, %f, %f, %f, %f, %f, %f\n", baja.omega_p, baja.tau_s, d_p, d_s, r_s/r_p, F_sp, F_flyarm, F_ss, F_helix, eq1);
    }
    return S;
}

BajaDynamicsResult solve_dynamics(const BajaState &baja, double dt) {
    BajaDynamicsResult res;
    double N_p = (baja.F_flyarm() - baja.F_sp())/cos(baja.phi);
    res.F_f = (baja.tau_e()/baja.r_p + baja.tau_s/baja.r_s);
    res.slipping = abs(res.F_f) > N_p*baja.mu_b;
    res.F_f = res.slipping ? (sign(res.F_f)*std::max(0.0, N_p*baja.mu_b_k)) : res.F_f; // ensuring sign is preserved when constraining value
    res.a = (res.F_f*baja.r_s*baja.N_g/(baja.r_wheel) - baja.m_car*baja.g*sin(baja.theta_hill))
        /(baja.m_car + (baja.I_e + baja.I_p)*((baja.N_g*baja.r_s)*(baja.N_g*baja.r_s)/(baja.r_wheel*baja.r_p*baja.r_wheel*baja.r_p)) +
            baja.I_s*(baja.N_g*baja.N_g/(baja.r_wheel*baja.r_wheel)) + baja.I_w*(1.0/(baja.r_wheel*baja.r_wheel)));
    // To integrate velocity and position, we can use this definition: u = [x, x'], where u' = [x', x'']
    auto integrand = [&res](Eigen::Vector2d u, double t){
        auto du = Eigen::Vector2d();
        du(0) = u(1);
        du(1) = res.a;
        return du;
    };
    Eigen::Vector2d u = runge_kutta_4_step<double, 2>(integrand, Eigen::Vector2d(baja.x, baja.v), 0, dt);
    res.x = u(0);
    res.v = u(1);

    if (!res.slipping) {
        // Correct angular velocity of primary and secondary to enforce pulley ratio
        double omega_s = res.v * baja.N_g / baja.r_wheel;
        double H = (baja.I_p + baja.I_e)*baja.omega_p + baja.I_s*omega_s + baja.I_w*(omega_s/baja.N_g) + baja.m_car*(omega_s/baja.N_g)*baja.r_wheel*baja.r_wheel;
        double omega_p_adj = H / (baja.I_p + baja.I_e + baja.I_s*baja.r_p/baja.r_s + baja.I_w*baja.r_p/(baja.r_s*baja.N_g) + baja.m_car*baja.r_p*baja.r_wheel*baja.r_wheel/(baja.r_s*baja.N_g));
        res.omega_p = omega_p_adj;
        res.v = omega_p_adj * (baja.r_p/baja.r_s) * (baja.r_wheel / baja.N_g);
    } else {
        // Integrate omega_p by treating the primary as a separate rigid body
        res.omega_p = baja.omega_p + dt*(baja.tau_e() - res.F_f*baja.r_p)/(baja.I_e + baja.I_p); // Explicit euler b/c this acceleration is constant
    }
    
    return res;
}

void apply_dynamics_result(BajaState &state, const BajaDynamicsResult &dynamics) {
    state.x = dynamics.x;
    state.v = dynamics.v;
    state.omega_p = dynamics.omega_p;
}

BajaDynamicsResult trb_sim_step(BajaState &baja, double dt) {
    baja.tau_s = baja.calc_tau_s();
    auto S_fly = solve_flyweight_position(baja.theta1, baja.theta2, baja.cvt_tune.p_ramp_fn, baja.L_arm, baja.r_roller, baja.d_p, baja.x_ramp, baja.r_cage, baja.r_shoulder);
    baja.theta1 = S_fly.x(0);
    baja.theta2 = S_fly.x(1);
    auto d_p = solve_cvt_shift(baja);
    baja.set_ratio_from_d_p(d_p);
    auto res = solve_dynamics(baja, dt);
    apply_dynamics_result(baja, res);
    return res;
}