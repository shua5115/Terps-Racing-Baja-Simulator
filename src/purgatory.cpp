// This is what I call "text purgatory", text that goes unused but may be useful in the future

#include "trb.hpp"

// void cvt_shift_solver_belt_stretch() {
//     BajaState state = TR24_GAGED_GX9;
//     state.omega_p = 3000*RPM2RADPS;
//     state.tau_s = 0*18.5*LBF2N/FT2M; // 18.5 lb-ft, max output from engine
//     state.F_f = 0;
//     state.controls.throttle = 1;

//     auto start = std::chrono::high_resolution_clock::now();
//     auto S = solve_cvt_shift_unstable(state);
//     auto finish = std::chrono::high_resolution_clock::now();
//     auto dur = finish - start;
//     printf("Shift solver finished in %f us\n", dur.count()*0.001);
//     printf("Solver %s after %llu iterations\n", S.converged ? "converged" : "did not converge", S.iterations);
//     printf("Error: %e = %e^2\n", S.f_of_x, sqrt(S.f_of_x));
//     double d_p = S.x(0);
//     double d_s = S.x(1);
//     double r_p = d_p/tan(state.phi) + state.r_p_min();
//     double r_s = d_s/tan(state.phi) + state.r_s_min();
//     printf("Solution: d_p=%f m, d_s=%f m, r_p=%f, r_s=%f, ratio=%f\n", d_p, d_s, r_p, r_s, r_s/r_p);
//     ASSERT(S.converged && (S.f_of_x < 1e-3));
// }

// void cvt_shift_solver_no_stretch() {
//     BajaState state = TR24_GAGED_GX9;
//     state.omega_p = 1800*RPM2RADPS;
//     state.tau_s = 0*18.5*LBF2N/FT2M; // 18.5 lb-ft, max output from engine
//     state.controls.throttle = 1;

//     auto start = std::chrono::high_resolution_clock::now();
//     auto S = solve_cvt_shift_no_stretch(state);
//     auto finish = std::chrono::high_resolution_clock::now();
//     auto dur = finish - start;
//     printf("Shift solver finished in %f us\n", dur.count()*0.001);
//     printf("Solver %s after %llu iterations\n", S.converged ? "converged" : "did not converge", S.iterations);
//     printf("Error: %e = %e^2\n", S.f_of_x, sqrt(S.f_of_x));
//     double d_p = S.x(0);
//     double r_p = d_p/tan(state.phi) + state.r_p_min();
//     double r_s = solve_r_s(r_p, state.r_s_min(), state.r_s_min() + state.d_s_max/tan(state.phi), state.L, state.L_b0, 15);
//     double d_s = state.d_s_max - (r_s - state.r_s_min())*tan(state.phi);
//     printf("Solution: d_p=%f m, d_s=%f m, r_p=%f, r_s=%f, ratio=%f\n", d_p, d_s, r_p, r_s, r_s/r_p);

//     double tau_e = matrix_linear_lookup(state.engine_torque_curve, std::max(RADPS2RPM*state.omega_p, state.rpm_idle));
//     tau_e = state.throttle_scale(tau_e, state.controls.throttle);


//     ASSERT(S.converged && (S.f_of_x < 1e-3));
// }

OptResults<6> solve_cvt_shift_unstable(const BajaState &baja) {
    constexpr double F_f_scale = 1;
    constexpr double F_f_unscale = 1;

    // Initial guess
    Eigen::Vector<double, 6> x0(
        baja.d_p,
        baja.d_s,
        baja.F_f*F_f_scale, // convert to similar order of magnitude to radius
        baja.theta1,
        baja.theta2,
        baja.d_r
    );

    // Precalculate derived constants
    double tau_e = matrix_linear_lookup(baja.engine_torque_curve, std::max(RADPS2RPM*baja.omega_p, baja.rpm_idle));
    tau_e = baja.throttle_scale(tau_e, baja.controls.throttle);
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
    double rho_b = baja.rho_b();
    double tan_helix = tan(baja.cvt_tune.theta_helix);
    auto ramp = baja.cvt_tune.p_ramp_fn;

    Eigen::Vector<double, 6> x_lb(0, 0, -INFINITY, -PI*0.5, -PI*0.5, 0);
    Eigen::Vector<double, 6> x_ub(d_p_max, d_s_max, INFINITY, PI*0.5, PI*0.5, INFINITY);

    ScalarFn<6> objective = [&](Eigen::Vector<double, 6> x){
        double d_p = clamp(x(0), x_lb(0), x_ub(0));
        double d_s = clamp(x(1), x_lb(1), x_ub(1));
        double F_f = clamp(x(2)*F_f_unscale, x_lb(2), x_ub(2));
        double theta1 = clamp(x(3), x_lb(3), x_ub(3));
        double theta2 = clamp(x(4), x_lb(4), x_ub(4));
        double d_r = clamp(x(5), x_lb(5), x_ub(5));

        // CVT system values
        double r_p = d_p*inv_tan_phi + r_p_min;
        double r_s = (d_s_max - d_s)*inv_tan_phi + r_s_min;
        double r_fly = baja.r_shoulder + L1*sin(theta1);
        double alpha = 2*acos(clamp((r_s-r_p)/L, -1, 1));
        double beta = 2*PI-alpha;
        double theta_s = d_s/(baja.r_helix*tan_helix);
        double ramp_eval = d_r - L2*cos(theta2);

        double L_b = r_p*alpha + r_s*beta + 2*sqrt(std::max(0.0, L*L - (r_p - r_s)*(r_p - r_s)));
        double T0 = baja.E_b*baja.A_b*(L_b/L0 - 1) - rho_b*baja.A_b*(r_s*r_s*baja.omega_s*baja.omega_s);
        double N_b = ((F_f == 0) ? T0 : (abs(F_f)/log(abs(F_f) + 1) - (T0 - 1))); // normal force from belt on either sheave
        double N_p = N_b*alpha/sin_phi;

        // Primary subsystem forces
        double F_sp = baja.cvt_tune.k_p*(baja.d_p_0 + d_p);
        double F_bp = N_p*cos_phi; // equivalent to full equation
        // whooo boy
        double F_flyarm = (baja.cvt_tune.m_fly*(baja.r_shoulder + L1*sin(theta1))*baja.omega_p*baja.omega_p*L1*cos(theta1)*cos(theta2))
            /(L1*sin(theta1 + theta2) + L2*sin(2*theta2));
        double F1 = std::abs(baja.F1);
        double F2 = std::abs(baja.F2);

        // Secondary subsystem forces
        double F_ss = baja.cvt_tune.k_s*(baja.d_s_0 + d_s);
        double F_bs = N_b*beta/tan_phi;
        double F_helix = (baja.cvt_tune.kappa_s*(baja.cvt_tune.theta_s_0 + theta_s) + baja.tau_s*0.5)*baja.r_helix/(tan_helix);

        // Check if slipping
        bool no_slip = tau_e/r_p <= N_p*baja.mu_b;
        double tau_p = no_slip ? (tau_e) : (N_p*baja.mu_b*r_p);

        // System of equations
        double eq1 = tau_p/r_p + baja.tau_s/r_s - F_f; // Sum of forces in belt section
        double eq2 = F_sp + F_bp - F_flyarm; // Sum of forces in primary
        double eq3 = F_ss - F_bs + F_helix; // Sum of forces in secondary
        double eq4 = L1*cos(theta1) - baja.x_ramp - d_p + d_r;
        double eq5 = L1*sin(theta1) + L2*sin(theta2) - baja.r_cage + baja.r_shoulder + ramp(ramp_eval);
        double eq6 = atan(diff_central(ramp, ramp_eval, 1e-6)) + PI*0.5 - theta2;

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

        return 100*eq1*eq1 + eq2*eq2 + eq3*eq3 + eq4*eq4 + eq5*eq5 + 100*eq6*eq6;
    };
    auto S = minimize_gradient_golden(objective, x0,
        x_lb, x_ub,
        1e-12, 1e-8, 1e-6, 100000, false
    );
    S.x(2) *= F_f_unscale; // correct F_f scale
    return S;
}

OptResults<2> solve_cvt_shift_no_stretch(const BajaState &baja) {
    Eigen::Vector<double, 2> x0(
        baja.d_p,
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

    Eigen::Vector<double, 2> x_lb(0, 0);
    Eigen::Vector<double, 2> x_ub(d_p_max, baja.x_ramp + d_p_max + 2*L1);

    std::cout << "Lower bound: " << Eigen::Transpose(x_lb) << ", Upper bound: " << Eigen::Transpose(x_ub) << "\n";

    bool print_final = false;
    ScalarFn<2> objective = [&](Eigen::Vector<double, 2> x){
        double d_p = clamp(x(0), x_lb(0), x_ub(0));
        double d_r = clamp(x(1), x_lb(1), x_ub(1));
        double r_p = d_p*inv_tan_phi + r_p_min;

        // approximating slope at ramp by sampling while ignoring diameter of roller
        // this makes theta1's and theta2's equation explicit instead of implicit
        double theta1 = acos(clamp((baja.x_ramp + d_p - d_r)/L1, -1, 1));
        double theta2 = atan(diff_central(ramp, d_r, 1e-6)) + PI*0.5;

        // Constrain r_s with belt length relation independently of this numerical solution.
        // This makes the solution more stable.
        double r_s = solve_r_s(r_p, r_s_min, r_s_max, L, L0, 15);

        double d_s = d_s_max - (r_s - r_s_min)*tan_phi;

        double r_fly = baja.r_cage - ramp(d_r) - L2; // baja.r_shoulder + baja.L_arm*sin(theta1);
        // CVT system values
        double alpha = 2*acos(clamp((r_s-r_p)/baja.L, -1, 1));
        double beta = 2*PI-alpha; //2*acos(clamp((r_p-r_s)/baja.L, -1, 1));
        
        double theta_s = d_s/(baja.r_helix*tan_helix);

        // Primary subsystem forces
        double F_sp = baja.cvt_tune.k_p*(baja.d_p_0 + d_p);
        // whooo boy
        double F_flyarm = (baja.cvt_tune.m_fly*(r_fly)*baja.omega_p*baja.omega_p*L1*cos(theta1)*cos(theta2))
            /(L1*sin(theta1 + theta2) + L2*sin(2*theta2));
        double F1 = std::abs(baja.F1);

        // Secondary subsystem forces
        double F_ss = baja.cvt_tune.k_s*(baja.d_s_0 + d_s);
        //N = ((N-m/rad * rad) + N-m) * 1/m
        double F_helix = (baja.cvt_tune.kappa_s*(baja.cvt_tune.theta_s_0 + theta_s) + baja.tau_s*0.5)/(baja.r_helix*tan_helix);

        double eq1 = (F_ss + F_helix)*beta/alpha - F_flyarm + F_sp; // Sum of forces between primary and secondary
        double eq2 = L1*sin(theta1) - baja.r_cage + baja.r_shoulder + ramp(d_r);

        if (print_final) {
            printf("eqs: [%f, %f], thetas: [%f, %f]\n", eq1, eq2, theta1, theta2);
        }

        // F1 always resists shifting, bringing forces closer to equilibrium, or reaching equilibrium if close enough
        if (eq1 > F1) {
            eq1 -= F1;
        } else if(eq1 < -F1) {
            eq1 += F1;
        } else {
            eq1 = 0;
        }

        return eq1*eq1 + eq2*eq2;
    };
    auto S = minimize_gradient_golden(objective, x0,
        x_lb, x_ub,
        1e-10, 1e-10, 0.005, 5000
    );

    print_final = true;
    objective(S.x);

    return S;
}