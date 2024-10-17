#pragma once

#include <cmath>
#include "Eigen/Dense"
#include "util.hpp"
#include "opt.hpp"

struct VehicleControls {
    double steering = 0;    // -1 to 1
    double throttle = 0;    // 0 to 1
    double brake_pedal = 0; // 0 to 1
    double awd_lever = 0;   // 0 to 1
};

struct CVTTune {
    // Primary
    double (*p_ramp_fn)(double x); // Function for ramp height vs. x, where x is 0 at the beginning of shift, all units in meters.
    double k_p; // N/m, linear spring constant of primary spring
    double m_fly; // kg, mass of all flyweights in the primary
    // Secondary
    double k_s; // N/m, linear spring constant of secondary spring
    double kappa_s; // N-m/rad, torsional spring constant of secondary spring
    double theta_s_0; // rad, angular pretension of secondary spring
    double theta_helix; // rad, angle of secondary helix ramp
};

struct BajaState {
    // Vehicle
    Eigen::MatrixX2d engine_torque_curve; // columns are: RPM, torque (N-m)
    VehicleControls controls;
    CVTTune cvt_tune;
    double rpm_gov;     // rpm, max rpm of engine
    double rpm_idle;    // rpm, min rpm of engine
    double N_g;         // Fixed gear ratio between CVT secondary and rear wheels
    double m_car;       // kg, total mass of the car (including the driver)
    double r_wheel;     // m, radius of the tires
    double wheelbase;   // m, distance between front and rear wheels
    double com_x;       // m, forward distance from rear wheel to center of mass
    double com_height;  // m, distance from ground to center of mass 
    // Powertrain
    double phi;         // rad, half of angle between sheaves
    double L;           // m, distance between primary and secondary
    double L_b0;        // m, unstretched belt length
    double b_min;       // m, minimum width of belt
    double b_max;       // m, maximum width of belt
    double h_v;         // m, height of v-shaped part of belt
    double h;           // m, total height of belt
    double A_b;         // m^2, belt cross section area, but only the part that stretches
    double m_b;         // kg, mass of belt
    double E_b;         // Pa, Young's modulus of belt
    double mu_b;        // Belt coefficient of friction with sheaves
    double I_e;         // kg-m^2, total moment of inertia of spinning engine components
    double I_p;         // kg-m^2, total moment of inertia of primary components
    double I_s;         // kg-m^2, total moment of inertia of secondary components
    double I_w;         // kg-m^2, total moment of inertia of all four wheels
    // Primary
    uint8_t N_fly;      // Number of flyweight linkages in primary
    double r_p_inner;   // m, radius where bottom of primary sheaves touch
    double d_p_max;     // m, max linear gap between primary sheaves
    double d_p_0;       // m, primary spring initial displacement
    double r_cage;      // m, radius from primary axis to outer edge of primary ramp
    double r_shoulder;  // m, radius from primary axis to flyweight arm pivot
    double L_arm;       // m, length of flyweight arm
    double r_roller;    // m, radius of rollers in primary
    double x_ramp;      // m, offset from flyweight arm pivot to furthest outwards edge of ramp
    // Secondary
    double r_s_inner;   // m, radius where bottom of secondary sheaves touch
    double d_s_max;     // m, max linear gap between secondary sheaves
    double d_s_0;       // m, Secondary spring initial displacement
    double r_helix;     // m, Radius of secondary helix ramp

    // Time-Dependent

    double tau_s;       // N-m, torque applied to secondary from gearbox
    double tau_p;       // N-m, torque applied to primary from engine
    double omega_p;     // rad/s, angular velocity of primary
    double omega_s;     // rad/s, angular velocity of secondary
    double T0;          // N, slack side belt tension
    double T1;          // N, taut side belt tension
    double L_b;         // m, current belt length
    // Primary
    double d_p;         // m, linear displacement of primary sheave during shift
    double d_r;         // m, linear displacement of roller from outermost edge of ramp
    double theta1;      // rad, angle beween flyweight arm and primary axis
    double theta2;      // rad, angle between primary ramp surface normal and primary axis
    // Secondary
    double d_s;         // rad, linear displacement of secondary sheave during shift

    // Derived Constants

    // Average width of belt
    double b() {
        return (b_min + b_max)/2; 
    }

    // Minimum radius for belt on primary sheave
    double r_p_min() {
        return (std::min(d_p_max, b_min)*0.5)/(tan(phi)) + r_p_inner + h_v*0.5;
    }

    // Minimum radius for belt on secondary sheave
    double r_s_min() {
        return (std::min(d_s_max, b_min)*0.5)/(tan(phi)) + r_s_inner + h_v*0.5;
    }

    double rho_b() {
        return m_b/(A_b*L_b0);
    }

    double theta_s_max() {
        return d_s_max/(r_helix*tan(cvt_tune.theta_helix));
    }

    // Derived Variables

    double r_p() {
        return clamp(d_p, 0, d_p_max)/tan(phi) + r_p_min();
    }

    double r_s() {
        return clamp(d_s_max - d_s, 0, d_s_max)/tan(phi) + r_s_min();
    }

    double alpha() {
        return 2*acos((r_s()-r_p())/L);
    }

    double beta() {
        return 2*acos((r_p()-r_s())/L);
    }

    double theta_s() {
        return (d_s)/(r_helix*tan(cvt_tune.theta_helix));
    }

    double r_fly() {
        return r_shoulder + L_arm*sin(theta1);
    }

    // Misc

    double throttle_scale(double torque, double u_gas) {
        return torque*(sin(u_gas*PI*0.5)*(1 - rpm_idle/rpm_gov) + rpm_idle/rpm_gov);
    }
};

OptResults<3> solve_flyweight_position(
    double theta1_guess, double theta2_guess, double d_r_guess,
    double (*ramp)(double x),
    double L1, double L2,
    double d_p, double x_ramp,
    double r_cage, double r_shoulder
);

constexpr size_t BAJA_SIZE = sizeof(BajaState);

// Column units: RPM, N-m
const Eigen::MatrixX2d CH440_TORQUE_CURVE = Eigen::MatrixX2d({
    {0   , 25.4},
    {1800, 25.4},
    {2000, 25.8},
    {2200, 26.4},
    {2400, 26.8},
    {2600, 26.7},
    {2800, 26.6},
    {3000, 26.4},
    {3200, 25.7},
    {3400, 24.8},
    {3600, 23.9},
    {3800, 23.9}, // Governor prevents torque increase above 3800, need to verify behavior
});

// const BajaState TR24_STATE = BajaState{
//     .engine_torque_curve=CH440_TORQUE_CURVE,
//     .N_g=8.32,
//     .m_car=500*LBF2KG,
//     .r_wheel=11.5*IN2M,
// };

/*
constexpr CVTConfig GAGED_GX9_CONFIG = {
    .r_p = 0.975*IN2M,          // m, from WVU
    .R_p = 2.775*IN2M,          // m, from WVU
    .r_s = 2.498*IN2M,          // m, from WVU
    .R_s = 3.803*IN2M,          // m, from WVU
    .v_angle = 23*DEG2RAD,      // m, from WVU
    .g_p = 0.8*IN2M,          // m, from WVU
    .g_s = 0.8*IN2M,          // m, from WVU
    .L = 9.1*IN2M,              // m, as designed
    .l_b = 28*IN2M,             // m, measured
    .A_b = 0.41*IN2M*IN2M,      // m^2, from WVU
    .m_b = 0.249,               // kg, measured
    .E_b = 95147650.646,        // Pa, =13.8 ksi, from WVU
    .G_b = 35025367.049,        // Pa, =5.08 ksi, from WVU
    .mu_b = 0.13,               // from WVU
    .p_spring_x0 = 0.975*IN2M,  // m, from WVU
    .s_spring_x0 = 1.44*IN2M,   // m, from WVU
    .r_helix = 1.455*IN2M, // m, from WVU
    .m_fly_arm = 0.058375612,// kg, from WVU
    // .r_fly_min = 1.875*IN2M,    // m, from WVU
    // .r_fly_max = 2.345*IN2M,    // m, from WVU
};

double gaged_gx9_primary_ramp_curved(double x) {
    // TODO, use function fit from excel
    return 0.0;
}

double gaged_gx9_primary_ramp_flat(double x) {
    // TODO, use function fit from excel
    return 0.0;
}
*/