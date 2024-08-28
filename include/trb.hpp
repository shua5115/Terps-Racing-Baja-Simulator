#pragma once

#include <vector>
#include "Eigen/Dense"
#include "util.hpp"

namespace TRB {

struct VehicleConfig {
    // columns are: rad/s, torque (N-m)
    Eigen::MatrixX2d engine_torque_curve;
    double fixed_gear_ratio;
    double wheel_radius;

    VehicleConfig(Eigen::MatrixX2d torque_curve, double ratio) : engine_torque_curve(torque_curve), fixed_gear_ratio(ratio) { }
};

struct VehicleControls {
    double steering = 0; // -1 to 1
    double throttle = 0; // 0 to 1
    double brake_pedal = 0; // 0 to 1
    double awd_lever = 0; // 0 or 1
};

// Variables about the CVT that are constant
struct CVTConfig {
    double r_p; // m, inner radius of primary
    double R_p; // m, outer radius of primary
    double r_s; // m, inner radius of secondary
    double R_s; // m, outer radius of secondary
    double v_angle; // radians, angle between sheaves
    double p_gap; // m, max gap between primary sheaves
    double s_gap; // m, max gap between secondary sheaves
    double L; // m, distance between primary and secondary
    double l_b; // m, belt length
    double A_b; // m^2, belt cross section area
    double m_b; // kg, mass of belt
    double E_b; // Pa, Young's modulus of belt
    double G_b; // Pa, shear modulus of belt
    double mu_b; // Belt coefficient of friction
    double p_spring_x0; // m, Primary spring initial displacement
    double s_spring_x0; // m, Secondary spring initial displacement
    double helix_radius; // m, Radius of secondary helix ramp
    double fly_arm_mass; // kg, Mass of arm connected to flyweight
    double fly_r_min;    // kg, Minimum radius of flyweights
    double fly_r_max;    // kg, Maximum radius of flyweights
};

// CVT parameters which are easily changed in reality
// that influence vehicle performance
struct CVTTune {
    double (*p_ramp_fn)(double x); // Function for ramp height vs. x, where x and height are 0 at the beginning of shift, all units in meters.
    double p_spring_k; // N/m
    double p_fly_mass_total; // kg
    double s_spring_k; // N/m
    double s_torsion_k; // N-m/rad
    double s_pretension; // radians
    double s_helix_angle; // radians
};

// Continuously Variable Transmission information
struct CVT {
    CVTConfig config;
    CVTTune tune;
    // realtime data
    double w_p; // rad/s of primary
    double w_s; // rad/s of secondary
};

struct Vehicle {
    VehicleConfig config;
    VehicleControls controls;
    CVT cvt;
    
    double w_wheel() {
        return cvt.w_s/config.fixed_gear_ratio;
    }
};

constexpr size_t sizeof_Vehicle = sizeof(Vehicle);

// Column units: rad/s, N-m
const Eigen::MatrixX2d CH440_TORQUE_CURVE = Eigen::MatrixX2d({
    {0             , 0.00},
    {1800*RPM2RADPS, 25.4},
    {2000*RPM2RADPS, 25.8},
    {2200*RPM2RADPS, 26.4},
    {2400*RPM2RADPS, 26.8},
    {2600*RPM2RADPS, 26.7},
    {2800*RPM2RADPS, 26.6},
    {3000*RPM2RADPS, 26.4},
    {3200*RPM2RADPS, 25.7},
    {3400*RPM2RADPS, 24.8},
    {3600*RPM2RADPS, 23.9},
});

constexpr CVTConfig GAGED_GX9_CONFIG = {
    .r_p = 0.975*IN2M,          // m, from WVU
    .R_p = 2.775*IN2M,          // m, from WVU
    .r_s = 2.498*IN2M,          // m, from WVU
    .R_s = 3.803*IN2M,          // m, from WVU
    .v_angle = 23*DEG2RAD,      // m, from WVU
    .p_gap = 0.8*IN2M,          // m, from WVU
    .s_gap = 0.8*IN2M,          // m, from WVU
    .L = 9.1*IN2M,              // m, as designed
    .l_b = 28*IN2M,             // m, measured
    .A_b = 0.41*IN2M*IN2M,      // m^2, from WVU
    .m_b = 0.249,               // kg, measured
    .E_b = 95147650.646,        // Pa, =13.8 ksi, from WVU
    .G_b = 35025367.049,        // Pa, =5.08 ksi, from WVU
    .mu_b = 0.13,               // from WVU
    .p_spring_x0 = 0.975*IN2M,  // m, from WVU
    .s_spring_x0 = 1.44*IN2M,   // m, from WVU
    .helix_radius = 1.455*IN2M, // m, from WVU
    .fly_arm_mass = 0.058375612,// kg, from WVU
    .fly_r_min = 1.875*IN2M,    // m, from WVU
    .fly_r_max = 2.345*IN2M,    // m, from WVU
};

// TODO
const Eigen::MatrixX2d GAGED_RAMP_POINTS = Eigen::MatrixX2d({
    {0, 0}
});

} // namespace TRB