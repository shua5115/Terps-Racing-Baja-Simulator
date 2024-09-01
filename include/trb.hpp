#pragma once

#include <cmath>
#include <vector>
#include "Eigen/Dense"
#include "util.hpp"

namespace TRB {

struct VehicleConfig {
    // columns are: RPM, torque (N-m)
    Eigen::MatrixX2d engine_torque_curve;
    double fixed_gear_ratio;
    double tire_radius;
    double mass;
};

struct VehicleControls {
    double steering = 0;    // -1 to 1
    double throttle = 0;    // 0 to 1
    double brake_pedal = 0; // 0 to 1
    double awd_lever = 0;   // 0 or 1
};

// Variables about the CVT that are constant
struct CVTConfig {
    double r_p;          // m, inner radius of primary for belt
    double R_p;          // m, outer radius of primary for belt
    double r_s;          // m, inner radius of secondary for belt
    double R_s;          // m, outer radius of secondary fot belt
    double v_angle;      // radians, angle between sheaves
    double g_p;          // m, max gap between primary sheaves
    double g_s;          // m, max gap between secondary sheaves
    double L;            // m, distance between primary and secondary
    double l_b;          // m, belt length
    double A_b;          // m^2, belt cross section area
    double m_b;          // kg, mass of belt
    double E_b;          // Pa, Young's modulus of belt
    double G_b;          // Pa, shear modulus of belt
    double mu_b;         // Belt coefficient of friction
    double p_spring_x0;  // m, Primary spring initial displacement
    double s_spring_x0;  // m, Secondary spring initial displacement
    double r_helix;      // m, Radius of secondary helix ramp
    double m_fly_arm;    // kg, Mass of arm connected to flyweight
    // double r_fly_min;    // kg, Minimum radius of flyweights
    // double r_fly_max;    // kg, Maximum radius of flyweights

    // Dependent values
    
    // Effective coefficient of friction
    double mu_e() {
        return mu_b/std::sin(0.5*v_angle);
    }
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

// Continuously Variable Transmission real-time information
struct CVTState {
    double w_p; // rad/s of primary
    double w_s; // rad/s of secondary
    double r_p; // current radius of primary, based on sheave position
    double r_s; // current radius of secondary, based on sheave position

    // Belt wrap angle around the primary, in radians
    double belt_wrap_primary(double cvt_center_to_center_dist) {
        return 2.0*std::acos((r_s-r_p)/cvt_center_to_center_dist);
    }

    // Belt wrap angle around the secondary, in radians
    double belt_wrap_secondary(double cvt_center_to_center_dist) {
        return 2.0*std::acos((r_p-r_s)/cvt_center_to_center_dist);
    }

    // Current gear ratio of the CVT based on sheave position
    double shift_ratio() {
        return r_s/r_p;
    }

    double primary_clamp_force(CVTConfig &config, CVTTune &tune) {
        // TODO include 
        return 0;
    }

    // Equation 39 from WVU, removing *12 constant due to unecessary unit conversion (thanks metric!)
    double secondary_clamp_force(CVTConfig &config, CVTTune &tune, double plunge, double secondary_torque_load) {
        double shift_angle = plunge/(std::tan(tune.s_helix_angle)*config.r_helix);
        double from_linear = tune.s_spring_k*(config.s_spring_x0 + plunge);
        double from_torsion = tune.s_torsion_k*(tune.s_pretension + shift_angle) + secondary_torque_load/config.r_helix;
        from_torsion /= 2*std::tan(tune.s_helix_angle);
        return (from_linear + from_torsion)/belt_wrap_secondary(config.L);
    }

    // Equation 33 from WVU (literally equivalent to eq. 34)
    double secondary_slack_tension(double r, double w, double f_clamp, double wrap_angle, double v_angle, double belt_mass, double mu_e) {
        return (2*std::sin(wrap_angle*0.5)*(2*f_clamp*std::tan(v_angle*0.5) + belt_mass*r*r*w*w/12))/
            (std::cos(0.5*(wrap_angle-PI)) * (std::exp(mu_e*wrap_angle) + 1));
    }

    
};

struct Vehicle {
    VehicleConfig config;
    VehicleControls controls;
    CVTConfig cvt_config;
    CVTTune cvt_tune;
    CVTState cvt;

    double torque_from_vehicle_on_secondary(double slope_angle) {
        return config.mass*std::sin(slope_angle)*config.tire_radius/config.fixed_gear_ratio;
    }

    double torque_from_engine_on_primary(double engine_rpm) {
        return matrix_linear_lookup(config.engine_torque_curve, engine_rpm, false);
    }

    
};

constexpr size_t sizeof_Vehicle = sizeof(Vehicle);

// Column units: RPM, N-m
const Eigen::MatrixX2d CH440_TORQUE_CURVE = Eigen::MatrixX2d({
    {0   , 0},
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
    {3800, 0} // Governor prevents going above 3800, need to verify behavior
});

const VehicleConfig TR24_CONFIG = VehicleConfig{
    .engine_torque_curve=CH440_TORQUE_CURVE,
    .fixed_gear_ratio=8.32,
    .tire_radius=11.5*IN2M,
    .mass=500*LBF2KG,
};

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

} // namespace TRB