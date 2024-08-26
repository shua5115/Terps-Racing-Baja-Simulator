#pragma once

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

// Variables about the CVT that do not change over time
struct CVTConfig {
    double r_p; // Meters, inner radius of primary
    double R_p; // Meters, outer radius of primary
    double r_s; // Meters, inner radius of secondary
    double R_s; // Meters, outer radius of secondary
    double v_angle; // radians, angle between sheaves
    double p_gap; // meters, max gap between primary sheaves
    double s_gap; // meters, max gap between secondary sheaves
    double l_b; // Meters, belt length
    double mu_b; // Belt coefficient of friction
};

struct CVTTune {
    double (*p_ramp_fn)(double x); // Function for ramp height vs. x, where x and height are 0 at the beginning of shift, all units in meters.
    double p_spring_k; // N/m
    double p_fly_mass_total; // kg
    double s_spring_k; // N/m, N-m/radian
    double s_pretension; // radians
    double s_helix_angle; // radians
};

struct CVT {
    CVTConfig config;
    CVTTune tune;
    double w_p; // rad/s of primary
    double w_s; // rad/s of secondary

    void update(std::vector<CVT> prev_state, double dt);
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
    {0             , 0},
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

// TODO
const Eigen::MatrixX2d GAGED_RAMP_POINTS = Eigen::MatrixX2d({
    {0, 0}
});

} // namespace TRB