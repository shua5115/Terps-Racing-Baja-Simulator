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
    double C_d;         // Drag coefficient of the entire vehicle, when moving forward over flat ground
    double A_front;     // Frontal area of the vehicle, for calculating drag force
    double x;           // m, position of the vehicle in 1D space
    double v;           // m/s, velocity of the vehicle in 1D space
    double shift_speed; // 1/s, time constant for how fast the CVT shifts
    // Rear Brakes (assuming only rear brakes)
    double r_rotor;     // m, radius of the rear brake rotor
    double mu_brake;    // coefficient of friction for brake pads on rotor
    double max_brake_clamp; // N, maximum clamping force the brake caliper can apply to the rotor
    // Environment
    double g;           // m/s^2, acceleration of gravity
    double theta_hill;  // rad, angle of hill slope when driving
    double mu_ground;   // coefficient of friction between tires and ground surface
    double rho_air;     // kg/m^3, density of ambient air
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
    double mu_b;        // Belt static coefficient of friction with sheaves
    double I_e;         // kg-m^2, total moment of inertia of spinning engine components
    double I_p;         // kg-m^2, total moment of inertia of primary components
    double I_s;         // kg-m^2, total moment of inertia of secondary components
    double I_w;         // kg-m^2, total moment of inertia of all four wheels
    double F_resist;    // N, Friction force which opposes vehicle movement
    double F1;          // N, Friction force which opposes shifting
    // Primary
    double r_p_inner;   // m, radius where bottom of primary sheaves touch
    double d_p_max;     // m, max linear gap between primary sheaves
    double d_p_0;       // m, primary spring initial displacement
    double r_cage;      // m, radius from primary axis to outer edge of primary ramp
    double r_shoulder;  // m, radius from primary axis to flyweight arm pivot
    double L_arm;       // m, length of flyweight arm
    double r_roller;    // m, radius of rollers in primary
    double x_ramp;      // m, offset from flyweight arm pivot to furthest outwards edge of ramp, at low shift displacement
    // Secondary
    double r_s_inner;   // m, radius where bottom of secondary sheaves touch
    double d_s_max;     // m, max linear gap between secondary sheaves
    double d_s_0;       // m, Secondary spring initial displacement
    double r_helix;     // m, Radius of secondary helix ramp

    // Time-Dependent

    double tau_s;       // N-m, torque applied to secondary from gearbox
    double omega_p;     // rad/s, angular velocity of primary
    // Primary
    double r_p;         // m, radius of belt on primary sheave
    double d_p;         // m, linear displacement of primary sheave during shift
    double theta1;      // rad, angle beween flyweight arm and primary axis
    double theta2;      // rad, angle between primary ramp surface normal and primary axis
    // Secondary
    double r_s;         // m, radius of belt on secondary sheave
    double d_s;         // m, linear displacement of secondary sheave during shift

    // Misc

    void set_ratio_from_d_p(double d);
    
    double throttle_scale(double torque, double u_gas) const {
        return torque*(sin(u_gas*PI*0.5)*(1 - rpm_idle/rpm_gov) + rpm_idle/rpm_gov);
    }

    // Derived Constants

    // Average width of belt
    double b() const {
        return (b_min + b_max)/2; 
    }

    double r_p_min() const {
        return r_p_inner + 0.3825*h_v;
    }

    double r_s_min() const {
        return r_s_inner + 0.3825*h_v;
    }

    double r_p_max() const {
        return r_p_min() + d_p_max/tan(phi);
    }

    double r_s_max() const {
        return r_s_min() + d_s_max/tan(phi);
    }

    double r_fly() const {
        return r_shoulder + L_arm*sin(theta1);
    }

    double theta_s() const {
        return d_s/(r_helix*tan(cvt_tune.theta_helix));
    }

    double alpha() const {
        return 2*acos(clamp((r_s-r_p)/L, -1, 1));
    }

    double rho_b() const {
        return m_b/(A_b*L_b0);
    }

    double theta_s_max() const {
        return d_s_max/(r_helix*tan(cvt_tune.theta_helix));
    }

    // Forces and Moments

    // Effective mass of the vehicle (includes spinning mass of all components)
    double M_effective(bool belt_slipping = false) const;

    double calc_tau_s() const;

    double tau_e() const {
        double tau = matrix_linear_lookup(engine_torque_curve, RADPS2RPM*omega_p);
        return throttle_scale(tau, controls.throttle);
    }

    double F_sp() const {
        return cvt_tune.k_p*(d_p_0 + d_p);
    }

    double F_flyarm() const {
        return (cvt_tune.m_fly*r_fly()*omega_p*omega_p*L_arm*cos(theta1)*cos(theta2))/(L_arm*sin(theta1 + theta2) + r_roller*sin(2*theta2));
    }

    double F_ss() const {
        return cvt_tune.k_s*(d_s_0 + d_s);
    }

    double tau_ss() const {
        return cvt_tune.kappa_s*(cvt_tune.theta_s_0 + theta_s());
    }

    double F_helix() const {
        return (tau_s*0.5 + tau_ss())/(r_helix*tan(cvt_tune.theta_helix));
    }

    double F_wheel(double F_f) const {
        return F_f*r_s*N_g/r_wheel;
    }

    double F_f_max() const;

    // Sum of all resistance forces on the car:
    // - Rolling resistance
    // - Aerodynamic drag
    // - Applied brake force
    // Does NOT include force from gravity.
    double F_total_resist() const;
};

inline double pretension_hole_to_theta_0_s(int pretension) {
    return (0.5 + (pretension - 4))*24*DEG2RAD;
}

double primary_ramp_TR24(double x);

constexpr size_t BAJASTATE_SIZE_CHECK = sizeof(BajaState);

// BajaState representing the initial state of the TR24 baja car with an unmodified Gaged GX9 CVT
extern const BajaState TR24_GAGED_GX9;

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
    {3800, 0}, // Governor prevents torque increase above 3800, need to verify behavior
});

// Solves for theta1 and theta2 of the flyweight assembly.
OptResults<2> solve_flyweight_position(
    double theta1_guess, double theta2_guess,
    double (*ramp)(double x),
    double L1, double L2,
    double d_p, double x_ramp,
    double r_cage, double r_shoulder
);

// Solve secondary radius given primary radius using root finding.
// Per testing, numerical error is less than floating point precision after N=8 for this function for all reasonable values of r_p.
double solve_r_s(double r_p, double r_s_min, double r_s_max, double L, double L0, unsigned int N);

// Solves d_p, displacement of the primary, for the current state
double solve_cvt_shift(const BajaState &baja_state, int debuglevel = 0);

struct BajaDynamicsResult {
    double a; // Forward acceleration of the vehicle
    double v; // New forward velocity of the vehicle
    double x; // New position of the vehicle
    double omega_p; // New angular velocity of the primary
    double F_f; // Net belt tension between CVT sheaves, difference between slack and taut side tension
    bool slipping; // Whether the belt is slipping on the CVT primary
};

// Solves dynamics of the vehicle based on the current state with the following idealized assumptions:
// - Wheels cannot slip on the ground
BajaDynamicsResult solve_dynamics(const BajaState &baja_state, double dt);

// Changes values in state to reflect results from dynamics
void apply_dynamics_result(BajaState &state, const BajaDynamicsResult &dynamics);

// Updates the entire vehicle for one time step.
// Returns the dynamics result that it applies.
BajaDynamicsResult trb_sim_step(BajaState &state, double dt);