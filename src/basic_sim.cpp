#include <iostream>
#include "trb.hpp"

int main() {
    BajaState baja = TR24_GAGED_GX9;
    // Set initial condition and cvt tune
    baja.controls.throttle = 1; // this one is required for the car to move forward

    // Default CVT tune
    baja.cvt_tune.k_p = 60*LBF2N/IN2M;
    baja.cvt_tune.m_fly = 0.536;
    baja.cvt_tune.k_s = 20*LBF2N/IN2M;
    baja.cvt_tune.kappa_s = 0.4644*LBF2N*IN2M*RAD2DEG;
    baja.cvt_tune.theta_s_0 = (0.5)*24*DEG2RAD;
    baja.cvt_tune.theta_helix = DEG2RAD*33;

    double dt = 1e-4; // simulation timestep
    
    double t = 0;
    // CSV header
    printf("t,x,v,omega_p,ratio\n");
    // while(t < 10) { // time condition
    while(baja.x < 150*FT2M) { // distance condition
        double ratio = baja.r_s/baja.r_p;
        // CSV row
        printf("%f,%f,%f,%f,%f\n",
            t, baja.x, baja.v, baja.omega_p, ratio
        );
        // Step simulation
        trb_sim_step(baja, dt);
        t += dt;
    }
}