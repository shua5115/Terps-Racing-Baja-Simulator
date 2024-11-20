#include <iostream>
#include <vector>
#include <array>
#include <thread>
#include "trb.hpp"
#include "opt.hpp"
#include "util.hpp"
#include "csv.hpp"

#define sim_dt (1e-3)

double accel_run_time(CVTTune tune, double hill_deg) {
    // Set up vehicle state
    BajaState baja = TR24_GAGED_GX9;
    baja.omega_p = RPM2RADPS*baja.rpm_idle;
    baja.theta_hill = hill_deg*DEG2RAD;
    baja.controls.throttle = 1;
    baja.cvt_tune = tune;

    double t = 0;
    while(baja.x < 150*FT2M) {
        t += sim_dt;
        trb_sim_step(baja, sim_dt);
    }
    return t;
}

struct AccelResult {
    CVTTune tune;
    double accel_time;
    double hill_time;
};

int main() {
    using V = Eigen::Vector<double, 5>;
    
    std::array<const char *, 3> k_p_names = {"black", "orange", "purple"};
    std::array<double, 3> k_p = {45, 60, 80}; // lbf/in

    // We have 11 x 35g, 5 x 56g, 12 x 67g weights
    constexpr std::array<double, 7> m_fly = {8*35, 4*35 + 4*56, 4*56 + 4*67, 503, 514, 525, 8*67};

    std::array<const char *, 2> k_s_names = {"yellow", "red"};
    std::array<double, 2> k_s = {22.75, 14.65}; // lbf/in
    std::array<double, 2> kappa_s = {0.4644, 0.592}; // lbf-in/deg

    std::array<double, 6> theta_helix = {28,32,33,36,38,40}; // deg
    std::array<double, 4> pretension = {4,5,6,7};

    constexpr size_t combos = k_p.size()*m_fly.size()*k_s.size()*theta_helix.size()*pretension.size();
    std::vector<AccelResult> results;
    results.resize(combos); // populate vector so we can use references to the items in indivual threads
    std::vector<std::thread> threads;

    size_t combo = 0;
    for(int i_k_p = 0; i_k_p < k_p.size(); i_k_p++) {
        for(int i_m_fly = 0; i_m_fly < m_fly.size(); i_m_fly++) {
            for(int i_k_s = 0; i_k_s < k_s.size(); i_k_s++) {
                threads.emplace_back([combo, i_k_p, i_m_fly, i_k_s, &k_p, &m_fly, &k_s, &kappa_s, &theta_helix, &pretension, &results](){
                    size_t subcombo = 0;
                    for(int i_pretension = 0; i_pretension < pretension.size(); i_pretension++) {
                        for(int i_helix = 0; i_helix < theta_helix.size(); i_helix++) {
                            CVTTune tune;
                            tune.p_ramp_fn = primary_ramp_TR24;
                            tune.k_p = k_p[i_k_p]*LBF2N/IN2M;
                            tune.m_fly = m_fly[i_m_fly]*1e-3; // g to kg
                            tune.k_s = k_s[i_k_s]*LBF2N/IN2M;
                            tune.kappa_s = kappa_s[i_k_s]*LBF2N*IN2M*RAD2DEG;
                            tune.theta_helix = theta_helix[i_helix]*DEG2RAD;
                            tune.theta_s_0 = (0.5 + (pretension[i_pretension]-4))*24*DEG2RAD;
                            printf("Starting run %llu/%llu\n", combo+subcombo, combos);
                            AccelResult &result = results.at(combo+subcombo);
                            result.tune = tune;
                            result.accel_time = accel_run_time(tune, 0);
                            result.hill_time = accel_run_time(tune, 30);
                            subcombo++;
                        }
                    }
                });
                combo += pretension.size()*theta_helix.size(); // must increment by product of all threaded for loop sizes
            }
        }
    }

    combo = 0;
    for(auto &thread : threads) {
        if(thread.joinable()) thread.join();
        printf("Finished run %llu/%llu\n", combo, combos);
        combo++;
    }

    constexpr double accel_importance = 10.0;
    constexpr double hill_importance = 1.0;
    auto objective = [=](AccelResult res){
        return res.accel_time*accel_importance + res.hill_time*hill_importance;
    };
    double best_val = INFINITY;
    double worst_val = 0;
    AccelResult best;
    AccelResult worst;
    for(const auto &result : results) {
        double f = objective(result);
        if(f < best_val) {
            best = result;
            best_val = f;
        }
        if(f > worst_val) {
            worst = result;
            worst_val = f;
        }
    }

    printf(
        "Worst result:\n"
        "k_p: %f lbf/in\n"
        "m_fly: %f g\n"
        "k_s: %f lbf/in\n"
        "helix: %f deg\n"
        "pretension: %f\n"
        "accel time: %f s\n"
        "hill time: %f s\n\n",
        worst.tune.k_p/(LBF2N/IN2M),
        worst.tune.m_fly*1e3,
        worst.tune.k_s/(LBF2N/IN2M),
        worst.tune.theta_helix*RAD2DEG,
        (worst.tune.theta_s_0/(24*DEG2RAD) - 0.5 + 4),
        worst.accel_time,
        worst.hill_time
    );

    printf(
        "Best result:\n"
        "k_p: %f lbf/in\n"
        "m_fly: %f g\n"
        "k_s: %f lbf/in\n"
        "helix: %f deg\n"
        "pretension: %f\n"
        "accel time: %f s\n"
        "hill time: %f s\n\n",
        best.tune.k_p/(LBF2N/IN2M),
        best.tune.m_fly*1e3,
        best.tune.k_s/(LBF2N/IN2M),
        best.tune.theta_helix*RAD2DEG,
        (best.tune.theta_s_0/(24*DEG2RAD) - 0.5 + 4),
        best.accel_time,
        best.hill_time
    );
}