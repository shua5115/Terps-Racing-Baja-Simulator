/*
 * Usage: opt_accuracy <tunes_folder> <tunes_csv_filename>
 * - tunes_folder: folder containing CSV data for accel runs
 * - tunes_csv_filename: cssv file within the tunes_folder with metadata about each of the runs (see the example file for layout)
 * tunes_folder by default = "../../data/tunes", which selects the tunes in this repository if building into a folder
 * - This works because CMake will place the executable into a folder "build/<build-type>/opt_accuracy.exe"
 *   Going back two folders reaches the root of the repository, where the data/tunes folders exist.
 * tunes_csv_filename by default = "tunes.csv"
 */

#include <iostream>
#include <vector>
#include <fstream>
#include <filesystem>
#include <string.h>
#include <direct.h>
#include <thread>
#include "trb.hpp"
#include "opt.hpp"
#include "util.hpp"
#include "csv.hpp"

#define sim_dt (1e-3)

int main(int argc, char **argv) {
    using namespace std::filesystem;

    path folderpath;
    path tunes_filename = "tunes.csv";
    {
        char folderpath_buf[4096] = {0};
        // getcwd(folderpath_buf, 4095);
        strcpy(folderpath_buf, "../../data/tunes/");
        if (argc > 1) {
            strncpy(folderpath_buf, argv[1], 4095);
        }
        folderpath = folderpath_buf;
        if (argc > 2) {
            tunes_filename = argv[2];
        }
    }

    path tunes_filepath = folderpath / tunes_filename;
    std::cout << "Reading tunes from: " << tunes_filepath << "\n";

    auto csv_file = std::ifstream(tunes_filepath);

    if (!csv_file) {
        std::cerr << "Could not open tunes file!\n";
        exit(1);
    }
    std::map<std::string, size_t> header;
    auto tunes = read_csv(csv_file, &header);
    csv_file.close();
    std::vector<std::string> index_to_header;
    index_to_header.resize(header.size());

    for(const auto &elems : header) {
        index_to_header[elems.second] = elems.first;
    }

    std::cout << "Got tunes:\n";
    for(const auto &name : index_to_header) {
        std::cout << name << " | ";
    }
    std::cout << "\n";
    for(const auto &row : tunes) {
        for(const auto &elem : row) {
            std::cout << elem << " | ";
        }
        std::cout << "\n";
    }
    size_t N_tunes = tunes.size();
    size_t logfile_index = header.at("logfile");
    size_t engine_rpm_header_index = header.at("engine_rpm_header");
    size_t wheel_rpm_header_index = header.at("wheel_rpm_header");
    size_t k_p_index = header.at("k_p");
    size_t m_fly_index = header.at("m_fly");
    size_t k_s_index = header.at("k_s");
    size_t kappa_s_index = header.at("kappa_s");
    size_t theta_helix_index = header.at("theta_helix");
    size_t pretension_index = header.at("pretension");
    size_t t_start_index = header.at("t_start");
    size_t t_end_index = header.at("t_end");

    // Load captured data and initialize each BajaStates's cvt tune
    std::vector<BajaState> states(N_tunes, TR24_GAGED_GX9);
    std::vector<double> initial_engine_rpms(N_tunes, 0.0);
    std::vector<double> initial_wheel_rpms(N_tunes, 0.0);
    std::vector<double> initial_d_p(N_tunes, 0.0);
    std::vector<double> t_start_ms(N_tunes, 0.0);
    std::vector<double> t_end_ms(N_tunes, 0.0);
    std::vector<std::vector<std::array<double, 3>>> datasets;
    datasets.reserve(N_tunes);

    for (size_t i = 0; i < N_tunes; i++) {
        const auto &row = tunes.at(i);
        auto &state = states.at(i);

        state.controls.throttle = 1.0;
        state.cvt_tune.k_p = atof(row.at(k_p_index).c_str())*(LBF2N/IN2M); // lb/in to N/m
        state.cvt_tune.m_fly = atof(row.at(m_fly_index).c_str())*1e-3; // g to kg
        state.cvt_tune.k_s = atof(row.at(k_s_index).c_str())*(LBF2N/IN2M); // lb/in to N/m
        state.cvt_tune.kappa_s = atof(row.at(kappa_s_index).c_str())*(LBF2N*IN2M*RAD2DEG); // lb-in/deg to N-m/rad
        state.cvt_tune.theta_helix = atof(row.at(theta_helix_index).c_str())*DEG2RAD;
        state.cvt_tune.theta_s_0 = pretension_hole_to_theta_0_s((int) atof(row.at(pretension_index).c_str()));
        t_start_ms[i] = atof(row.at(t_start_index).c_str());
        t_end_ms[i] = atof(row.at(t_end_index).c_str());

        path datapath = folderpath / row.at(logfile_index);
        std::ifstream logfile(datapath);
        if (!logfile) {
            std::cerr << "Could not open rc log file " << datapath << "\n";
            exit(1);
        }
        // Adding a dataset with columns: time, engine_rpm, wheel_rpm
        const std::array<const char *, 3> colnames = {"Interval", row.at(engine_rpm_header_index).c_str(), row.at(wheel_rpm_header_index).c_str()};
        auto &dataset = datasets.emplace_back(read_rc_log(logfile, colnames));
        logfile.close();

        // find first row where the time is greater than or equal to start time
        auto start_row_it = std::find_if(dataset.begin(), dataset.end(), [&t_start_ms, i](const std::array<double, 3> &elem){ return elem[0] >= t_start_ms[i]; });

        // set the initial conditions of the state based on collected data
        if (start_row_it != dataset.end()) {
            auto &start_row = *start_row_it;
            initial_engine_rpms[i] = start_row.at(1);
            initial_wheel_rpms[i] = start_row.at(2);
            double omega_p = initial_engine_rpms[i]*RPM2RADPS;
            double omega_s = initial_wheel_rpms[i]*RPM2RADPS*state.N_g;
            double target_ratio = omega_p/omega_s;
            // this root bisection calls another root bisection, it's biception
            initial_d_p[i] = root_bisection([&](double d_p){
                state.set_ratio_from_d_p(d_p);
                return state.r_s/state.r_p - target_ratio;
            }, 0, state.d_p_max, 20);
        } else {
            std::cerr << "File " << row.at(logfile_index) << " does not contain data at its indicated start time of " << t_start_ms[i] << " ms." << "\n";
            exit(1);
        }

        std::cout << "Reading " << datapath.filename() << ": " << dataset.size() << " data points, initial engine rpm is " << initial_engine_rpms[i] 
            << ", initial wheel rpm is " << initial_wheel_rpms[i] 
            << ", initial d_p is " << initial_d_p[i] << "\n";

        // Filter dataset for bad values
        while(start_row_it != dataset.end()) {
            auto &row = *start_row_it;
            if (row.at(0) > t_end_ms[i]) break;
            // If the engine rpm is below idle, this is a sensor error and we should interpolate using the average of the adjacent values instead.
            if (row.at(1) < state.rpm_idle && i > 0 && i < (dataset.size()-1)) {
                row[1] = 0.5*(dataset.at(i-1).at(1) + dataset.at(i+1).at(1));
            }
            start_row_it++;
        }
    }

    Eigen::Vector2d lower_bound(0.0, 0.0);
    Eigen::Vector2d upper_bound(100.0, 20.0);
    Eigen::Vector2d x0(TR24_GAGED_GX9.F_resist, TR24_GAGED_GX9.shift_speed);
    const char *name1 = "F_resist";
    const char *name2 = "shift_speed";

    // The objective function measures the difference between the measured data and sim data.
    // The optimizer will attempt to reduce the difference as much as possible.
    auto objective = [&](Eigen::Vector2d x){
        std::vector<std::thread> threads;
        // Threads will modify a single entry in this vector,
        // so it must be pre-allocated to avoid invalidating references
        std::vector<double> error_vals(N_tunes, 0.0);
        // Create a new thread to process each dataset
        for(size_t i = 0; i < N_tunes; i++) {
            auto &baja = states.at(i);
            // Adjust state tune based on variables
            baja.F_resist = clamp(x(0), lower_bound(0), upper_bound(0));
            baja.m_car = clamp(x(1), lower_bound(1), upper_bound(1));
            baja.shift_speed = 0.5;
            // Initialize sim state from data
            baja.omega_p = initial_engine_rpms[i]*RPM2RADPS;
            baja.v = initial_wheel_rpms[i]*RPM2RADPS*baja.r_wheel;
            baja.set_ratio_from_d_p(initial_d_p[i]);
            // Isolate this thread's variables
            auto &dataset = datasets.at(i);
            auto &sim_error = error_vals.at(i); // this is a reference to the error val this thread will overwrite
            double t_start = t_start_ms[i];
            double t_end = t_end_ms[i];
            // double t_final = durations.at(i);
            // create a thread to simulate the vehicle for <duration> seconds
            threads.emplace_back([&baja, &sim_error, &dataset, t_start, t_end](){
                double sum_errors = 0;
                size_t r = 0;
                while (dataset.at(r).at(0) < t_start) {
                    r++;
                }
                while(dataset.at(r+1).at(0) < t_end) {
                    auto row = dataset.at(r);
                    auto row_next = dataset.at(r+1);
                    double row_dt = (row_next.at(0) - row.at(0))*1e-3;
                    size_t N_iters = (size_t) (row_dt/sim_dt) + 1;
                    double t = 0, t_prev = 0;
                    for(size_t iters = 0; iters < N_iters; iters++) {
                        t = std::min(t+sim_dt, row_dt);
                        double dt = t - t_prev; // the last timestep may not equal sim_dt
                        if (dt == 0) break;
                        auto step = trb_sim_step(baja, dt);
                        t_prev = t;
                    }
                    double real_omega_p = row_next.at(1)*RPM2RADPS;
                    double real_v = row_next.at(2)*RPM2RADPS*baja.r_wheel;
                    double err = abs(real_omega_p - baja.omega_p) + abs(real_v - baja.v);
                    err *= row_dt; // scale error by the timestep for timestep-independent error scales
                    sum_errors += err;
                    r++;
                }
                sim_error = sum_errors;
            });
        }
        // Wait for all threads to complete
        double error_grand_total = 0.0;
        for(size_t i = 0; i < N_tunes; i++) {
            auto &thread = threads.at(i);
            if (thread.joinable()) thread.join();
            error_grand_total += error_vals.at(i);
        }
        printf("%s=%.2f, %s=%.2e -> error=%.6e\n", name1, x(0), name2, x(1), error_grand_total);
        return error_grand_total;
    };

    auto opt_results = minimize_gradient_golden<2>(objective, x0, lower_bound, upper_bound, 1e-3, 0.1, 1.0, 1000);

    std::cout << "Solver " << (opt_results.converged ? "converged" : "did not converge") << " after " << opt_results.iterations << " iterations\n";
    std::cout << "Optimal configuration: " << opt_results.x << "\n";
    std::cout << "Resulting error: " << opt_results.f_of_x << "\n";

    path results_path = folderpath / "opt_results.csv";
    auto results_file = std::ofstream(results_path);

    if (results_file) {
        results_file << name1 << "," << opt_results.x(0) << "\n";
        results_file << name2 << "," << opt_results.x(1) << "\n";
        results_file << "error," << opt_results.f_of_x << "\n";
        results_file << "converged," << (opt_results.converged ? "True" : "False") << "\n";
        results_file << "iterations," << opt_results.iterations << "\n";
        results_file.close();
    }

    return EXIT_SUCCESS;
}