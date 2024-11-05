#include <tuple>
#include <thread>
#include "trb.hpp"
#include "raylib.h"
#include "raymath.h"
#include "orbitcam.hpp"
#include "raygraph.hpp"

std::vector<Eigen::MatrixXd> gen_data(unsigned int N) {
    std::vector<std::thread> threads;
    threads.reserve(N);
    std::vector<Eigen::MatrixXd> data;
    data.resize(8, Eigen::MatrixXd());
    for(auto &d : data) {
        d.resize(N, N);
    }
    for(unsigned int i = 0; i < N; i++) {
        threads.emplace_back([&data, N, i](){
            for(unsigned int j = 0; j < N; j++) {
                BajaState baja = TR24_GAGED_GX9;
                baja.d_p = baja.d_p_max/2;
                auto S_fly = solve_flyweight_position(baja.theta1, baja.theta2, baja.cvt_tune.p_ramp_fn, baja.L_arm, baja.r_roller,
                    baja.d_p, baja.x_ramp, baja.r_cage, baja.r_shoulder
                );
                baja.theta1 = S_fly.x(0);
                baja.theta2 = S_fly.x(1);
                baja.controls.throttle = 1;
                baja.omega_p = remap(i, 0, N-1, 1800, 3800)*RPM2RADPS;
                baja.tau_s = remap(j, 0, N-1, 0*18.5*LBF2N/FT2M, 1*18.5*LBF2N/FT2M);
                double S = solve_cvt_shift(baja);
                baja.set_ratio_from_d_p(S);
                double ratio = baja.r_s/baja.r_p;
                double F_sp = baja.F_sp();
                double F_flyarm = baja.F_flyarm();
                double F_ss = baja.F_ss();
                double F_helix = baja.F_helix();
                double alpha = baja.alpha();
                double beta = 2*PI-alpha;
                double F_p = (F_flyarm - F_sp);
                double F_s = (F_ss + F_helix);
                double eq1 = -F_p + (alpha/beta)*F_s;
                data[0](i, j) = baja.omega_p;
                data[1](i, j) = baja.tau_s;
                data[2](i, j) = baja.d_p;
                data[3](i, j) = baja.d_s;
                data[4](i, j) = ratio;
                data[5](i, j) = F_p;
                data[6](i, j) = F_s;
                data[7](i, j) = F_p - (alpha/beta)*F_s;
                // printf("%f, %f, %f, %f, %f, %f, %f, %f\n", baja.omega_p, baja.tau_s, v(0), v(1), v(2), F_p, F_s, eq1);
            }
        });
    }
    for(auto &run : threads) {
        if (run.joinable()) run.join();
    }
    return data;
}

int main() {
    Eigen::initParallel();

    float mouse_sensitivity = 0.005;
    float key_sensitivity = PI/12;

    Camera3D cam {
        .position={-10,10,-10},
        .target={0,0,0},
        .up={0,1,0},
        .fovy=20,
        .projection=CAMERA_ORTHOGRAPHIC
    };
    float pitch = PI*0.125, yaw = -PI*0.25, orbit_dist = 25;

    auto data = gen_data(128);

    int graph = 0;

    SetConfigFlags(FLAG_VSYNC_HINT | FLAG_WINDOW_RESIZABLE);

    InitWindow(720, 720, "CVT Graphs");
    
    while(!WindowShouldClose()) {
        float dt = GetFrameTime();
        BeginDrawing();
        ClearBackground(WHITE);

        // Orbital cam
        auto dmouse = IsMouseButtonDown(MOUSE_BUTTON_LEFT) ? GetMouseDelta() : Vector2Zero();
        Vector3 orbitcam_input;
        orbitcam_input.x = mouse_sensitivity*dmouse.x + key_sensitivity*(IsKeyPressed(KEY_RIGHT) + IsKeyPressedRepeat(KEY_RIGHT) - IsKeyPressed(KEY_LEFT) - IsKeyPressedRepeat(KEY_LEFT));
        orbitcam_input.y = mouse_sensitivity*dmouse.y + key_sensitivity*(IsKeyPressed(KEY_DOWN) + IsKeyPressedRepeat(KEY_DOWN) - IsKeyPressed(KEY_UP) - IsKeyPressedRepeat(KEY_UP));
        orbitcam_input.z = -GetMouseWheelMove();
        UpdateOrbitCamera(&cam, orbitcam_input, &yaw, &pitch, &orbit_dist, {1, 100}, true);
        
        for(int key = KEY_ONE; key <= KEY_NINE; key++) {
            if (IsKeyPressed(key)) graph = key - KEY_ONE;
        }

        BeginMode3D(cam);
        
        std::function<void()> draw_labels = [](){};
        auto config = PlotConfig<3>("", {"x", "y", "z"}, 40);
        switch(graph) {
            case 0:
                config.title = "ratio vs. w_p, tau_s";
                config.axis_labels = {"w_p", "tau_s", "ratio"};
                draw_labels = DrawScatterPlot3D(data[0], data[1], data[4], {0, 0, 0}, {10, 4, 10}, cam, config);
                break;
            case 1:
                config.title = "F_p vs. w_p, tau_s";
                config.axis_labels = {"w_p", "tau_s", "F_p"};
                draw_labels = DrawScatterPlot3D(data[0], data[1], data[5], {0, 0, 0}, {10, 4, 10}, cam, config);
                break;
            case 2:
                config.title = "F_s vs. w_p, tau_s";
                config.axis_labels = {"w_p", "tau_s", "F_s"};
                draw_labels = DrawScatterPlot3D(data[0], data[1], data[6], {0, 0, 0}, {10, 4, 10}, cam, config);
                break;
            case 3:
                config.title = "F_p vs. d_p, w_p";
                config.axis_labels = {"d_p", "w_p", "F_p"};
                draw_labels = DrawScatterPlot3D(data[2], data[0], data[5], {0, 0, 0}, {10, 4, 10}, cam, config);
                break;
            case 4:
                config.title = "F_s vs. d_s, tau_s";
                config.axis_labels = {"d_s", "tau_s", "F_s"};
                draw_labels = DrawScatterPlot3D(data[3], data[1], data[6], {0, 0, 0}, {10, 4, 10}, cam, config);
                break;
            case 5:
                config.title = "F_eq vs. w_p, tau_s";
                config.axis_labels = {"w_p", "tau_s", "F_eq"};
                draw_labels = DrawScatterPlot3D(data[0], data[1], data[7], {0, 0, 0}, {10, 4, 10}, cam, config);
                break;
            case 6:
                config.title = "d_p vs. w_p, tau_s";
                config.axis_labels = {"w_p", "tau_s", "d_p"};
                draw_labels = DrawScatterPlot3D(data[0], data[1], data[2], {0, 0, 0}, {10, 4, 10}, cam, config);
                break;
            case 7:
                config.title = "d_s vs. w_p, tau_s";
                config.axis_labels = {"w_p", "tau_s", "d_s"};
                draw_labels = DrawScatterPlot3D(data[0], data[1], data[3], {0, 0, 0}, {10, 4, 10}, cam, config);
                break;
            case 8:
                config.title = "F_eq vs. d_p, d_s";
                config.axis_labels = {"d_p", "d_s", "F_eq"};
                draw_labels = DrawScatterPlot3D(data[2], data[3], data[7], {0, 0, 0}, {10, 4, 10}, cam, config);
                break;
        }
        EndMode3D();
        
        draw_labels();

        EndDrawing();
    }

    CloseWindow();
}