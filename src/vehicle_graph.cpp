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
    data.resize(6, Eigen::MatrixXd());

    double dt = 0.02;
    double tf = 10;
    size_t M = floor(tf/dt) + 1;

    for(auto &d : data) {
        d.resize(N, M);
    }
    for(unsigned int i = 0; i < N; i++) {
        threads.emplace_back([&data, N, M, i, dt](){
            BajaState baja = TR24_GAGED_GX9;
            baja.controls.throttle = 1;
            baja.omega_p = 1800*RPM2RADPS;
            double t = 0;

            baja.theta_hill = remap(i, 0, N-1, 0, 45*DEG2RAD);

            for(unsigned int j = 0; j < M; j++) {
                t = j*dt;
                auto dyn = trb_sim_step(baja, dt);
                data[0](i, j) = t;
                data[1](i, j) = dyn.x;
                data[2](i, j) = dyn.v;
                data[3](i, j) = dyn.omega_p;
                data[4](i, j) = baja.r_s/baja.r_p;
                data[5](i, j) = baja.theta_hill;
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
                config.title = "omega_p vs. time, hill angle";
                config.axis_labels = {"time", "hill angle", "omega_p"};
                draw_labels = DrawScatterPlot3D(data[0], data[5], data[3], {0, 0, 0}, {10, 4, 10}, cam, config);
                break;
            case 1:
                config.title = "ratio vs. time, hill angle";
                config.axis_labels = {"time", "hill angle", "ratio"};
                draw_labels = DrawScatterPlot3D(data[0], data[5], data[4], {0, 0, 0}, {10, 4, 10}, cam, config);
                break;
        }
        EndMode3D();
        
        draw_labels();

        EndDrawing();
    }

    CloseWindow();
}