#include <tuple>
#include <thread>
#include "trb.hpp"
#include "raylib.h"
#include "raymath.h"

double cvt_ratio(double omega_p, double tau_s) {
    BajaState state = TR24_GAGED_GX9;
    state.controls.throttle = 1;
    state.d_p = state.d_p_max/2; // set so flyweight solver is in medium section

    auto S_fly = solve_flyweight_position(state.theta1, state.theta2,
        state.d_r, state.cvt_tune.p_ramp_fn, state.L_arm, state.r_roller,
        state.d_p, state.x_ramp, state.r_cage, state.r_shoulder
    );
    state.theta1 = S_fly.x(0);
    state.theta2 = S_fly.x(1);
    state.d_r = S_fly.x(2);
    
    state.omega_p = omega_p;
    state.tau_s = tau_s;
    double d_p = solve_cvt_shift(state);
    state.set_ratio_from_d_p(d_p);
    return state.r_s/state.r_p;
}

struct ShiftRatioData {
    Eigen::MatrixXd omega_p, tau_s, ratio;
    ShiftRatioData() : omega_p(), tau_s(), ratio() {}
};

ShiftRatioData gen_data(unsigned int N) {
    std::vector<std::thread> threads;
    threads.reserve(N);
    ShiftRatioData data;
    data.omega_p.resize(N, N);
    data.tau_s.resize(N, N);
    data.ratio.resize(N, N);
    for(unsigned int i = 0; i < N; i++) {
        threads.emplace_back([&data, N, i](){
            for(unsigned int j = 0; j < N; j++) {
                double omega_p = remap(i, 0, N-1, 1800, 3800)*RPM2RADPS;
                double tau_s = remap(j, 0, N-1, 0, 4*18.5*LBF2N/FT2M);
                double ratio = cvt_ratio(omega_p, tau_s);
                data.omega_p(i, j) = omega_p;
                data.tau_s(i, j) = tau_s;
                data.ratio(i, j) = ratio;
            }
        });
    }
    for(auto &run : threads) {
        if (run.joinable()) run.join();
    }
    return data;
}

// Draws a 3D bar plot.
// Axis labels are drawn using the returned lambda while in the default draw mode.
std::function<void()> DrawBarPlot3D(const Eigen::MatrixXd &x, const Eigen::MatrixXd &y, const Eigen::MatrixXd &f, Vector3 position, Vector3 size, Camera cam, float fontsize = 24, float font_spacing = 1, Font font = GetFontDefault()) {
    int N = f.rows();
    int M = f.cols();
    std::array<std::tuple<Vector2, double>, 6> label_data; // position, value
    auto [min_x, max_x] = minmax(x);
    auto [min_y, max_y] = minmax(f);
    auto [min_z, max_z] = minmax(y);
    float step_x = (max_x-min_x)/N;
    float step_z = (max_z-min_z)/N;
    Vector3 scale = {
        (float) (size.x/(max_x-min_x)),
        (float) (size.y/(max_y-min_y)),
        (float) (size.z/(max_z-min_z))
    };
    position = Vector3Subtract(position, Vector3Scale(size, 0.5f));
    // Axes
    auto x_axis_start = Vector3Add({size.x + step_x*scale.x, 0, -step_z*scale.z}, position);
    auto x_axis_end = Vector3Add({0, 0, -step_z*scale.z}, position);
    auto y_axis_start = Vector3Add({size.x + step_x*scale.x, 0, -step_z*scale.z}, position);
    auto y_axis_end = Vector3Add({size.x + step_x*scale.x, size.y, -step_z*scale.z}, position);
    auto z_axis_start = Vector3Add({size.x + step_x*scale.x, 0, -step_z*scale.z}, position);
    auto z_axis_end = Vector3Add({size.x + step_x*scale.x, 0, size.z}, position);
    DrawLine3D(x_axis_start, x_axis_end, RED);
    DrawLine3D(y_axis_start, y_axis_end, GREEN);
    DrawLine3D(z_axis_start, z_axis_end, BLUE);
    
    for(unsigned int i = 0; i < N; i++) {
        for(unsigned int j = 0; j < N; j++) {
            Vector3 pos, cubesize;
            Color col = WHITE;
            double val = f(i, j);
            pos.x = scale.x * (float) remap(i, N-1, 0, 0, max_x-min_x);
            pos.y = 0.5 * scale.y * (float) val; // /2 b/c box is centered, we want it top aligned
            pos.z = scale.z * (float) remap(j, 0, N-1, 0, max_z-min_z);
            col.r = (char) remap(val, min_y, max_y, 0, 255);
            col.g = (char) remap(val, min_y, max_y, 128, 128);
            col.b = (char) remap(val, min_y, max_y, 255, 0);
            cubesize.x = step_x*scale.x;
            cubesize.y = 2*pos.y;
            cubesize.z = step_z*scale.z;
            pos = Vector3Add(pos, position);
            DrawPlane(Vector3Add(pos, {0, 0.5f*cubesize.y+0.01f, 0}), {cubesize.x, cubesize.z}, col);
            col.r *= 0.75;
            col.g *= 0.75;
            col.b *= 0.75;
            DrawCube(pos, cubesize.x, cubesize.y, cubesize.z, col);
        }
    }

    label_data[0] = std::tuple(Vector2Add(GetWorldToScreen(x_axis_start, cam), {MeasureTextEx(font, "_", fontsize, font_spacing).x, 0}), min_x);
    label_data[1] = std::tuple(Vector2Add(GetWorldToScreen(x_axis_end, cam), {MeasureTextEx(font, "_", fontsize, font_spacing).x, 0}), max_x);
    label_data[2] = std::tuple(Vector2Add(GetWorldToScreen(y_axis_start, cam), {-MeasureTextEx(font, TextFormat("_%.2f", min_y), fontsize, font_spacing).x, 0}), min_y);
    label_data[3] = std::tuple(Vector2Add(GetWorldToScreen(y_axis_end, cam), {-MeasureTextEx(font, TextFormat("_%.2f", max_y), fontsize, font_spacing).x, 0}), max_y);
    label_data[4] = std::tuple(Vector2Add(GetWorldToScreen(z_axis_start, cam), {-MeasureTextEx(font, TextFormat("_%.2f", min_z), fontsize, font_spacing).x, 0}), min_z);
    label_data[5] = std::tuple(Vector2Add(GetWorldToScreen(z_axis_end, cam), {-MeasureTextEx(font, TextFormat("_%.2f", max_z), fontsize, font_spacing).x, 0}), max_z);

    auto callback = [label_data, font, fontsize, font_spacing](){
        for(auto &data : label_data) {
            const auto [pos, val] = data;
            DrawTextEx(font, TextFormat("%.2f", val), pos, fontsize, 1, BLACK);
        }
    };

    return callback;
}

int main() {
    Eigen::initParallel();

    float mouse_sensitivity = 0.01;
    float key_sensitivity = 1;

    Camera3D cam {
        .position={-10,10,-10},
        .target={0,0,0},
        .up={0,1,0},
        .fovy=40,
        .projection=CAMERA_ORTHOGRAPHIC
    };
    float pitch = PI*0.125, yaw = -PI*0.25, orbit_dist = 25;

    auto data = gen_data(128);

    SetConfigFlags(FLAG_VSYNC_HINT | FLAG_WINDOW_RESIZABLE);

    InitWindow(720, 720, "CVT Ratio vs. omega_p, tau_s");

    while(!WindowShouldClose()) {
        float dt = GetFrameTime();
        BeginDrawing();
        ClearBackground(WHITE);

        // Orbital cam
        auto dmouse = IsMouseButtonDown(MOUSE_BUTTON_LEFT) ? GetMouseDelta() : Vector2Zero();
        yaw = yaw + mouse_sensitivity*dmouse.x + key_sensitivity*(IsKeyDown(KEY_RIGHT) - IsKeyDown(KEY_LEFT))*dt;
        pitch = clamp(pitch + mouse_sensitivity*dmouse.y + key_sensitivity*(IsKeyDown(KEY_UP) - IsKeyDown(KEY_DOWN))*dt, -PI*0.5, PI*0.5);
        if (cam.projection == CAMERA_ORTHOGRAPHIC) {
            cam.fovy = clamp(cam.fovy - GetMouseWheelMove(), 1, 100);
        } else {
            orbit_dist = clamp(orbit_dist - GetMouseWheelMove(), 1, 100);
        }
        cam.position.x = orbit_dist*cosf(yaw)*cosf(pitch);
        cam.position.y = orbit_dist*sinf(pitch);
        cam.position.z = orbit_dist*sinf(yaw)*cosf(pitch);
        cam.up.x = -cosf(yaw)*sinf(pitch);
        cam.up.y = cosf(pitch);
        cam.up.z = -sinf(yaw)*sinf(pitch);
        BeginMode3D(cam);
        
        auto draw_axes = DrawBarPlot3D(data.omega_p, data.tau_s, data.ratio, {0, 0, 0}, {10, 4, 10}, cam);

        EndMode3D();
        
        draw_axes();

        EndDrawing();
    }

    CloseWindow();
}