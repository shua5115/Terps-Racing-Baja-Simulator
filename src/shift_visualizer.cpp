#include <tuple>
#include <thread>
#include "trb.hpp"
#include "raylib.h"
#include "raymath.h"

void cvt_ratio(BajaState &state) {
    state.controls.throttle = 1;
    // state.d_p = state.d_p_max/2;
    // auto S_fly = solve_flyweight_position(state.theta1, state.theta2, state.cvt_tune.p_ramp_fn, state.L_arm, state.r_roller,
    //     state.d_p, state.x_ramp, state.r_cage, state.r_shoulder
    // );
    // state.theta1 = S_fly.x(0);
    // state.theta2 = S_fly.x(1);
    
}

struct ShiftRatioData {
    Eigen::MatrixXd omega_p, tau_s, ratio;
    ShiftRatioData() : omega_p(), tau_s(), ratio() {}
};

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

// Draws a 3D bar plot.
// Axis labels are drawn using the returned lambda while in the default draw mode.
std::function<void()> DrawScatterPlot3D(const Eigen::MatrixXd &x, const Eigen::MatrixXd &y, const Eigen::MatrixXd &f, Vector3 position, Vector3 size, Camera cam, const char *xlabel="x", const char *ylabel="y", const char *flabel="f", float fontsize = 24, float font_spacing = 1, Font font = GetFontDefault()) {
    int N = f.rows();
    int M = f.cols();
    std::array<std::tuple<Vector2, Vector2, double>, 6> label_data; // position, offset, value
    auto [min_x, max_x] = minmax(x);
    auto [min_y, max_y] = minmax(f);
    auto [min_z, max_z] = minmax(y);
    auto check_zero = [](double v, double alt){ return (v!=0) ? v : alt; };
    float step_x = check_zero(max_x-min_x, 1)/N;
    float step_y = 2*check_zero(max_y-min_y, 1)/(N+M);
    float step_z = check_zero(max_z-min_z, 1)/M;
    Vector3 scale = {
        (float) (size.x/check_zero(max_x-min_x, 1)),
        (float) (size.y/check_zero(max_y-min_y, 1)),
        (float) (size.z/check_zero(max_z-min_z, 1))
    };
    position = Vector3Subtract(position, Vector3Scale(size, 0.5f));
    // Axes
    float GAP = 4;
    auto x_axis_start = Vector3Add({size.x + GAP*step_x*scale.x, 0, -GAP*step_z*scale.z}, position);
    auto x_axis_end = Vector3Add({0, 0, -GAP*step_z*scale.z}, position);
    auto y_axis_start = Vector3Add({size.x + GAP*step_x*scale.x, 0, -GAP*step_z*scale.z}, position);
    auto y_axis_end = Vector3Add({size.x + GAP*step_x*scale.x, size.y, -GAP*step_z*scale.z}, position);
    auto z_axis_start = Vector3Add({size.x + GAP*step_x*scale.x, 0, -GAP*step_z*scale.z}, position);
    auto z_axis_end = Vector3Add({size.x + GAP*step_x*scale.x, 0, size.z}, position);
    DrawCubeV(Vector3Lerp(x_axis_start, x_axis_end, 0.5f), {Vector3Distance(x_axis_start, x_axis_end), step_z*scale.z, step_z*scale.z}, RED);
    DrawCubeV(Vector3Lerp(y_axis_start, y_axis_end, 0.5f), {step_x*scale.x, Vector3Distance(y_axis_start, y_axis_end), step_z*scale.z}, GREEN);
    DrawCubeV(Vector3Lerp(z_axis_start, z_axis_end, 0.5f), {step_x*scale.x, step_x*scale.x, Vector3Distance(x_axis_start, x_axis_end)}, BLUE);
    
    for(unsigned int i = 0; i < N; i++) {
        for(unsigned int j = 0; j < M; j++) {
            Vector3 pos, cubesize;
            Color col;
            double val = f(i, j);
            cubesize = {(float) check_zero(step_x*scale.x, size.x)*1.025f, 0, (float)check_zero(step_z*scale.z, size.z)*1.025f};
            cubesize.y = std::min(cubesize.x, cubesize.z);
            pos.x = (max_x==min_x) ? size.x : (float) remap(x(i, j), min_x, max_x, size.x, 0);
            pos.y = (max_y==min_y) ? 0 : (float) remap(val, min_y, max_y, 0, size.y);
            pos.z = (max_z==min_z) ? 0 : (float) remap(y(i, j), min_z, max_z, 0, size.z);
            col = ColorFromHSV(remap(val, min_y, max_y, 240, 0), 0.9, (abs(val) < 2*step_y) ? 0.0 : 1.0);
            pos = Vector3Add(pos, position);
            DrawCubeV(pos, cubesize, col);
        }
    }
    
    Vector2 center = GetWorldToScreen(Vector3Add(position, Vector3Scale(size, 0.5f)), cam);
    label_data[0] = std::tuple<Vector2, Vector2, double>(GetWorldToScreen(x_axis_start, cam), {MeasureTextEx(font, "_", fontsize, font_spacing).x, 0}, min_x);
    label_data[1] = std::tuple<Vector2, Vector2, double>(GetWorldToScreen(x_axis_end, cam), {MeasureTextEx(font, "_", fontsize, font_spacing).x, 0}, max_x);
    label_data[2] = std::tuple<Vector2, Vector2, double>(GetWorldToScreen(y_axis_start, cam), Vector2Multiply(MeasureTextEx(font, TextFormat("_%.2f", min_y), fontsize, font_spacing), {-0.5, -1}), min_y);
    label_data[3] = std::tuple<Vector2, Vector2, double>(GetWorldToScreen(y_axis_end, cam), Vector2Multiply(MeasureTextEx(font, TextFormat("_%.2f", max_y), fontsize, font_spacing), {-0.5, -1}), max_y);
    label_data[4] = std::tuple<Vector2, Vector2, double>(GetWorldToScreen(z_axis_start, cam), {-MeasureTextEx(font, TextFormat("_%.2f", min_z), fontsize, font_spacing).x, 0}, min_z);
    label_data[5] = std::tuple<Vector2, Vector2, double>(GetWorldToScreen(z_axis_end, cam), {-MeasureTextEx(font, TextFormat("_%.2f", max_z), fontsize, font_spacing).x, 0}, max_z);

    auto callback = [center, label_data, font, fontsize, font_spacing, xlabel, ylabel, flabel](){
        const char *labels[] = {xlabel, flabel, ylabel};
        for(size_t i = 0; i < label_data.size(); i++) {
            auto [p1, off1, val1] = label_data[(i/2)*2];
            auto [p2, off2, val2] = label_data[(i/2)*2 + 1];
            if (Vector2Equals(p1, p2)) continue;
            if (i%2 == 0) {
                const char *label = labels[i/2];
                Vector2 textsize = MeasureTextEx(font, label, fontsize, font_spacing);
                Vector2 labelpos = Vector2Lerp(p1, p2, 0.5f);
                if (i/2 != 1) {
                    labelpos = Vector2Subtract(labelpos, center);
                    labelpos = Vector2Add(labelpos, Vector2Scale(Vector2Normalize(labelpos), 0.5f*Vector2Length(textsize)));
                    labelpos = Vector2Add(labelpos, center);
                }
                labelpos = Vector2Add(labelpos, Vector2Scale(textsize, -0.5f));
                DrawTextEx(font, label, labelpos, fontsize, font_spacing, BLACK);
            }

            const auto [src, off, val] = label_data.at(i);
            Vector2 pos = Vector2Add(src, off);
            pos = {floorf(pos.x), floorf(pos.y)};
            DrawTextEx(font, TextFormat("%.2f", val), pos, fontsize, 1, BLACK);
        }
    };

    return callback;
}

int main() {
    Eigen::initParallel();

    float mouse_sensitivity = 0.005;
    float key_sensitivity = 1;

    Camera3D cam {
        .position={-10,10,-10},
        .target={0,0,0},
        .up={0,1,0},
        .fovy=20,
        .projection=CAMERA_ORTHOGRAPHIC
    };
    float pitch = PI*0.125, yaw = -PI*0.25, orbit_dist = 25;

    auto data = gen_data(256);

    int graph = 0;

    SetConfigFlags(FLAG_VSYNC_HINT | FLAG_WINDOW_RESIZABLE);

    InitWindow(720, 720, "CVT Ratio vs. omega_p, tau_s");
    
    while(!WindowShouldClose()) {
        float dt = GetFrameTime();
        BeginDrawing();
        ClearBackground(WHITE);

        // Orbital cam
        auto dmouse = IsMouseButtonDown(MOUSE_BUTTON_LEFT) ? GetMouseDelta() : Vector2Zero();
        yaw = yaw + mouse_sensitivity*dmouse.x + key_sensitivity*(IsKeyPressed(KEY_RIGHT) + IsKeyPressedRepeat(KEY_RIGHT) - IsKeyPressed(KEY_LEFT) - IsKeyPressedRepeat(KEY_LEFT))*dt;
        yaw = mod_euclid(yaw, 2*PI);
        pitch = clamp(pitch + mouse_sensitivity*dmouse.y + key_sensitivity*(IsKeyPressed(KEY_DOWN) + IsKeyPressedRepeat(KEY_DOWN) - IsKeyPressed(KEY_UP) - IsKeyPressedRepeat(KEY_UP))*dt, -PI*0.5, PI*0.5);
        if (cam.projection == CAMERA_ORTHOGRAPHIC) {
            cam.fovy = clamp(cam.fovy - GetMouseWheelMove(), 1, 100);
        } else {
            orbit_dist = clamp(orbit_dist - GetMouseWheelMove(), 1, 100);
        }
        double dist_from_45 = mod_euclid(yaw+PI*0.125, PI*0.25) - PI*0.125;
        double yaw_adj = yaw;
        double pitch_adj = pitch;
        if (abs(dist_from_45) < 3*DEG2RAD) {
            yaw_adj -= dist_from_45;
        }
        if (abs(pitch) < 3*DEG2RAD) {
            pitch_adj = 0;
        }

        cam.position.x = orbit_dist*cosf(yaw_adj)*cosf(pitch_adj);
        cam.position.y = orbit_dist*sinf(pitch_adj);
        cam.position.z = orbit_dist*sinf(yaw_adj)*cosf(pitch_adj);
        cam.up.x = -cosf(yaw_adj)*sinf(pitch_adj);
        cam.up.y = cosf(pitch_adj);
        cam.up.z = -sinf(yaw_adj)*sinf(pitch_adj);

        for(int key = KEY_ONE; key <= KEY_NINE; key++) {
            if (IsKeyPressed(key)) graph = key - KEY_ONE;
        }

        BeginMode3D(cam);
        
        std::function<void()> draw_axes = [](){};

        switch(graph) {
            case 0:
                draw_axes = DrawScatterPlot3D(data[0], data[1], data[4], {0, 0, 0}, {10, 4, 10}, cam, "w_p", "tau_s", "ratio", 40);
                break;
            case 1:
                draw_axes = DrawScatterPlot3D(data[0], data[1], data[5], {0, 0, 0}, {10, 4, 10}, cam, "w_p", "tau_s", "F_p", 40);
                break;
            case 2:
                draw_axes = DrawScatterPlot3D(data[0], data[1], data[6], {0, 0, 0}, {10, 4, 10}, cam, "w_p", "tau_s", "F_s", 40);
                break;
            case 3:
                draw_axes = DrawScatterPlot3D(data[2], data[0], data[5], {0, 0, 0}, {10, 4, 10}, cam, "d_p", "w_p", "F_p", 40);
                break;
            case 4:
                draw_axes = DrawScatterPlot3D(data[3], data[1], data[6], {0, 0, 0}, {10, 4, 10}, cam, "d_s", "tau_s", "F_s", 40);
                break;
            case 5:
                draw_axes = DrawScatterPlot3D(data[0], data[1], data[7], {0, 0, 0}, {10, 4, 10}, cam, "w_p", "tau_s", "F_eq", 40);
                break;
            case 6:
                draw_axes = DrawScatterPlot3D(data[0], data[1], data[2], {0, 0, 0}, {10, 4, 10}, cam, "w_p", "tau_s", "d_p", 40);
                break;
            case 7:
                draw_axes = DrawScatterPlot3D(data[0], data[1], data[3], {0, 0, 0}, {10, 4, 10}, cam, "w_p", "tau_s", "d_s", 40);
                break;
            case 8:
                draw_axes = DrawScatterPlot3D(data[2], data[3], data[7], {0, 0, 0}, {10, 4, 10}, cam, "d_p", "d_s", "F_eq", 40);
                break;
        }
        EndMode3D();
        
        draw_axes();

        EndDrawing();
    }

    CloseWindow();
}