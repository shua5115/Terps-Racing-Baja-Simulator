#include "raygraph.hpp"
#include "raymath.h"
#include "util.hpp"

// Draws a 3D bar plot.
// Axis labels are drawn using the returned lambda while in the default draw mode.
std::function<void()> DrawScatterPlot3D(const Eigen::MatrixXd &x, const Eigen::MatrixXd &y, const Eigen::MatrixXd &f, Vector3 position, Vector3 size, Camera cam, PlotConfig<3> config) {
    const char *title=config.title;
    const char *xlabel=config.axis_labels[0];
    const char *ylabel=config.axis_labels[1];
    const char *flabel=config.axis_labels[2];
    float fontsize = config.fontsize;
    float font_spacing = config.fontspacing;
    Font font = config.font;
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
    Vector2 title_pos = Vector2Add(center, {0,-0.5f*Vector2Distance(GetWorldToScreen(size, cam), GetWorldToScreen(Vector3Zero(), cam))});
    label_data[0] = std::tuple<Vector2, Vector2, double>(GetWorldToScreen(x_axis_start, cam), {MeasureTextEx(font, "_", fontsize, font_spacing).x, 0}, min_x);
    label_data[1] = std::tuple<Vector2, Vector2, double>(GetWorldToScreen(x_axis_end, cam), {MeasureTextEx(font, "_", fontsize, font_spacing).x, 0}, max_x);
    label_data[2] = std::tuple<Vector2, Vector2, double>(GetWorldToScreen(y_axis_start, cam), Vector2Multiply(MeasureTextEx(font, TextFormat("_%.2f", min_y), fontsize, font_spacing), {-0.5, -1}), min_y);
    label_data[3] = std::tuple<Vector2, Vector2, double>(GetWorldToScreen(y_axis_end, cam), Vector2Multiply(MeasureTextEx(font, TextFormat("_%.2f", max_y), fontsize, font_spacing), {-0.5, -1}), max_y);
    label_data[4] = std::tuple<Vector2, Vector2, double>(GetWorldToScreen(z_axis_start, cam), {-MeasureTextEx(font, TextFormat("_%.2f", min_z), fontsize, font_spacing).x, 0}, min_z);
    label_data[5] = std::tuple<Vector2, Vector2, double>(GetWorldToScreen(z_axis_end, cam), {-MeasureTextEx(font, TextFormat("_%.2f", max_z), fontsize, font_spacing).x, 0}, max_z);

    auto callback = [center, label_data, font, fontsize, font_spacing, title, title_pos, xlabel, ylabel, flabel](){
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
        Vector2 titlepos = title_pos; // redefine to remove const
        Vector2 titlesize = MeasureTextEx(font, title, fontsize, font_spacing);
        titlepos = Vector2Add(titlepos, Vector2Scale(titlesize, -0.5f));
        titlepos = {floorf(titlepos.x), floorf(titlepos.y)};
        DrawTextEx(font, title, titlepos, fontsize, font_spacing, BLACK);
    };

    return callback;
}