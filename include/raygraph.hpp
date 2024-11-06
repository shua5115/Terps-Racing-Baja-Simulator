#include "raylib.h"
#include "Eigen/Dense"

template<int N>
struct PlotConfig {
    const char *title;
    std::array<const char *, N> axis_labels;
    float fontsize;
    float fontspacing;
    Font font; // function which returns the font to use
    
    PlotConfig() : title(""), fontsize(24.0f), fontspacing(1.0f), font(GetFontDefault()) {
        axis_labels.fill("");
    }

    PlotConfig(const char *Title, std::array<const char *, N> Labels, float FontSize = 24.0f, float FontSpacing = 1.0f, Font font = {0})
        : title(Title), axis_labels(Labels), fontsize(FontSize), fontspacing(FontSpacing) {
        this->font = IsFontReady(font) ? font : GetFontDefault();
    }
};

Shader ScatterShader();

std::function<void()> DrawScatterPlot3D(const Eigen::MatrixXd &x, const Eigen::MatrixXd &y, const Eigen::MatrixXd &f, Vector3 position, Vector3 size, Camera cam, PlotConfig<3> config);
