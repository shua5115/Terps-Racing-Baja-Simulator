#include "util.hpp"

// 1-d linear interpolation lookup
// Matrix should have first column sorted in increasing order
double matrix_linear_lookup(Eigen::MatrixX2d mat, double val, bool extrapolate) {
    auto len = mat.rows();
    if (len < 1) {
        return 0;
    }
    if (len < 2) {
        return mat(0, 1);
    }
    if (val < mat(0, 0)) {
        if (extrapolate) {
            return remap(val, mat(0, 0), mat(1, 0), mat(0, 1), mat(1, 1));
        }
        return mat(0, 1);
    }
    for (size_t i = 0; i < len-1; i++) {
        double a = mat(i, 0);
        double b = mat(i+1, 0);
        if (val >= a && val <= b) {
            return remap(val, a, b, mat(i,1), mat(i+1,1));
        }
    }
    if (extrapolate) {
        return remap(val, mat(len-2, 0), mat(len-1, 0), mat(len-2, 1), mat(len-1, 1));
    }
    return mat(len-1, 1);
}

std::tuple<double, double> minmax(const Eigen::MatrixXd &m) {
    double min = INFINITY;
    double max = -INFINITY;
    for(Eigen::Index i = 0; i < m.rows(); i++) {
        for(Eigen::Index j = 0; j < m.cols(); j++) {
            double val = m(i, j);
            if (val < min) min = val;
            if (val > max) max = val;
        }
    }
    return std::tuple(min, max);
};