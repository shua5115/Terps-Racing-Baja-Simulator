#include "orbitcam.hpp"
#include "raymath.h"
#include <cmath>

static double mod_euclid(double i, double n) {
    return fmod(n + fmod(i, n), n);
}

static float clamp(float v, float lo, float hi) {
    return std::min(std::max(v, lo), hi);
}

void UpdateOrbitCamera(Camera *cam, Vector3 input, float *yaw_ptr, float *pitch_ptr, float *orbit_dist_ptr, Vector2 zoom_bounds, bool snapping) {
    if (cam == NULL || yaw_ptr == NULL || pitch_ptr == NULL || orbit_dist_ptr == NULL) return;
    float yaw = *yaw_ptr;
    float pitch = *pitch_ptr;
    float orbit_dist = *orbit_dist_ptr;
    yaw = mod_euclid(yaw + input.x, 2*PI);
    pitch = clamp(pitch + input.y, -PI*0.5, PI*0.5);
    if (cam->projection == CAMERA_ORTHOGRAPHIC) {
        cam->fovy = clamp(cam->fovy + input.z, zoom_bounds.x, zoom_bounds.y);
    } else {
        orbit_dist = clamp(orbit_dist + input.z, zoom_bounds.x, zoom_bounds.y);
    }
    float yaw_adj = yaw;
    float pitch_adj = pitch;
    if (snapping) {
        float dist_from_45 = mod_euclid(yaw+PI*0.125, PI*0.25) - PI*0.125;
        if (fabs(dist_from_45) < 3*DEG2RAD) {
            yaw_adj -= dist_from_45;
        }
        if (fabs(pitch) < 3*DEG2RAD) {
            pitch_adj = 0;
        }
    }
    (*yaw_ptr) = yaw;
    (*pitch_ptr) = pitch;
    (*orbit_dist_ptr) = orbit_dist;

    cam->position.x = orbit_dist*cosf(yaw_adj)*cosf(pitch_adj);
    cam->position.y = orbit_dist*sinf(pitch_adj);
    cam->position.z = orbit_dist*sinf(yaw_adj)*cosf(pitch_adj);
    cam->up.x = -cosf(yaw_adj)*sinf(pitch_adj);
    cam->up.y = cosf(pitch_adj);
    cam->up.z = -sinf(yaw_adj)*sinf(pitch_adj);
}