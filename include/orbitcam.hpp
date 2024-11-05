#pragma once
#include "raylib.h"

void UpdateOrbitCamera(Camera *cam, Vector3 input, float *yaw, float *pitch, float *orbit_dist, Vector2 zoom_bounds = {1, 100}, bool snapping = false);