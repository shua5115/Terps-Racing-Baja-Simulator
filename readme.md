# Terps Racing Baja Simulator

This repo contains a powertrain simulation for the Terps Racing Baja vehicle.
The sim is intended to be a tool to predict vehicle performance when changing design parameters, such as CVT tune and gear ratio.

## Theory

All equations used for this project are derived in the document [CVT Derivation.md](CVT%20Derivation.md).
Free body diagrams, rigid body equations of motion, algebra, and calculus are used to derive the governing equations for the simulation.

## Building

There are no external code dependencies, so all you need is CMake and a C++ compiler.

### Windows

You can get a C++ compiler for windows by downloading [Visual Studio](https://visualstudio.microsoft.com/downloads/).
The Community Edition is free, and installing it will install the C++ compiler.

You can download CMake from [here](https://cmake.org/download/).

### Linux

For debian-based distributions, you can install these apt packages for a C++ compiler and CMake:
```bash
sudo apt install build-essential cmake
```

### Building with CMake

#### Using Command Prompt

Open a terminal in the terps-racing-baja-simulator directory and run the following commands:
```bash
mkdir build                                 # Create build directory
cd build                                    # 
cmake .. -DCMAKE_BUILD_TYPE=RelWithDebInfo  # Configure cmake
cmake --build . --config RelWithDebInfo     # Compile code
```
RelWithDebInfo is the preferred configuration because it allows the compiler to optimize math operations,
making the programs run much faster, while still allowing the program to be inspected with a debugger.

#### Using Visual Studio Code

Open this folder in VS code and install the C++ and CMake extensions from Microsoft.
With these extensions, you can build the project using the GUI.

## Applications

### basic_sim.cpp
Outputs simulation data in CSV format.

1. Set cvt tune variables in `basic_sim.cpp`
```c++
baja.cvt_tune.k_p; // N/m
baja.cvt_tune.m_fly; // kg
baja.cvt_tune.k_s; // N/m
baja.cvt_tune.kappa_s; // N-m/rad
baja.cvt_tune.theta_s_0; // rad
baja.cvt_tune.theta_helix; // rad
```
2. Build executables
```sh
# From root of repository
cd build
cmake .. -DCMAKE_BUILD_TYPE=RelWithDebInfo
cmake --build . --config RelWithDebInfo --target ALL_BUILD
```
3. Run `basic_sim` or `basic_sim.exe` and output to file `out.csv`
```sh
# From root of repository
cd build/RelWithDebInfo/
./basic_sim > out.csv
# or on windows
.\basic_sim.exe > out.csv
```

### test.cpp
Unit tests for various functions of the simulator, ensuring accuracy and correctness

### opt_accuracy.cpp
Optimizes simulation accuracy by adjusting "free" variables, and comparing to collected data.

Add racecapture CSV data to `data/tunes` and add row to `data/tunes.csv` to include file in optimization.

In `data/tunes.csv`, these headers determine what data is read from each racecapture file:
- logfile (string): Name of the file in the `data` directory
- engine_rpm_header (string): Name of the header in `logfile` for engine rpm
- wheel_rpm_header (string): Name of header in `logfile` for wheel rpm
- idle (RPM): Idle RPM of the engine during the run
- k_p (lbs/in): Primary spring constant
- m_fly (g): Flyweight mass in grams
- k_s (lbs/in): Secondary spring constant
- kappa_s (lb-in/deg): Secondary torsional spring constant
- theta_helix (deg): Angle of helix ramp in secondary
- pretension (integer): which hole in the secondary is the spring twisted to (counting from fully untwisted)?
- t_start (ms): Start of run in the `logfile`. Should be when the driver hits the throttle.
- t_end (ms): End of run in the `logfile`. Should be right before the car starts to slow down.

### shift_visualizer.cpp
Graphs CVT shift curves with respect to hill angle

### vehicle_graph.cpp
Graphs engine speed, vehicle speed, and position with respect to time

## Dependencies

This project uses these great open-source libraries. Check them out!

### [Eigen 3.4](https://eigen.tuxfamily.org)
Eigen is a linear algebra library.

### [Raylib 5.0](https://github.com/raysan5/raylib/)
Raylib is a library which simplifies making graphical applications in C/C++.