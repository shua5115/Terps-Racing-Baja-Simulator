# Terps Racing Baja Simulator

This repo contains a powertrain simulation for the Terps Racing Baja vehicle.
The sim is intended to be a tool to predict vehicle performance when changing design parameters, such as CVT tune and gear ratio.

## Theory

All equations used for this project are derived in the document [CVT Derivation.md](CVT%20Derivation.md).
Free body diagrams, rigid body equations of motion, algebra, and calculus are used to derive the governing equations for the simulation.

## Applications

shift_visualizer: Graphs CVT shift curves with respect to hill angle

vehicle_graph: Graphs engine speed, vehicle speed, and position with respect to time

test: Unit tests for various functions of the simulator, ensuring accuracy and correctness

opt_accuracy: Optimizes simulation accuracy by adjusting "free" variables, and comparing to collected data

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

## Dependencies

This project uses these great open-source libraries. Check them out!

### [Eigen 3.4](https://eigen.tuxfamily.org)
Eigen is a linear algebra library.

### [Raylib 5.0](https://github.com/raysan5/raylib/)
Raylib is a library which simplifies making graphical applications in C/C++.