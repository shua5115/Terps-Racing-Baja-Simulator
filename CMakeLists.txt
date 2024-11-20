cmake_minimum_required(VERSION 3.20)
project(trb_sim)

add_subdirectory("raylib-5.0")

include_directories(./include)
include_directories(.)

set(CMAKE_CXX_STANDARD 20)

# if(WIN32)
#     add_compile_options(/bigobj)
# endif()

# These files are the minimum necessary to run the simulator
set(model_src
src/trb.cpp
src/util.cpp
)

# Runs a single CVT tune and prints out results as CSV
add_executable(basic_sim
src/basic_sim.cpp
${model_src}
)

# Unit tests for various functions of the simulator
add_executable(test
src/test.cpp
src/csv.cpp
${model_src}
)

# Testing reading CSV data
add_executable(test_csv
src/test_csv.cpp
src/csv.cpp
)

# Optimizes accuracy of the sim by comparing to real world data
add_executable(opt_accuracy
src/opt_accuracy.cpp
src/csv.cpp
${model_src}
)

# Optimizes cvt tune of vehicle for accel and hill climb events
add_executable(opt_cvt_tune
src/opt_cvt_tune.cpp
${model_src}
)

# Optimizes gear ratio assuming ideal CVT operation
add_executable(opt_gear_ratio
src/opt_gear_ratio.cpp
${model_src}
)

# Plots CVT shift behavior vs vehicle conditions
add_executable(shift_visualizer
src/shift_visualizer.cpp
src/orbitcam.cpp
src/raygraph.cpp
${model_src}
)
target_link_libraries(shift_visualizer raylib)

# Plots vehicle performance under various conditions
add_executable(vehicle_graph
src/vehicle_graph.cpp
src/orbitcam.cpp
src/raygraph.cpp
${model_src}
)
target_link_libraries(vehicle_graph raylib)