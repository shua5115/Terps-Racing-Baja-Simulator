cmake_minimum_required(VERSION 3.20)
project(trb_sim)

add_subdirectory("raylib-5.0")

include_directories(./include)
include_directories(.)

set(CMAKE_CXX_STANDARD 20)

# if(WIN32)
#     add_compile_options(/bigobj)
# endif()

set(model_src
src/trb.cpp
src/util.cpp
)

add_executable(test
src/test.cpp
${model_src}
)

add_executable(shift_visualizer
src/shift_visualizer.cpp
src/orbitcam.cpp
src/raygraph.cpp
${model_src}
)
target_link_libraries(shift_visualizer raylib)