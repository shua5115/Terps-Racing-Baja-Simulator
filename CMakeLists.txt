cmake_minimum_required(VERSION 3.20)
project(trb_sim)

include_directories(./include)
include_directories(.)

add_executable(baja_sim
src/main.cpp
src/util.cpp
src/opt.cpp
)
set_property(TARGET baja_sim PROPERTY CXX_STANDARD 20)

add_executable(test
src/test.cpp
src/util.cpp
src/opt.cpp
)
set_property(TARGET test PROPERTY CXX_STANDARD 20)