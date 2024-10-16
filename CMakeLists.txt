cmake_minimum_required(VERSION 3.20)
project(trb_sim)

include_directories(./include)
include_directories(.)

set(CMAKE_CXX_STANDARD 20)

if(WIN32)
    add_compile_options(/bigobj)
endif()

set(model_src
src/util.cpp
src/opt.cpp
)

add_executable(test
src/test.cpp
${model_src}
)