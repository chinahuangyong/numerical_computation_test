cmake_minimum_required(VERSION 3.10)
project(nonlinear_pendulum)

set(CMAKE_CXX_STANDARD 17)

add_executable(nonlinear_pendulum main.cpp)
target_link_libraries(nonlinear_pendulum glog)
