cmake_minimum_required(VERSION 3.12)
project(SchwarzchildSim)

set(CMAKE_CXX_STANDARD 17)
find_package(OpenMP REQUIRED)

add_executable(SchwarzchildSim main.cpp)
target_link_libraries(SchwarzchildSim OpenMP::OpenMP_CXX)
target_include_directories(SchwarzchildSim PRIVATE ${CMAKE_SOURCE_DIR}/libs/stb)
