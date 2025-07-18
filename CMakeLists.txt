cmake_minimum_required(VERSION 3.12)
project(SchwarzschildSim)

set(CMAKE_CXX_STANDARD 17)
find_package(OpenMP REQUIRED)

add_executable(SchwarzschildSim main.cpp)
target_link_libraries(SchwarzschildSim OpenMP::OpenMP_CXX)
target_include_directories(SchwarzschildSim PRIVATE 
        ${CMAKE_SOURCE_DIR}/libs/stb/
        ${CMAKE_SOURCE_DIR}/math/core/
        ${CMAKE_SOURCE_DIR}/math/relativity/
        ${CMAKE_SOURCE_DIR}/renderer/)

enable_testing()
include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/refs/tags/v1.14.0.zip
)
# For Windows: prevents overriding our compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

add_subdirectory(tests)