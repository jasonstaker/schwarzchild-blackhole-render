# tests/CMakeLists.txt
cmake_minimum_required(VERSION 3.14)
include(CTest)

# Where to put test binaries (so CTest can find them easily)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/$<CONFIG>)

# Helper to define a GoogleTest‐based test
function(add_gtest_test name src)
  add_executable(${name} ${src})
  target_link_libraries(${name}
    PRIVATE
      gtest_main
      OpenMP::OpenMP_CXX
  )
  target_include_directories(${name} PRIVATE
    ${CMAKE_SOURCE_DIR}/math/core
    ${CMAKE_SOURCE_DIR}/math/relativity
  )
  # Register with CTest
  add_test(NAME ${name}
           COMMAND $<TARGET_FILE:${name}>
           WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/$<CONFIG>)
endfunction()

# Declare your three test executables
add_gtest_test(test_relativity test_relativity.cpp)
add_gtest_test(test_vec4 test_vec4.cpp)