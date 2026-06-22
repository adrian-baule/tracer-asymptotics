cmake_minimum_required(VERSION 3.14)
project(aoup_tracer VERSION 1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Release flags
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native -DNDEBUG")

# OpenMP
find_package(OpenMP REQUIRED)

# Main executable
add_executable(aoup_tracer
    src/main.cpp
    src/simulation.cpp
    src/aoup.cpp
    src/force.cpp
    src/sampling.cpp
)

target_include_directories(aoup_tracer PRIVATE include)
target_link_libraries(aoup_tracer PRIVATE OpenMP::OpenMP_CXX)

# Tests
enable_testing()

add_executable(test_force tests/test_force.cpp src/force.cpp)
target_include_directories(test_force PRIVATE include)
add_test(NAME test_force COMMAND test_force)

add_executable(test_aoup tests/test_aoup.cpp src/aoup.cpp)
target_include_directories(test_aoup PRIVATE include)
add_test(NAME test_aoup COMMAND test_aoup)

add_executable(test_sampling tests/test_sampling.cpp src/sampling.cpp)
target_include_directories(test_sampling PRIVATE include)
add_test(NAME test_sampling COMMAND test_sampling)
