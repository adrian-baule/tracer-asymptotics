cmake_minimum_required(VERSION 3.14)
project(aoup_tracer VERSION 1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Release flags
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native -DNDEBUG")

# OpenMP
find_package(OpenMP REQUIRED)

# One binary per swimmer process
foreach(PROC aoup abp rtp levy)
    string(TOUPPER ${PROC} PROC_UPPER)
    add_executable(tracer_${PROC}
        src/main.cpp
        src/simulation.cpp
        src/force.cpp
        src/sampling.cpp
    )
    target_include_directories(tracer_${PROC} PRIVATE include)
    target_compile_definitions(tracer_${PROC} PRIVATE PROCESS_${PROC_UPPER})
    target_link_libraries(tracer_${PROC} PRIVATE OpenMP::OpenMP_CXX)
endforeach()

# Tests
enable_testing()

add_executable(test_force tests/test_force.cpp src/force.cpp)
target_include_directories(test_force PRIVATE include)
add_test(NAME test_force COMMAND test_force)

add_executable(test_sampling tests/test_sampling.cpp src/sampling.cpp)
target_include_directories(test_sampling PRIVATE include)
add_test(NAME test_sampling COMMAND test_sampling)

add_executable(test_swimmers tests/test_swimmers.cpp src/force.cpp src/sampling.cpp)
target_include_directories(test_swimmers PRIVATE include)
target_compile_definitions(test_swimmers PRIVATE PROCESS_AOUP)
add_test(NAME test_swimmers COMMAND test_swimmers)
