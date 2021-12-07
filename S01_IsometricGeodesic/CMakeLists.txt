cmake_minimum_required(VERSION 2.8)

project(S01_IsometricGeodesic)

find_package(Ceres REQUIRED PATHS $ENV{Ceres_DIR})
include_directories(${CERES_INCLUDE_DIRS})
include_directories("../nlohmann")

add_compile_definitions(_ENABLE_EXTENDED_ALIGNED_STORAGE)
include_directories($ENV{AC_DIR})
include_directories($ENV{MESHFRAME_DIRECTORY})

# Project
add_executable(S01_IsometricGeodesic main.cpp)
target_link_libraries(S01_IsometricGeodesic  ${CERES_LIBRARIES})
