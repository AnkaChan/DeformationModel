cmake_minimum_required(VERSION 2.8)

project(S04_IsometricGeodesicOnlyZ)

find_package(Ceres REQUIRED PATHS $ENV{Ceres_DIR})
include_directories(${CERES_INCLUDE_DIRS})
include_directories("../nlohmann")

add_compile_definitions(_ENABLE_EXTENDED_ALIGNED_STORAGE)
include_directories($ENV{AC_DIR})
include_directories($ENV{MESHFRAME_DIRECTORY})

# Project
add_executable(S04_IsometricGeodesicOnlyZ main.cpp)
target_link_libraries(S04_IsometricGeodesicOnlyZ  ${CERES_LIBRARIES})
