project(screwrobotics)
cmake_minimum_required(VERSION 3.10)
find_package(Eigen3 3.3 REQUIRED NO_MODULE)

add_definitions(-DTOL=1e-5)

include_directories(
    include
    ${EIGEN3_INCLUDE_DIR}
)

add_executable(main src/main.cpp)

target_link_libraries(main Eigen3::Eigen)

set_target_properties(main PROPERTIES CXX_STANDARD 17)