# Set the minimum version of CMake that can be used
# to 3.8 to support C++17 features.
cmake_minimum_required(VERSION 3.22)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -ggdb3")

# This should point to the location of FindOpenVDB.cmake on your system
list(APPEND CMAKE_MODULE_PATH "/usr/local/lib/cmake/OpenVDB/")

# Set the project name
project(fluidsim)

find_package(GTest REQUIRED)
find_package(OpenVDB REQUIRED)
find_package(Boost COMPONENTS program_options REQUIRED)
find_package(Eigen3 REQUIRED)

# Define the executable
add_executable(${PROJECT_NAME} main.cpp grid.cpp)
target_link_libraries(${PROJECT_NAME} OpenVDB::openvdb ${Boost_LIBRARIES})
include_directories(${OpenVDB_INCLUDE_DIR})
include_directories(${EIGEN3_INCLUDE_DIR})
include_directories(${BOOST_INCLUDE_DIR})
include_directories(${GTEST_INCLUDE_DIRS})

# the tests
add_executable(tests test.cpp grid.cpp)
target_link_libraries(tests GTest::gtest GTest::gtest_main)
target_link_libraries(tests OpenVDB::openvdb)

# Set C++ standard to C++ 17
set_target_properties(${PROJECT_NAME} PROPERTIES
                      CXX_STANDARD 17
                      CXX_STANDARD_REQUIRED ON)
