cmake_minimum_required(VERSION 3.10)

project(scf)

# Set C++14 standard
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Define source and include directories
set(SOURCE_DIR src)
set(INCLUDE_DIR include)

# Include directories
include_directories(${INCLUDE_DIR})

# Collect source files
file(GLOB SOURCES "${SOURCE_DIR}/*.cpp")

# Option to enable/disable DIIS
option(ENABLE_DIIS "Enable DIIS" ON)
if(NOT ENABLE_DIIS)
  add_definitions(-DNDIIS)
endif()

# Find required packages
find_package(Boost REQUIRED COMPONENTS math_c99)
find_package(LAPACKE REQUIRED)
find_package(LAPACK REQUIRED)

# Set output directories
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/build)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/build)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_OBJECT_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/build)

# Add executable
add_executable(scf ${SOURCES})

# Link libraries
target_link_libraries(scf
  Boost::math_c99
  ${LAPACK_LIBRARIES}
  ${LAPACKE_LIBRARIES}
)
