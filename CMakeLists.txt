cmake_minimum_required(VERSION 3.10)

project(scf)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(SOURCE_DIR src)
set(INCLUDE_DIR include)

include_directories(${INCLUDE_DIR})

file(GLOB SOURCES "${SOURCE_DIR}/*.cpp")

option(ENABLE_DIIS "Enable DIIS" ON)
if(NOT ENABLE_DIIS)
  add_definitions(-DNDIIS)
endif()

find_package(Boost REQUIRED COMPONENTS math_c99)
find_package(LAPACK REQUIRED)
find_package(LAPACKE REQUIRED)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY
${CMAKE_BINARY_DIR}/build)

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY
${CMAKE_BINARY_DIR}/build)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY
${CMAKE_BINARY_DIR}/bin)

set(CMAKE_OBJECT_OUTPUT_DIRECTORY
${CMAKE_BINARY_DIR}/build)

add_executable(scf ${SOURCES})

target_link_libraries(scf Boost::math_c99
${LAPACK_LIBRARIES} ${LAPACKE_LIBRARIES})
