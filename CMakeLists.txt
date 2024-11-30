# Minimum CMake version
cmake_minimum_required(VERSION 3.15)

# Project name and C++ standard
project(scf LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Include directories
include_directories(${CMAKE_SOURCE_DIR}/include)

# Source files
file(GLOB SRC_FILES ${CMAKE_SOURCE_DIR}/src/*.cpp)

# Add the executable
add_executable(${PROJECT_NAME} ${SRC_FILES})

# Find LAPACK
find_package(LAPACK REQUIRED)

# Find Boost Math
find_package(Boost REQUIRED COMPONENTS math_c99)

# Link LAPACK
if (LAPACK_FOUND)
  target_link_libraries(${PROJECT_NAME} PRIVATE LAPACK::LAPACK)
else()
    message(FATAL_ERROR "LAPACK not found!")
endif()

# Link Boost Math
if (Boost_FOUND)
    target_include_directories(${PROJECT_NAME} PRIVATE ${Boost_INCLUDE_DIRS})
    target_link_libraries(${PROJECT_NAME} PRIVATE ${Boost_LIBRARIES})
else()
    message(FATAL_ERROR "Boost Math not found!")
endif()

# Add extra warnings (optional but recommended)
if (CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(${PROJECT_NAME} PRIVATE -Wall -Wextra -pedantic -Werror)
elseif (MSVC)
    target_compile_options(${PROJECT_NAME} PRIVATE /W4 /WX)
endif()
