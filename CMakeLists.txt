cmake_minimum_required(VERSION 2.8)

# Project name
project(mc-mini)

# Path to installation directory
set(ROOT_DIR .)
 
# Compiler flags
SET(CMAKE_CXX_FLAGS "-std=c++11 -Wall -Wextra -Weffc++ -O3 -funroll-all-loops") 

# Path to hdf5 library 
find_library(hdf5 libhdf5_serial.so.8 PATHS /usr/lib/x86_64-linux-gnu)
find_library(hdf5_cpp libhdf5_cpp.so.8.0.2 PATHS /usr/lib/x86_64-linux-gnu)

# Include requred header files
SET(INCLUDE_DIR /usr/include)
include_directories(${INCLUDE_DIR}/hdf5/serial ${ROOT_DIR}/include) 

# Location of source files
SET(SOURCE_DIR ${ROOT_DIR}/source)
file(GLOB SOURCES ${SOURCE_DIR}/*.cpp ${SOURCE_DIR}/*/*.cpp)

# Name the executable file that will run "main.cpp"
add_executable(mc-mini ${SOURCES})

# Link the executable with the required libraries
target_link_libraries(mc-mini ${hdf5} ${hdf5_cpp})

# Build type; "Release" means enable all optimizaitons
set(CMAKE_BUILD_TYPE Release)