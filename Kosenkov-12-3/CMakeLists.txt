cmake_minimum_required(VERSION 3.17)
project(lab3_ai)
set(CMAKE_CXX_STANDARD 17)

set(SRCS src/EllipseDataGenerator.cpp src/SigmoidalNeuron.cpp src/config.cpp)
include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(lab3_ai main.cpp ${SRCS})
configure_file(config.cfg config.cfg COPYONLY)