cmake_minimum_required(VERSION 3.17)
project(lab4_ai)
set(CMAKE_CXX_STANDARD 17)

set(SRCS src/SigmoidalMLP.cpp src/SigmoidalNeuron.cpp src/config.cpp)
include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(lab4_ai main.cpp ${SRCS})
configure_file(config.cfg config.cfg COPYONLY)